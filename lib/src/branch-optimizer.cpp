/* Copyright (C) 2018-2021 INRAE
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "debug.hpp"
#include "itm.hpp"
#include "private.hpp"
#include "problem.hpp"
#include "utils.hpp"

#include <fstream>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <optional>
#include <set>
#include <utility>

/**
 * @brief Return true if lhs is better to rhs.
 */
static bool
is_best(const baryonyx::result& lhs,
        const baryonyx::result& rhs,
        baryonyx::objective_function_type type)
{
    if (lhs) {
        if (rhs) {
            if (type == baryonyx::objective_function_type::maximize) {
                return lhs.solutions.back().value >
                       rhs.solutions.back().value;
            } else {
                return lhs.solutions.back().value <
                       rhs.solutions.back().value;
            }
        } else {
            return true;
        }
    } else {
        if (rhs) {
            return false;
        } else {
            return lhs.remaining_constraints < rhs.remaining_constraints;
        }
    }
}

static bool
is_best(const baryonyx::context& ctx,
        const baryonyx::result& current,
        const baryonyx::result& best,
        baryonyx::objective_function_type type)
{
    if (is_best(current, best, type)) {
        if (current)
            baryonyx::notice(ctx,
                             "  - branch optimization found solution {:G}\n",
                             current.solutions.back().value);
        else
            baryonyx::notice(
              ctx,
              "  - branch optimization remaining_constraint {}\n",
              current.remaining_constraints);

        return true;
    }

    return false;
}

struct node
{
    node() = default;

    node(baryonyx::problem problem_,
         double solution_,
         int remaining_constraints_,
         int split_variable_)
      : problem(std::move(problem_))
      , solution(solution_)
      , remaining_constraints(remaining_constraints_)
      , split_variable(split_variable_)
    {}

    baryonyx::problem problem;
    double solution = 0.0;
    int remaining_constraints = -1;
    int split_variable = -1;
};

struct node_affected_compare
{
    bool operator()(const node& lhs, const node& rhs) const
    {
        return lhs.problem.affected_vars.values.size() >
               rhs.problem.affected_vars.values.size();
    }
};

struct node_result_compare
{
    bool operator()(const node& lhs, const node& rhs) const
    {
        if (lhs.remaining_constraints == 0) {
            if (rhs.remaining_constraints == 0) {
                return lhs.problem.type ==
                           baryonyx::objective_function_type::maximize
                         ? lhs.solution > rhs.solution
                         : lhs.solution < rhs.solution;
            }

            return true;
        } else {
            if (rhs.remaining_constraints == 0) {
                return false;
            } else {
                return lhs.remaining_constraints < rhs.remaining_constraints;
            }
        }
    }
};

using problem_set = std::multiset<node, node_result_compare>;

static bool
is_violated(const baryonyx::constraint& cst,
            const std::vector<std::int8_t>& x,
            int sense)
{
    int sum = 0;
    for (const auto& elem : cst.elements)
        sum += elem.factor * x[elem.variable_index];

    return sense == 0 ? sum != cst.value
                      : (sense > 0 ? sum < cst.value : sum > cst.value);
}

// Infeasible result: the variable shared by the most violated constraints,
// fixing it either way removes the conflict the solver keeps oscillating on.
// Feasible result: the selected variable that costs the most, its 0 branch
// looks for a cheaper solution. Ties go to the largest absolute cost.
static int
select_split_variable(const baryonyx::problem& pb, const baryonyx::result& r)
{
    const auto n = pb.vars.values.size();
    if (n == 0 || r.solutions.empty() ||
        r.solutions.front().variables.size() != n)
        return -1;

    const auto& x = r.solutions.front().variables;

    std::vector<double> cost(n, 0.0);
    for (const auto& elem : pb.objective.elements)
        cost[static_cast<std::size_t>(elem.variable_index)] += elem.factor;

    const double sign =
      pb.type == baryonyx::objective_function_type::maximize ? -1.0 : 1.0;

    std::vector<int> score(n, 0);
    bool feasible = true;

    const auto count = [&](const auto& constraints, int sense) {
        for (const auto& cst : constraints) {
            if (!is_violated(cst, x, sense))
                continue;

            feasible = false;
            for (const auto& elem : cst.elements)
                ++score[static_cast<std::size_t>(elem.variable_index)];
        }
    };

    count(pb.equal_constraints, 0);
    count(pb.greater_constraints, 1);
    count(pb.less_constraints, -1);

    int best = -1;
    double best_key = 0.0;
    double best_tie = 0.0;

    for (std::size_t i = 0; i != n; ++i) {
        double key, tie = std::abs(cost[i]);

        if (!feasible) {
            if (score[i] == 0)
                continue;
            key = score[i];
        } else {
            if (x[i] == 0)
                continue;
            key = sign * cost[i];
        }

        if (best < 0 || key > best_key ||
            (key == best_key && tie > best_tie)) {
            best = static_cast<int>(i);
            best_key = key;
            best_tie = tie;
        }
    }

    return best;
}

// Fixing a variable can make the sub-problem infeasible: that branch holds no
// solution and is pruned instead of aborting the search.
static std::optional<std::pair<baryonyx::problem, baryonyx::result>>
try_branch(const baryonyx::context& ctx,
           const baryonyx::problem& pb,
           int variable,
           bool value)
{
    try {
        auto sub = baryonyx::affect(ctx, pb, variable, value);
        auto ret = baryonyx::itm::optimize(ctx, sub);
        return std::make_pair(std::move(sub), std::move(ret));
    } catch (const baryonyx::solver_failure& e) {
        if (e.failure() != baryonyx::solver_error_tag::unrealisable_constraint)
            throw;
    }

    return std::nullopt;
}

static baryonyx::result
optimize(const baryonyx::context& ctx, const baryonyx::problem& pb)
{
    baryonyx::context internal_ctx(ctx);
    internal_ctx.log_priority = baryonyx::context::message_type::warning;
    internal_ctx.start = nullptr;
    internal_ctx.update = nullptr;
    internal_ctx.finish = nullptr;

    if (ctx.start)
        ctx.start(ctx.parameters);

    // The split tree can hold 2^n nodes: the time limit bounds the whole
    // search. The root gets half of it, a good starting point matters more
    // than many nodes, and each node then gets a sixteenth.
    const double budget = ctx.parameters.time_limit;
    const bool bounded = std::isfinite(budget) && budget > 0.0;
    const auto start = std::chrono::steady_clock::now();
    const auto remaining = [&]() {
        return budget - std::chrono::duration<double>(
                          std::chrono::steady_clock::now() - start)
                          .count();
    };
    if (bounded)
        internal_ctx.parameters.time_limit = budget / 2.0;

    auto best = baryonyx::itm::optimize(internal_ctx, pb);

    if (best)
        baryonyx::notice(internal_ctx,
                         "  - branch optimization found solution {:f}\n",
                         best.solutions.front().value);

    problem_set jobs;
    jobs.emplace(pb,
                 best ? best.solutions.back().value : 0,
                 best.remaining_constraints,
                 select_split_variable(pb, best));

    do {
        if (bounded) {
            const double left = remaining();
            if (left <= 0.0)
                break;
            internal_ctx.parameters.time_limit = std::min(budget / 16.0, left);
        }

        auto it = jobs.begin();
        const int variable = it->split_variable;

        if (variable < 0 ||
            variable >= baryonyx::length(it->problem.vars.names)) {
            jobs.erase(it);
            continue;
        }

        baryonyx::notice(internal_ctx,
                         "  - branch-optimization splits on {} (id: {})\n",
                         it->problem.vars.names[variable],
                         variable);

        for (bool value : { true, false }) {
            auto branch =
              try_branch(internal_ctx, it->problem, variable, value);
            if (!branch)
                continue;

            auto& [sub, ret] = *branch;

            if (is_best(internal_ctx, ret, best, it->problem.type))
                best = ret;

            jobs.emplace(std::move(sub),
                         ret ? ret.solutions.back().value : 0,
                         ret.remaining_constraints,
                         select_split_variable(sub, ret));
        }

        jobs.erase(it);

        baryonyx::notice(internal_ctx, "    Jobs lists:\n");
        for (const auto& elem : jobs) {
            if (elem.remaining_constraints == 0)
                baryonyx::notice(
                  internal_ctx,
                  "    * solution: {} split: {} vars: {}/{}\n",
                  elem.solution,
                  elem.split_variable,
                  elem.problem.affected_vars.values.size(),
                  elem.problem.affected_vars.values.size() +
                    elem.problem.vars.values.size());
            else
                baryonyx::notice(
                  internal_ctx,
                  "    * remaining: {} split: {} vars: {}/{}\n",
                  elem.remaining_constraints,
                  elem.split_variable,
                  elem.problem.affected_vars.values.size(),
                  elem.problem.affected_vars.values.size() +
                    elem.problem.vars.values.size());
        }

    } while (!jobs.empty() && !baryonyx::is_stop_requested());

    if (ctx.finish)
        ctx.finish(best);

    return best;
}

namespace baryonyx {
namespace itm {

result
branch_optimize(const baryonyx::context& ctx, const baryonyx::problem& pb)
{
    baryonyx::notice(ctx, "- branch-optimization starts\n");
    return ::optimize(ctx, pb);
}

} // namespace itm
} // namespace baryonyx
