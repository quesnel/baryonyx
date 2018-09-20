/* Copyright (C) 2018 INRA
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
#include <set>

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
                return lhs.solutions.front().value >
                       rhs.solutions.front().value;
            } else {
                return lhs.solutions.front().value <
                       rhs.solutions.front().value;
            }
        } else {
            return true;
        }
    } else {
        if (rhs) {
            return false;
        } else {
            return lhs.constraints < rhs.constraints;
        }
    }
}

static bool
is_best(const baryonyx::context_ptr& ctx,
        const baryonyx::result& current,
        const baryonyx::result& best,
        baryonyx::objective_function_type type)
{
    if (is_best(current, best, type)) {
        if (current)
            baryonyx::notice(ctx,
                             "  - branch optimization found solution {:G}\n",
                             current.solutions.front().value);
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
         int annoying_variable_)
      : problem(std::move(problem_))
      , solution(solution_)
      , remaining_constraints(remaining_constraints_)
      , annoying_variable(annoying_variable_)
    {}

    baryonyx::problem problem;
    double solution;
    int remaining_constraints;
    int annoying_variable;
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

using problem_set = std::set<node, node_result_compare>;

static baryonyx::result
optimize(const baryonyx::context_ptr& ctx, const baryonyx::problem& pb)
{
    // NOTE 1 - add a solver parameter to limit the size of the list.

    auto old_log_priority = ctx->log_priority;
#ifndef BARYONYX_ENABLE_DEBUG
    ctx->log_priority = baryonyx::context::message_type::notice;
#endif

    auto best = baryonyx::itm::optimize(ctx, pb);

    if (best)
        baryonyx::notice(ctx,
                         "  - branch optimization found solution {:f}\n",
                         best.solutions.front().value);

    problem_set jobs;
    jobs.emplace(pb,
                 best ? best.solutions.back().value : 0,
                 best.remaining_constraints,
                 best.annoying_variable);

    int annoying_variable = best.annoying_variable;

    do {
        auto it = jobs.begin();

        baryonyx::notice(ctx,
                         "  - branch-optimization splits on {} (id: {})\n",
                         it->problem.vars.names[it->annoying_variable],
                         it->annoying_variable);

        auto sp = baryonyx::split(ctx, it->problem, annoying_variable);

        auto ret0 = baryonyx::itm::optimize(ctx, std::get<0>(sp));
        if (is_best(ctx, ret0, best, it->problem.type))
            best = ret0;

        jobs.emplace(std::get<0>(sp),
                     ret0 ? ret0.solutions.back().value : 0,
                     ret0.remaining_constraints,
                     ret0.annoying_variable);

        auto ret1 = baryonyx::itm::optimize(ctx, std::get<1>(sp));
        if (is_best(ctx, ret1, best, it->problem.type))
            best = ret1;

        jobs.emplace(std::get<1>(sp),
                     ret1 ? ret1.solutions.back().value : 0,
                     ret1.remaining_constraints,
                     ret1.annoying_variable);

        jobs.erase(it);

        baryonyx::notice(ctx, "    Jobs lists:\n");
        for (const auto& elem : jobs) {
            if (elem.remaining_constraints == 0)
                baryonyx::notice(
                  ctx,
                  "    * solution: {} annoying: {} vars: {}/{}\n",
                  elem.solution,
                  elem.annoying_variable,
                  elem.problem.affected_vars.values.size(),
                  elem.problem.affected_vars.values.size() +
                    elem.problem.vars.values.size());
            else
                baryonyx::notice(
                  ctx,
                  "    * remaining: {} annoying: {} vars: {}/{}\n",
                  elem.remaining_constraints,
                  elem.annoying_variable,
                  elem.problem.affected_vars.values.size(),
                  elem.problem.affected_vars.values.size() +
                    elem.problem.vars.values.size());
        }

    } while (!jobs.empty());

    ctx->log_priority = old_log_priority;

    return best;
}

namespace baryonyx {
namespace itm {

result
branch_optimize(const baryonyx::context_ptr& ctx, const baryonyx::problem& pb)
{
    baryonyx::notice(ctx, "- branch-optimization starts\n");

    return ::optimize(ctx, pb);
}

} // namespace itm
} // namespace baryonyx
