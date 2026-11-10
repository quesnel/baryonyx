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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "debug.hpp"
#include "itm.hpp"
#include "private.hpp"
#include "problem.hpp"
#include "random-distributions.hpp"
#include "random-engine.hpp"
#include "utils.hpp"

namespace {

namespace bx = baryonyx;

using clock_type = std::chrono::steady_clock;

bool
is_better(double lhs, double rhs, bx::objective_function_type type) noexcept
{
    return type == bx::objective_function_type::maximize ? lhs > rhs
                                                         : lhs < rhs;
}

/** Constraints of every kind, as indices into a single numbering. */
class constraint_graph
{
public:
    explicit constraint_graph(const bx::problem& pb)
      : m_variable_constraints(pb.vars.values.size())
    {
        for (const auto* csts : { &pb.equal_constraints,
                                  &pb.greater_constraints,
                                  &pb.less_constraints }) {
            for (const auto& cst : *csts) {
                const int id = bx::length(m_constraint_variables);
                auto& vars = m_constraint_variables.emplace_back();

                for (const auto& elem : cst.elements) {
                    vars.emplace_back(elem.variable_index);
                    m_variable_constraints[static_cast<std::size_t>(
                                             elem.variable_index)]
                      .emplace_back(id);
                }
            }
        }
    }

    int constraints() const noexcept
    {
        return bx::length(m_constraint_variables);
    }

    // Breadth-first walk from a random constraint through the constraints
    // sharing a variable: the freed variables interact, so the sub-problem
    // can actually move.
    void walk(bx::xoshiro256pp& rng,
              int wanted,
              std::vector<char>& is_free,
              std::vector<char>& visited,
              std::vector<int>& queue) const
    {
        std::fill(visited.begin(), visited.end(), 0);
        queue.clear();

        int freed = 0;
        std::size_t head = 0;

        while (freed < wanted) {
            if (head == queue.size()) {
                const auto start = static_cast<int>(bx::uniform_below(
                  rng, static_cast<std::uint64_t>(constraints())));
                if (visited[static_cast<std::size_t>(start)])
                    continue;
                visited[static_cast<std::size_t>(start)] = 1;
                queue.emplace_back(start);
            }

            const auto cst = static_cast<std::size_t>(queue[head++]);
            for (const int var : m_constraint_variables[cst]) {
                const auto v = static_cast<std::size_t>(var);
                if (!is_free[v]) {
                    is_free[v] = 1;
                    ++freed;
                }

                for (const int next : m_variable_constraints[v]) {
                    const auto n = static_cast<std::size_t>(next);
                    if (!visited[n]) {
                        visited[n] = 1;
                        queue.emplace_back(next);
                    }
                }
            }
        }
    }

private:
    std::vector<std::vector<int>> m_constraint_variables;
    std::vector<std::vector<int>> m_variable_constraints;
};

/** Writes the values of a sub-problem result into the full assignment. */
void
merge(const bx::problem& sub_pb,
      const bx::result& sub,
      const std::unordered_map<std::string_view, int>& index,
      std::vector<std::int8_t>& x)
{
    const auto set = [&](std::string_view name, std::int8_t value) {
        if (const auto it = index.find(name); it != index.end())
            x[static_cast<std::size_t>(it->second)] = value;
    };

    for (std::size_t i = 0, e = sub.affected_vars.names.size(); i != e; ++i)
        set(sub.affected_vars.names[i], sub.affected_vars.values[i] ? 1 : 0);

    const auto& values = sub.solutions.back().variables;
    if (values.size() == sub.variable_name.size()) {
        for (std::size_t i = 0, e = sub.variable_name.size(); i != e; ++i)
            set(sub.variable_name[i], values[i]);
        return;
    }

    // Without constraint left the optimizers return no variable value: each
    // free variable takes the value of its linear cost.
    std::vector<double> cost(sub_pb.vars.names.size(), 0.0);
    for (const auto& elem : sub_pb.objective.elements)
        cost[static_cast<std::size_t>(elem.variable_index)] += elem.factor;

    const bool maximize =
      sub_pb.type == bx::objective_function_type::maximize;
    for (std::size_t i = 0, e = cost.size(); i != e; ++i)
        set(sub_pb.vars.names[i], (maximize ? cost[i] > 0 : cost[i] < 0));
}

// Large neighbourhood search: most variables keep their value in the best
// solution, a neighbourhood of k variables is freed and the reduced problem
// is solved by short restart runs. Wedelin solves such small problems in
// milliseconds where a full run takes seconds on large instances. Without
// time limit, the search runs until the process is stopped and reports each
// improvement through the update callback.
bx::result
optimize(const bx::context& ctx, const bx::problem& pb)
{
    const double budget = ctx.parameters.time_limit;
    const bool unbounded = !std::isfinite(budget) || budget <= 0.0;

    const auto start = clock_type::now();
    const auto elapsed = [&start]() {
        return std::chrono::duration<double>(clock_type::now() - start)
          .count();
    };

    bx::context internal(ctx);
    internal.start = nullptr;
    internal.update = nullptr;
    internal.finish = nullptr;
    internal.log_priority = bx::context::message_type::warning;

    if (ctx.start)
        ctx.start(ctx.parameters);

    using optimizer_type = bx::solver_parameters::optimizer_type;

    bx::result best;
    const auto run = [&](optimizer_type optimizer, double time_limit) {
        if (bx::is_stop_requested())
            return;

        internal.parameters = ctx.parameters;
        internal.parameters.optimizer = optimizer;
        internal.parameters.time_limit = time_limit;

        auto ret = bx::itm::optimize(internal, pb);
        if (ret || ret.remaining_constraints < best.remaining_constraints)
            best = std::move(ret);
    };

    // A first solution with restart, then with tuned parameters, which find
    // one on parameter-sensitive problems. With a time limit: a quarter, then
    // half of it; a candidate of the tuned optimizer must not tune again and
    // keeps searching with restart instead. Without time limit: alternate
    // both with doubling durations until a solution is found.
    if (unbounded) {
        for (double t = 1.0; !best && !bx::is_stop_requested(); t *= 2.0) {
            run(optimizer_type::restart, t);
            if (!best)
                run(optimizer_type::tuned, 2.0 * t);
        }
    } else {
        run(optimizer_type::restart, 0.25 * budget);
        if (!best) {
            if (ctx.tuning)
                run(optimizer_type::restart, 0.75 * budget);
            else
                run(optimizer_type::tuned, 0.5 * budget);
        }
    }

    internal.parameters = ctx.parameters;
    internal.parameters.optimizer = optimizer_type::restart;

    const auto n = pb.vars.values.size();
    if (!best || n == 0 || best.solutions.empty() ||
        best.solutions.back().variables.size() != n) {
        if (ctx.finish)
            ctx.finish(best);
        return best;
    }

    std::vector<std::int8_t> x = best.solutions.back().variables;
    double value = bx::compute_solution(pb, x);

    if (ctx.update)
        ctx.update(0, value, 0, elapsed(), 0);

    std::unordered_map<std::string_view, int> index;
    index.reserve(n);
    for (std::size_t i = 0; i != n; ++i)
        index.emplace(pb.vars.names[i], static_cast<int>(i));

    const constraint_graph graph(pb);
    bx::xoshiro256pp rng(
      ctx.parameters.seed > 0
        ? static_cast<bx::xoshiro256pp::result_type>(ctx.parameters.seed)
        : 1u);

    const double min_free = std::min<double>(20.0, static_cast<double>(n));
    const double max_free = std::max(min_free, 0.5 * static_cast<double>(n));
    double wanted =
      std::clamp(0.1 * static_cast<double>(n), min_free, max_free);

    // Short runs: most neighbourhoods are solved within a few milliseconds
    // and more neighbourhoods beat longer runs on each.
    internal.parameters.time_limit = 0.02;
    internal.parameters.thread = 2;

    std::vector<char> is_free(n);
    std::vector<char> visited(static_cast<std::size_t>(graph.constraints()));
    std::vector<int> queue;
    std::vector<std::pair<int, bool>> assignments;
    std::vector<std::int8_t> candidate;
    long int iterations = 0, improvements = 0;

    const auto has_time = [&]() {
        return unbounded ||
               elapsed() + internal.parameters.time_limit < budget;
    };

    while (!bx::is_stop_requested() && has_time()) {
        ++iterations;
        std::fill(is_free.begin(), is_free.end(), 0);
        const int k = static_cast<int>(wanted);

        if (iterations % 2 == 0 || graph.constraints() == 0) {
            for (int freed = 0; freed < k;) {
                const auto v = bx::uniform_below(rng, n);
                if (!is_free[v]) {
                    is_free[v] = 1;
                    ++freed;
                }
            }
        } else {
            graph.walk(rng, k, is_free, visited, queue);
        }

        assignments.clear();
        for (std::size_t i = 0; i != n; ++i)
            if (!is_free[i])
                assignments.emplace_back(static_cast<int>(i), x[i] != 0);

        bx::problem sub;
        try {
            sub = bx::affect(internal, pb, assignments);
        } catch (const bx::solver_failure&) {
            continue;
        }

        if (sub.vars.values.empty())
            continue;

        const auto ret = bx::itm::optimize(internal, sub);
        if (!ret || ret.solutions.empty()) {
            wanted = std::max(min_free, 0.8 * wanted);
            continue;
        }

        candidate = x;
        merge(sub, ret, index, candidate);
        const double candidate_value = bx::compute_solution(pb, candidate);

        if (!is_better(value, candidate_value, pb.type) &&
            bx::is_valid_solution(pb, candidate)) {
            if (is_better(candidate_value, value, pb.type)) {
                ++improvements;
                if (ctx.update)
                    ctx.update(0, candidate_value, iterations, elapsed(), 0);
            } else {
                wanted = std::min(max_free, 1.2 * wanted);
            }

            x.swap(candidate);
            value = candidate_value;
        } else {
            wanted = std::min(max_free, 1.2 * wanted);
        }
    }

    bx::notice(ctx,
               "  - lns optimizer: {} neighbourhoods, {} improvements, "
               "{:.0f} free variables at the end\n",
               iterations,
               improvements,
               wanted);

    best.solutions.assign(1, bx::solution(std::move(x), value));
    best.remaining_constraints = 0;
    best.duration = elapsed();
    best.loop = iterations;

    if (ctx.finish)
        ctx.finish(best);

    return best;
}

} // namespace

namespace baryonyx::itm {

result
lns_optimize(const context& ctx, const problem& pb)
{
    notice(ctx, "- lns optimizer starts\n");
    return ::optimize(ctx, pb);
}

} // namespace baryonyx::itm
