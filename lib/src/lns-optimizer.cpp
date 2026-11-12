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
#include <optional>
#include <span>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "debug.hpp"
#include "itm-polish.hpp"
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
        int visited_count = 0;
        std::size_t head = 0;

        while (freed < wanted) {
            if (head == queue.size()) {
                // Every constraint is visited but variables outside any
                // constraint remain: complete with random variables.
                if (visited_count == constraints()) {
                    while (freed < wanted) {
                        const auto v = bx::uniform_below(rng, is_free.size());
                        if (!is_free[v]) {
                            is_free[v] = 1;
                            ++freed;
                        }
                    }
                    return;
                }

                const auto start = static_cast<int>(bx::uniform_below(
                  rng, static_cast<std::uint64_t>(constraints())));
                if (visited[static_cast<std::size_t>(start)])
                    continue;
                visited[static_cast<std::size_t>(start)] = 1;
                ++visited_count;
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
                        ++visited_count;
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

/**
 * The full problem in the form of the local searches: merged constraints and
 * a minimized objective, built once for the polish of each candidate and
 * for the feasibility jump phase.
 */
class full_problem
{
public:
    full_problem(const bx::context& ctx, const bx::problem& pb)
      : m_objective(bx::itm::minimization_objective(pb))
      , m_constraints(bx::itm::make_merged_constraints(ctx, pb))
      , m_polish(m_constraints, bx::length(pb.vars.values))
      , m_x(bx::length(pb.vars.values))
    {
        const auto n = bx::length(pb.vars.values);

        if (m_objective.qelements.empty())
            m_linear.emplace(m_objective, n);
        else
            m_quadratic.emplace(m_objective, n);
    }

    // The cost functions and the searches keep references to the members.
    full_problem(const full_problem&) = delete;
    full_problem& operator=(const full_problem&) = delete;

    /**
     * Descent to a local optimum: a variable that no recent neighbourhood
     * freed keeps a costly value the sub-problems never see.
     *
     * @param changed the variables changed since @c x was a local optimum,
     *     an empty span checks every variable.
     * @return true if @c x changed.
     */
    bool polish(std::vector<std::int8_t>& x, std::span<const int> changed)
    {
        load(x);

        const auto moves = m_linear ? m_polish(m_x, *m_linear, changed)
                                    : m_polish(m_x, *m_quadratic, changed);
        if (moves == 0)
            return false;

        store(x);
        return true;
    }

    /**
     * Weighted flip search from @c x until it is feasible or @c stop returns
     * true. The constraint weights persist between calls.
     *
     * @return true if @c x is now a solution.
     */
    template<typename Stop>
    bool jump(std::vector<std::int8_t>& x,
              bx::xoshiro256pp& rng,
              const Stop& stop)
    {
        if (m_linear && !m_linear_search)
            m_linear_search.emplace(m_constraints, *m_linear, m_x.size());
        if (m_quadratic && !m_quadratic_search)
            m_quadratic_search.emplace(
              m_constraints, *m_quadratic, m_x.size());

        load(x);

        const auto flip_limit = bx::itm::repair_flip_limit(
          m_x.size(), bx::length(m_constraints));
        int remaining = 1;
        while (remaining != 0 && !stop() && !bx::is_stop_requested())
            remaining =
              m_linear_search
                ? (*m_linear_search)(m_x, rng, flip_limit, stop)
                : (*m_quadratic_search)(m_x, rng, flip_limit, stop);

        store(x);
        return remaining == 0;
    }

private:
    void load(const std::vector<std::int8_t>& x)
    {
        for (int j = 0, e = m_x.size(); j != e; ++j)
            m_x.set(j, x[static_cast<std::size_t>(j)] != 0);
    }

    void store(std::vector<std::int8_t>& x) const
    {
        for (int j = 0, e = m_x.size(); j != e; ++j)
            x[static_cast<std::size_t>(j)] = static_cast<std::int8_t>(m_x[j]);
    }

    bx::objective_function m_objective;
    std::vector<bx::itm::merged_constraint> m_constraints;
    bx::itm::polisher m_polish;
    std::optional<bx::itm::default_cost_type> m_linear;
    std::optional<bx::itm::quadratic_cost_type> m_quadratic;
    std::optional<bx::itm::weighted_flip_search<bx::itm::default_cost_type>>
      m_linear_search;
    std::optional<bx::itm::weighted_flip_search<bx::itm::quadratic_cost_type>>
      m_quadratic_search;
    bx::bit_array m_x;
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

    internal.parameters = ctx.parameters;
    full_problem full(internal, pb);
    bx::xoshiro256pp rng(
      ctx.parameters.seed > 0
        ? static_cast<bx::xoshiro256pp::result_type>(ctx.parameters.seed)
        : 1u);

    const auto n = pb.vars.values.size();
    std::unordered_map<std::string_view, int> index;
    index.reserve(n);
    for (std::size_t i = 0; i != n; ++i)
        index.emplace(pb.vars.names[i], static_cast<int>(i));

    bx::result best;
    std::string_view first_by;
    double first_at = 0.0;

    const auto run = [&](optimizer_type optimizer, double time_limit) {
        if (bx::is_stop_requested())
            return;

        internal.parameters = ctx.parameters;
        internal.parameters.optimizer = optimizer;
        internal.parameters.time_limit = time_limit;

        auto ret = bx::itm::optimize(internal, pb);
        if (ret || ret.remaining_constraints < best.remaining_constraints)
            best = std::move(ret);
        if (best) {
            first_by = optimizer == optimizer_type::tuned ? "tuned" : "restart";
            first_at = elapsed();
        }
    };

    // Searches feasibility only, from the assignment of the best infeasible
    // run: Wedelin can end with violated constraints on problems where a
    // weighted flip search finds a solution in a fraction of a second.
    std::vector<std::int8_t> x(n, 0);
    const auto jump = [&](double until) {
        if (bx::is_stop_requested() || best || n == 0)
            return;

        if (!best.solutions.empty()) {
            std::fill(x.begin(), x.end(), 0);
            merge(pb, best, index, x);
        }

        if (!full.jump(x, rng, [&]() { return elapsed() >= until; }) ||
            !bx::is_valid_solution(pb, x))
            return;

        best.status = bx::result_status::success;
        best.strings = pb.strings;
        best.affected_vars = pb.affected_vars;
        best.variable_name = pb.vars.names;
        best.variables = static_cast<int>(n);
        best.remaining_constraints = 0;
        best.solutions.assign(1, bx::solution(x, bx::compute_solution(pb, x)));
        first_by = "feasibility jump";
        first_at = elapsed();
    };

    // A first solution with restart, then with feasibility jump, then with
    // tuned parameters, which find one on parameter-sensitive problems. With
    // a time limit: a quarter for restart, a tenth for feasibility jump,
    // which succeeds within a fraction of a second or not at all, and half
    // for tuned, so that lns keeps time to improve a late first solution. A
    // candidate of the tuned optimizer must not tune again and keeps
    // searching with restart instead, without feasibility jump which would
    // take its time and hide the quality of its parameters. Without time
    // limit: repeat the three with doubling durations until a solution is
    // found.
    if (unbounded) {
        for (double t = 1.0; !best && !bx::is_stop_requested(); t *= 2.0) {
            run(optimizer_type::restart, t);
            jump(elapsed() + t);
            if (!best)
                run(optimizer_type::tuned, 2.0 * t);
        }
    } else {
        run(optimizer_type::restart, 0.25 * budget);
        // A time limit of zero means no limit for the optimizers.
        if (ctx.tuning) {
            if (const double left = budget - elapsed(); !best && left > 0.0)
                run(optimizer_type::restart, left);
        } else {
            jump(elapsed() + 0.1 * budget);
            if (const double left = budget - elapsed(); !best && left > 0.0)
                run(optimizer_type::tuned, std::min(0.5 * budget, left));
        }
    }

    internal.parameters = ctx.parameters;
    internal.parameters.optimizer = optimizer_type::restart;

    if (!best || n == 0 || best.solutions.empty() ||
        best.solutions.back().variables.size() != n) {
        if (ctx.finish)
            ctx.finish(best);
        return best;
    }

    bx::notice(ctx,
               "  - lns optimizer: first solution at {:.3f}s by {}\n",
               first_at,
               first_by);

    x = best.solutions.back().variables;
    full.polish(x, {});
    double value = bx::compute_solution(pb, x);

    if (ctx.update)
        ctx.update(0, value, 0, elapsed(), 0);

    const constraint_graph graph(pb);

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
    std::vector<int> changed;
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

        if (!bx::is_valid_solution(pb, candidate)) {
            wanted = std::min(max_free, 1.2 * wanted);
            continue;
        }

        changed.clear();
        for (std::size_t i = 0; i != n; ++i)
            if (candidate[i] != x[i])
                changed.emplace_back(static_cast<int>(i));

        if (!changed.empty())
            full.polish(candidate, changed);

        const double candidate_value = bx::compute_solution(pb, candidate);

        if (is_better(value, candidate_value, pb.type)) {
            wanted = std::min(max_free, 1.2 * wanted);
            continue;
        }

        if (is_better(candidate_value, value, pb.type)) {
            ++improvements;
            if (ctx.update)
                ctx.update(0, candidate_value, iterations, elapsed(), 0);
        } else {
            wanted = std::min(max_free, 1.2 * wanted);
        }

        x.swap(candidate);
        value = candidate_value;
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
