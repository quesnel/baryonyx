/* Copyright (C) 2016-2026 INRAE
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
#include <array>
#include <chrono>
#include <cmath>
#include <functional>
#include <numeric>
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
        walk(rng, wanted, is_free, visited, queue, [&]() {
            return static_cast<int>(bx::uniform_below(
              rng, static_cast<std::uint64_t>(constraints())));
        });
    }

    /**
     * As walk, but each connected component starts from a constraint of a
     * variable drawn by pick.
     */
    template<typename PickVariable>
    void walk_from(bx::xoshiro256pp& rng,
                   int wanted,
                   std::vector<char>& is_free,
                   std::vector<char>& visited,
                   std::vector<int>& queue,
                   PickVariable&& pick) const
    {
        walk(rng, wanted, is_free, visited, queue, [&]() {
            for (;;) {
                const auto& csts =
                  m_variable_constraints[static_cast<std::size_t>(pick())];
                if (!csts.empty())
                    return csts[bx::uniform_below(rng, csts.size())];
            }
        });
    }

    template<typename PickConstraint>
    void walk(bx::xoshiro256pp& rng,
              int wanted,
              std::vector<char>& is_free,
              std::vector<char>& visited,
              std::vector<int>& queue,
              PickConstraint&& pick) const
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

                const auto start = pick();
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

/**
 * Key variables: set to 1, one of them saturates packing or partitioning rows
 * and forces many others to 0, as the tree of a land unit fixes its crops.
 * Changing a few of them together and solving every other variable again
 * links solutions that the neighbourhoods of free variables cannot: their
 * walk spreads through the large rows and seldom frees a key variable with
 * all the variables that depend on it.
 */
class key_variables
{
public:
    explicit key_variables(const bx::problem& pb)
    {
        const auto n = pb.vars.values.size();
        std::vector<int> reach(n, 0);
        std::vector<int> arities;

        for_each_constraint(pb, [&](const bx::constraint& cst, sense) {
            arities.emplace_back(bx::length(cst.elements));
        });

        if (arities.empty())
            return;

        std::ranges::nth_element(arities, arities.begin() + arities.size() / 2);
        const int small_arity =
          std::max(16, 4 * arities[arities.size() / 2]);

        // Only small rows of nonnegative factors: a variable whose factor
        // equals the right-hand side forces the others to 0. A large row
        // saturated by one variable is a global limit, not a local
        // dependency a region of key variables can move with.
        for_each_constraint(pb, [&](const bx::constraint& cst, sense s) {
            if (s == sense::greater ||
                bx::length(cst.elements) > small_arity)
                return;
            if (std::ranges::any_of(cst.elements,
                                    [](const auto& e) { return e.factor < 0; }))
                return;

            const int others = bx::length(cst.elements) - 1;
            for (const auto& e : cst.elements)
                if (e.factor == cst.value)
                    reach[static_cast<std::size_t>(e.variable_index)] += others;
        });

        std::vector<int> sorted(reach);
        std::ranges::nth_element(sorted, sorted.begin() + sorted.size() / 2);
        const int threshold = std::max(8, 4 * sorted[sorted.size() / 2]);

        m_position.assign(n, -1);
        for (std::size_t i = 0; i != n; ++i) {
            if (reach[i] >= threshold) {
                m_position[i] = bx::length(m_keys);
                m_keys.emplace_back(static_cast<int>(i));
            }
        }

        if (bx::length(m_keys) < 4 ||
            m_keys.size() > n / 4) {
            m_keys.clear();
            return;
        }

        m_neighbours.resize(m_keys.size());
        m_key_constraints_of.resize(m_keys.size());

        std::vector<int> in_row;
        for_each_constraint(pb, [&](const bx::constraint& cst, sense s) {
            in_row.clear();
            for (const auto& e : cst.elements)
                if (const int p = key(e.variable_index); p >= 0)
                    in_row.emplace_back(p);

            if (in_row.size() == cst.elements.size()) {
                const int id = bx::length(m_key_constraints);
                m_key_constraints.push_back({ &cst, s });
                for (const int p : in_row)
                    m_key_constraints_of[static_cast<std::size_t>(p)]
                      .emplace_back(id);
            }

            if (bx::length(cst.elements) > small_arity)
                return;

            for (const int a : in_row)
                for (const int b : in_row)
                    if (a != b)
                        m_neighbours[static_cast<std::size_t>(a)].emplace_back(
                          b);
        });

        for (auto& nb : m_neighbours) {
            std::ranges::sort(nb);
            const auto [first, last] = std::ranges::unique(nb);
            nb.erase(first, last);
        }
    }

    bool empty() const noexcept { return m_keys.empty(); }

    int size() const noexcept { return bx::length(m_keys); }

    int variable(int position) const noexcept
    {
        return m_keys[static_cast<std::size_t>(position)];
    }

    void disable() noexcept { m_keys.clear(); }

    /** Breadth-first walk of @c k key variables from a random one. */
    void region(bx::xoshiro256pp& rng, int k, std::vector<int>& out) const
    {
        out.clear();
        std::vector<char> seen(m_keys.size(), 0);
        const auto start = static_cast<int>(
          bx::uniform_below(rng, static_cast<std::uint64_t>(m_keys.size())));
        out.emplace_back(start);
        seen[static_cast<std::size_t>(start)] = 1;

        for (std::size_t head = 0; head != out.size() && bx::length(out) < k;
             ++head) {
            auto nb = m_neighbours[static_cast<std::size_t>(out[head])];
            for (std::size_t i = nb.size(); i > 1; --i)
                std::swap(nb[i - 1], nb[bx::uniform_below(rng, i)]);

            for (const int p : nb) {
                if (bx::length(out) == k)
                    break;
                if (!seen[static_cast<std::size_t>(p)]) {
                    seen[static_cast<std::size_t>(p)] = 1;
                    out.emplace_back(p);
                }
            }
        }
    }

    /**
     * The assignments of the key variables of @c region, the other variables
     * keeping their value in @c x, that satisfy the rows of key variables
     * only. At most @c limit of them, chosen uniformly by reservoir sampling.
     */
    std::vector<std::uint32_t> assignments(const std::vector<std::int8_t>& x,
                                           const std::vector<int>& region,
                                           int limit,
                                           bx::xoshiro256pp& rng) const
    {
        std::vector<std::int8_t> value(m_keys.size());
        for (std::size_t p = 0; p != m_keys.size(); ++p)
            value[p] = x[static_cast<std::size_t>(m_keys[p])];

        std::vector<char> open(m_keys.size(), 0);
        for (const int p : region)
            open[static_cast<std::size_t>(p)] = 1;

        std::vector<std::uint32_t> found;
        long int seen = 0;

        const auto satisfiable = [&](int p) {
            for (const int id : m_key_constraints_of[static_cast<std::size_t>(p)]) {
                const auto& [cst, s] =
                  m_key_constraints[static_cast<std::size_t>(id)];
                int low = 0, high = 0;
                for (const auto& e : cst->elements) {
                    const auto q = static_cast<std::size_t>(
                      key(e.variable_index));
                    if (open[q]) {
                        low += std::min(0, e.factor);
                        high += std::max(0, e.factor);
                    } else {
                        low += e.factor * value[q];
                        high += e.factor * value[q];
                    }
                }
                if ((s != sense::greater && low > cst->value) ||
                    (s != sense::less && high < cst->value))
                    return false;
            }
            return true;
        };

        const auto depth_first = [&](const auto& self,
                                     std::size_t depth,
                                     std::uint32_t bits) -> void {
            if (depth == region.size()) {
                ++seen;
                if (bx::length(found) < limit)
                    found.emplace_back(bits);
                else if (const auto r = bx::uniform_below(
                           rng, static_cast<std::uint64_t>(seen));
                         r < static_cast<std::uint64_t>(limit))
                    found[r] = bits;
                return;
            }

            const auto p = static_cast<std::size_t>(region[depth]);
            open[p] = 0;
            for (std::int8_t v = 0; v != 2; ++v) {
                value[p] = v;
                if (satisfiable(static_cast<int>(p)))
                    self(self,
                         depth + 1,
                         bits | (static_cast<std::uint32_t>(v) << depth));
            }
            open[p] = 1;
        };

        depth_first(depth_first, 0, 0u);
        return found;
    }

private:
    enum class sense { equal, greater, less };

    template<typename Function>
    static void for_each_constraint(const bx::problem& pb, Function&& f)
    {
        for (const auto& cst : pb.equal_constraints)
            f(cst, sense::equal);
        for (const auto& cst : pb.greater_constraints)
            f(cst, sense::greater);
        for (const auto& cst : pb.less_constraints)
            f(cst, sense::less);
    }

    int key(int variable) const noexcept
    {
        return m_position[static_cast<std::size_t>(variable)];
    }

    struct key_constraint
    {
        const bx::constraint* cst;
        sense s;
    };

    std::vector<int> m_keys;
    std::vector<int> m_position;
    std::vector<std::vector<int>> m_neighbours;
    std::vector<key_constraint> m_key_constraints;
    std::vector<std::vector<int>> m_key_constraints_of;
};

enum class neighbourhood : std::uint8_t
{
    random,
    walk,
    hesitant,
    hesitant_walk
};

constexpr std::size_t neighbourhood_count = 4;

constexpr std::array<std::string_view, neighbourhood_count>
  neighbourhood_names = { "random", "walk", "hesitant", "hesitant walk" };

/**
 * Adaptive large neighbourhood search: draws the kind of neighbourhood in
 * proportion to its recent rate of improvement. The kind that works differs
 * between problems (the hesitant variables find a third of the improvements
 * on some, none on others) and a fixed rotation wastes the iterations of
 * the useless ones.
 */
class neighbourhood_selector
{
public:
    explicit neighbourhood_selector(bool has_constraints) noexcept
    {
        m_enabled.fill(true);
        if (!has_constraints) {
            m_enabled[static_cast<std::size_t>(neighbourhood::walk)] = false;
            m_enabled[static_cast<std::size_t>(
              neighbourhood::hesitant_walk)] = false;
        }

        // Improvements are rare: the prior of 5% keeps a kind with one lucky
        // success from taking every iteration.
        m_successes.fill(0.5);
        m_tries.fill(10.0);
    }

    neighbourhood draw(bx::xoshiro256pp& rng) const noexcept
    {
        std::array<double, neighbourhood_count> weight{};
        double rate_sum = 0.0;
        int enabled = 0;
        for (std::size_t i = 0; i != neighbourhood_count; ++i) {
            if (m_enabled[i]) {
                weight[i] = m_successes[i] / m_tries[i];
                rate_sum += weight[i];
                ++enabled;
            }
        }

        // A tenth of the draws ignore the rates: a kind can become useful
        // later, once the others have exhausted their improvements.
        double total = 0.0;
        for (std::size_t i = 0; i != neighbourhood_count; ++i) {
            if (m_enabled[i])
                weight[i] = 0.1 / enabled + 0.9 * weight[i] / rate_sum;
            total += weight[i];
        }

        bx::uniform_real_distribution<double> unit(0.0, total);
        double u = unit(rng);
        std::size_t last = 0;
        for (std::size_t i = 0; i != neighbourhood_count; ++i) {
            if (!m_enabled[i])
                continue;
            last = i;
            if (u < weight[i])
                return static_cast<neighbourhood>(i);
            u -= weight[i];
        }

        return static_cast<neighbourhood>(last);
    }

    void update(neighbourhood kind, bool improved) noexcept
    {
        // Forgets the rates over about a hundred uses: the useful kind
        // changes as the search goes deeper.
        constexpr double decay = 0.99;

        const auto i = static_cast<std::size_t>(kind);
        m_successes[i] = decay * m_successes[i] + (improved ? 1.0 : 0.0);
        m_tries[i] = decay * m_tries[i] + 1.0;
        ++m_uses[i];
        if (improved)
            ++m_improvements[i];
    }

    long int uses(std::size_t i) const noexcept { return m_uses[i]; }

    long int improvements(std::size_t i) const noexcept
    {
        return m_improvements[i];
    }

private:
    std::array<double, neighbourhood_count> m_successes;
    std::array<double, neighbourhood_count> m_tries;
    std::array<long int, neighbourhood_count> m_uses{};
    std::array<long int, neighbourhood_count> m_improvements{};
    std::array<bool, neighbourhood_count> m_enabled;
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
    bx::affecter affect(internal, pb);

    const double min_free = std::min<double>(20.0, static_cast<double>(n));
    const double max_free = std::max(min_free, 0.5 * static_cast<double>(n));
    double wanted =
      std::clamp(0.1 * static_cast<double>(n), min_free, max_free);

    // Short runs: most neighbourhoods are solved within a few milliseconds
    // and more neighbourhoods beat longer runs on each.
    internal.parameters.time_limit = 0.02;
    internal.parameters.thread = 2;

    std::vector<char> is_free(n);
    // Hesitation of each variable in [0, 1], learned from the Wedelin margins
    // of the sub-problems where it was free: 1 when it had the smallest
    // margin of its sub-problem. A variable never freed keeps the middle.
    std::vector<double> margins;
    internal.margins = &margins;
    std::vector<float> hesitation(n, 0.5f);
    std::vector<int> order;
    const auto learn_hesitation = [&](const bx::problem& sub) {
        const auto m = margins.size();
        if (m < 2 || m != sub.vars.names.size())
            return;

        order.resize(m);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return margins[static_cast<std::size_t>(a)] <
                   margins[static_cast<std::size_t>(b)];
        });

        for (std::size_t r = 0; r != m; ++r) {
            const auto it =
              index.find(sub.vars.names[static_cast<std::size_t>(order[r])]);
            if (it == index.end())
                continue;
            const auto score = 1.0f - static_cast<float>(r) /
                                        static_cast<float>(m - 1);
            auto& h = hesitation[static_cast<std::size_t>(it->second)];
            h = 0.5f * (h + score);
        }
    };

    // The variables Wedelin hesitates on are drawn with probabilities rising
    // with their hesitation: RINS frees the variables where the incumbent
    // and the relaxation disagree. Scattered, they suit problems of loosely
    // linked variables, grown into a walk, problems of tight constraints.
    bx::uniform_real_distribution<float> unit(0.0f, 1.0f);
    const auto hesitation_weight = [&](std::size_t v) {
        const auto h = hesitation[v];
        return 0.05f + h * h;
    };
    const auto pick_hesitant = [&]() {
        for (;;) {
            const auto v = bx::uniform_below(rng, n);
            if (unit(rng) * 1.05f <= hesitation_weight(v))
                return static_cast<int>(v);
        }
    };

    // Weighted sampling without replacement (Efraimidis and Spirakis).
    std::vector<std::pair<float, int>> keyed;
    const auto free_hesitant = [&](int k) {
        keyed.resize(n);
        for (std::size_t i = 0; i != n; ++i) {
            const auto u = std::max(1e-12f, unit(rng));
            keyed[i] = { std::log(u) / hesitation_weight(i),
                         static_cast<int>(i) };
        }
        const auto kk = std::min(static_cast<std::size_t>(k), n);
        std::nth_element(keyed.begin(),
                         keyed.begin() + static_cast<std::ptrdiff_t>(kk),
                         keyed.end(),
                         std::greater<>());
        for (std::size_t i = 0; i != kk; ++i)
            is_free[static_cast<std::size_t>(keyed[i].second)] = 1;
    };

    neighbourhood_selector selector(graph.constraints() != 0);

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

    const auto key_assignments = [&](const key_variables& keys) {
        assignments.clear();
        for (int p = 0, e = keys.size(); p != e; ++p) {
            const int v = keys.variable(p);
            assignments.emplace_back(v, x[static_cast<std::size_t>(v)] != 0);
        }
    };

    // The key variables are worth moving only if fixing them propagates to
    // many other variables: otherwise each assignment of a region leaves
    // almost the whole problem to solve again.
    key_variables keys(pb);
    if (!keys.empty()) {
        key_assignments(keys);
        try {
            const auto rest = affect(assignments).vars.values.size();
            const auto propagated = n - rest - keys.size();
            if (propagated < 4 * static_cast<std::size_t>(keys.size()))
                keys.disable();
        } catch (const bx::solver_failure&) {
            keys.disable();
        }
    }

    // The best solution: a key move may accept a worse x whose key
    // variables promise better once the other variables are improved.
    std::vector<std::int8_t> incumbent = x;
    double incumbent_value = value;
    const auto keep_best = [&]() {
        if (!is_better(value, incumbent_value, pb.type))
            return;
        incumbent = x;
        incumbent_value = value;
        if (ctx.update)
            ctx.update(0, value, iterations, elapsed(), 0);
    };

    // Fixes the key variables to x, except the region set to bits, and
    // solves the other variables again in candidate.
    std::vector<int> region;
    const auto solve_keys = [&](std::uint32_t bits) -> std::optional<double> {
        key_assignments(keys);
        for (std::size_t d = 0; d != region.size(); ++d)
            assignments[static_cast<std::size_t>(region[d])].second =
              ((bits >> d) & 1u) != 0;

        bx::problem sub;
        try {
            sub = affect(assignments);
        } catch (const bx::solver_failure&) {
            return std::nullopt;
        }

        if (sub.vars.values.empty())
            return std::nullopt;

        const auto ret = bx::itm::optimize(internal, sub);
        if (!ret || ret.solutions.empty())
            return std::nullopt;

        candidate = x;
        merge(sub, ret, index, candidate);
        if (!bx::is_valid_solution(pb, candidate))
            return std::nullopt;

        changed.clear();
        for (std::size_t i = 0; i != n; ++i)
            if (candidate[i] != x[i])
                changed.emplace_back(static_cast<int>(i));
        if (!changed.empty())
            full.polish(candidate, changed);

        return bx::compute_solution(pb, candidate);
    };

    // Enumerates the assignments of a region of key variables and solves
    // the other variables again for each. The reference is the current
    // assignment of the region solved the same way, not x itself: the other
    // variables of x had many neighbourhoods to improve, a candidate only one
    // short run, and a better assignment of the key variables would lose.
    std::vector<std::int8_t> best_candidate;
    const auto key_move = [&]() {
        constexpr int region_size = 16;
        constexpr int candidate_limit = 64;

        keys.region(rng, std::min(region_size, keys.size()), region);
        const auto candidates =
          keys.assignments(x, region, candidate_limit, rng);

        std::uint32_t current = 0;
        for (std::size_t d = 0; d != region.size(); ++d)
            if (x[static_cast<std::size_t>(keys.variable(region[d]))])
                current |= 1u << d;

        const auto reference = solve_keys(current);
        if (!reference)
            return false;

        double best_value = *reference;
        best_candidate.clear();

        for (const auto bits : candidates) {
            if (bx::is_stop_requested() || !has_time())
                break;
            if (bits == current)
                continue;

            if (const auto v = solve_keys(bits);
                v && is_better(*v, best_value, pb.type)) {
                best_value = *v;
                best_candidate = candidate;
            }
        }

        if (best_candidate.empty())
            return false;

        x.swap(best_candidate);
        value = best_value;
        keep_best();
        return true;
    };

    // Share of the search time given to the moves of key variables: raised
    // after an improvement and slowly reduced after a region without one,
    // so that the problems where they never help lose little time. They only
    // start once the neighbourhoods stall: while the lns still improves, the
    // key layout is not yet the bottleneck and moving it only adds variance.
    constexpr long int stall_limit = 100;
    double key_share = 0.25;
    double key_time = 0.0;
    const double loop_start = elapsed();
    long int key_moves = 0, key_improvements = 0;
    long int last_improvement = 0;

    while (!bx::is_stop_requested() && has_time()) {
        if (!keys.empty() && iterations - last_improvement >= stall_limit &&
            key_time <= key_share * (elapsed() - loop_start)) {
            const double move_start = elapsed();
            ++key_moves;
            if (key_move()) {
                ++key_improvements;
                last_improvement = iterations;
                key_share = std::min(0.6, 1.5 * key_share);
            } else {
                key_share = std::max(0.02, 0.95 * key_share);
            }
            key_time += elapsed() - move_start;
            continue;
        }

        ++iterations;
        std::fill(is_free.begin(), is_free.end(), 0);
        const int k = static_cast<int>(wanted);

        const auto kind = selector.draw(rng);
        switch (kind) {
        case neighbourhood::random:
            for (int freed = 0; freed < k;) {
                const auto v = bx::uniform_below(rng, n);
                if (!is_free[v]) {
                    is_free[v] = 1;
                    ++freed;
                }
            }
            break;
        case neighbourhood::walk:
            graph.walk(rng, k, is_free, visited, queue);
            break;
        case neighbourhood::hesitant:
            free_hesitant(k);
            break;
        case neighbourhood::hesitant_walk:
            graph.walk_from(rng, k, is_free, visited, queue, pick_hesitant);
            break;
        }

        assignments.clear();
        for (std::size_t i = 0; i != n; ++i)
            if (!is_free[i])
                assignments.emplace_back(static_cast<int>(i), x[i] != 0);

        bx::problem sub;
        try {
            sub = affect(assignments);
        } catch (const bx::solver_failure&) {
            selector.update(kind, false);
            continue;
        }

        if (sub.vars.values.empty()) {
            selector.update(kind, false);
            continue;
        }

        const auto ret = bx::itm::optimize(internal, sub);
        if (ret)
            learn_hesitation(sub);
        if (!ret || ret.solutions.empty()) {
            selector.update(kind, false);
            wanted = std::max(min_free, 0.8 * wanted);
            continue;
        }

        candidate = x;
        merge(sub, ret, index, candidate);

        if (!bx::is_valid_solution(pb, candidate)) {
            selector.update(kind, false);
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

        const bool improved = is_better(candidate_value, value, pb.type);
        selector.update(kind, improved);

        if (is_better(value, candidate_value, pb.type)) {
            wanted = std::min(max_free, 1.2 * wanted);
            continue;
        }

        if (improved) {
            ++improvements;
            last_improvement = iterations;
        } else
            wanted = std::min(max_free, 1.2 * wanted);

        x.swap(candidate);
        value = candidate_value;
        keep_best();
    }

    bx::notice(ctx,
               "  - lns optimizer: {} neighbourhoods, {} improvements, "
               "{:.0f} free variables at the end\n",
               iterations,
               improvements,
               wanted);

    for (std::size_t i = 0; i != neighbourhood_count; ++i)
        if (selector.uses(i) > 0)
            bx::notice(ctx,
                       "    {}: {} neighbourhoods, {} improvements\n",
                       neighbourhood_names[i],
                       selector.uses(i),
                       selector.improvements(i));

    if (!keys.empty())
        bx::notice(ctx,
                   "  - lns optimizer: {} key variables, {} key moves, {} "
                   "improvements, {:.0f}% of the time\n",
                   keys.size(),
                   key_moves,
                   key_improvements,
                   100.0 * key_time / std::max(1e-9, elapsed() - loop_start));

    best.solutions.assign(1,
                          bx::solution(std::move(incumbent), incumbent_value));
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
