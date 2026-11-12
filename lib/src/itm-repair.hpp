/* Copyright (C) 2016-2021 INRAE
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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_REPAIR_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_REPAIR_HPP

#include <algorithm>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include "bit-array.hpp"
#include "itm-common.hpp"
#include "random-distributions.hpp"
#include "utils.hpp"

namespace baryonyx::itm {

/** The constraints of each variable, with its coefficients. */
class constraint_columns
{
public:
    constraint_columns(const std::vector<merged_constraint>& constraints,
                       const int variables)
      : m_start(static_cast<std::size_t>(variables) + 1u, 0)
    {
        for (const auto& cst : constraints)
            for (const auto& elem : cst.elements)
                ++m_start[static_cast<std::size_t>(elem.variable_index) + 1u];

        for (std::size_t j = 1; j < m_start.size(); ++j)
            m_start[j] += m_start[j - 1];

        m_elements.resize(static_cast<std::size_t>(m_start.back()));
        auto next = m_start;
        for (int k = 0, e = length(constraints); k != e; ++k)
            for (const auto& elem : constraints[k].elements)
                m_elements[static_cast<std::size_t>(
                  next[static_cast<std::size_t>(elem.variable_index)]++)] = {
                    k, elem.factor
                };
    }

    /** @return the (constraint, factor) pairs of @c variable. */
    std::span<const std::pair<int, int>> operator()(
      const int variable) const noexcept
    {
        const auto first = m_start[static_cast<std::size_t>(variable)];
        const auto last = m_start[static_cast<std::size_t>(variable) + 1u];

        return { m_elements.data() + first,
                 static_cast<std::size_t>(last - first) };
    }

private:
    std::vector<int> m_start;
    std::vector<std::pair<int, int>> m_elements;
};

/** A run ending with so few violated constraints is worth a repair. */
inline bool
is_repairable(const int remaining, const int constraints) noexcept
{
    const auto threshold =
      std::max(10,
               static_cast<int>(static_cast<double>(constraints) *
                                solver_parameters::repair_ratio));

    return remaining > 0 && remaining <= threshold;
}

inline long int
repair_flip_limit(const int variables, const int constraints) noexcept
{
    return std::min(solver_parameters::repair_flips *
                      (static_cast<long int>(variables) + constraints),
                    solver_parameters::repair_max_flips);
}

/**
 * Weighted flip search for a feasible solution, as in Feasibility Jump or
 * weighted WalkSAT: each step flips the variable of a random violated
 * constraint that most reduces the weighted violation, and the weights of
 * the constraints that stay violated grow at each local minimum. Ties are
 * broken by the cost of the flip.
 *
 * The weights and the tabu list persist between calls: the weights remember
 * which constraints are hard for the next starting solutions.
 */
template<typename Cost>
class weighted_flip_search
{
public:
    weighted_flip_search(const std::vector<merged_constraint>& constraints,
                         const Cost& costs,
                         const int variables)
      : m_constraints(constraints)
      , m_costs(costs)
      , m_columns(constraints, variables)
      , m_lhs(constraints.size(), 0)
      , m_position(constraints.size(), -1)
      , m_weight(constraints.size(), 1)
      , m_tabu_until(static_cast<std::size_t>(variables), 0)
    {}

    /**
     * Searches from @c x until no constraint is violated, @c flip_limit
     * flips are done, or @c stop returns true (checked every 1024 flips).
     *
     * @return the number of constraints still violated, 0 if @c x is now a
     *     solution.
     */
    template<typename Stop>
    int operator()(bit_array& x,
                   random_engine& rng,
                   const long int flip_limit,
                   const Stop& stop)
    {
        constexpr long int tabu_tenure = 10;

        const int m = length(m_constraints);

        m_violated.clear();
        std::fill(m_position.begin(), m_position.end(), -1);
        for (int k = 0; k != m; ++k) {
            long long value = 0;
            for (const auto& elem : m_constraints[k].elements)
                if (x[elem.variable_index])
                    value += elem.factor;
            m_lhs[static_cast<std::size_t>(k)] = value;
            update_violated(k);
        }

        for (long int done = 0; done < flip_limit && !m_violated.empty();
             ++done) {
            const long int step = ++m_step;
            if ((done & 1023) == 1023 && (is_stop_requested() || stop()))
                break;

            const int k = m_violated[static_cast<std::size_t>(uniform_below(
              rng, static_cast<std::uint64_t>(m_violated.size())))];
            const auto& elements = m_constraints[k].elements;

            int best = -1;
            long long best_delta = std::numeric_limits<long long>::max();
            real best_cost = std::numeric_limits<real>::infinity();
            std::uint64_t ties = 0;

            for (const auto& elem : elements) {
                const int j = elem.variable_index;
                if (m_tabu_until[static_cast<std::size_t>(j)] >= step)
                    continue;

                const int sign = x[j] ? -1 : 1;
                long long delta = 0;

                for (const auto [r, factor] : m_columns(j)) {
                    const auto value = m_lhs[static_cast<std::size_t>(r)];
                    delta += m_weight[static_cast<std::size_t>(r)] *
                             (violation(r, value + sign * factor) -
                              violation(r, value));
                }

                const real cost = static_cast<real>(sign) * m_costs(j, x);

                if (delta < best_delta ||
                    (delta == best_delta && cost < best_cost)) {
                    best = j;
                    best_delta = delta;
                    best_cost = cost;
                    ties = 1;
                } else if (delta == best_delta && cost == best_cost) {
                    // Reservoir sampling keeps a uniform choice among ties.
                    if (uniform_below(rng, ++ties) == 0)
                        best = j;
                }
            }

            if (best < 0 || best_delta >= 0) {
                for (const int v : m_violated)
                    ++m_weight[static_cast<std::size_t>(v)];

                best = elements[static_cast<std::size_t>(uniform_below(
                                  rng,
                                  static_cast<std::uint64_t>(elements.size())))]
                         .variable_index;
            }

            flip(x, best);
            m_tabu_until[static_cast<std::size_t>(best)] = step + tabu_tenure;
        }

        return length(m_violated);
    }

private:
    long long violation(const int k, const long long value) const noexcept
    {
        const auto& cst = m_constraints[static_cast<std::size_t>(k)];
        return value < cst.min   ? cst.min - value
               : value > cst.max ? value - cst.max
                                 : 0LL;
    }

    void update_violated(const int k)
    {
        const bool is_violated =
          violation(k, m_lhs[static_cast<std::size_t>(k)]) > 0;
        auto& pos = m_position[static_cast<std::size_t>(k)];

        if (is_violated && pos < 0) {
            pos = length(m_violated);
            m_violated.emplace_back(k);
        } else if (!is_violated && pos >= 0) {
            const int last = m_violated.back();
            m_violated[static_cast<std::size_t>(pos)] = last;
            m_position[static_cast<std::size_t>(last)] = pos;
            m_violated.pop_back();
            pos = -1;
        }
    }

    void flip(bit_array& x, const int j)
    {
        const int sign = x[j] ? -1 : 1;
        x.invert(j);

        for (const auto [k, factor] : m_columns(j)) {
            m_lhs[static_cast<std::size_t>(k)] += sign * factor;
            update_violated(k);
        }
    }

    const std::vector<merged_constraint>& m_constraints;
    const Cost& m_costs;
    constraint_columns m_columns;
    std::vector<long long> m_lhs;
    std::vector<int> m_violated;
    std::vector<int> m_position;
    std::vector<long long> m_weight;
    std::vector<long int> m_tabu_until;
    long int m_step = 0;
};

/**
 * Tries to satisfy the few constraints a solver run left violated. Wedelin
 * oscillates on these last constraints: their multipliers and preferences
 * compensate each other, and a weighted flip search escapes this cycle.
 *
 * @return the number of constraints still violated, 0 if @c x is now a
 *     solution.
 */
template<typename Cost>
int
repair(bit_array& x,
       const std::vector<merged_constraint>& constraints,
       const Cost& costs,
       const int variables,
       random_engine& rng,
       const long int flip_limit)
{
    weighted_flip_search<Cost> search(constraints, costs, variables);

    return search(x, rng, flip_limit, []() { return false; });
}

} // namespace baryonyx::itm

#endif
