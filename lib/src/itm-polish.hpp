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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_POLISH_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_POLISH_HPP

#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

#include "bit-array.hpp"
#include "itm-common.hpp"
#include "itm-repair.hpp"
#include "utils.hpp"

namespace baryonyx::itm {

/**
 * Descent to a local optimum of a feasible solution, for a minimized
 * objective. Wedelin stops where the multipliers balance, not at a local
 * optimum: on large problems some variables stay at a costly value that a
 * single flip, or a flip with a partner of the same constraint, removes
 * without violating any constraint.
 *
 * The column view and the buffers are built once and reused by each call.
 */
class polisher
{
public:
    polisher(const std::vector<merged_constraint>& constraints,
             const int variables)
      : m_constraints(constraints)
      , m_columns(constraints, variables)
      , m_lhs(constraints.size())
      , m_listed(static_cast<std::size_t>(variables), 0)
      , m_variables(variables)
    {}

    /**
     * Applies improving moves to the feasible solution @c x until none
     * remains.
     *
     * @param changed the variables changed since @c x was a local optimum,
     *     an empty span checks every variable.
     * @return the number of moves applied.
     */
    template<typename Cost>
    long int operator()(bit_array& x,
                        const Cost& costs,
                        std::span<const int> changed = {})
    {
        for (int k = 0, e = length(m_constraints); k != e; ++k) {
            long long value = 0;
            for (const auto& elem : m_constraints[k].elements)
                if (x[elem.variable_index])
                    value += elem.factor;
            m_lhs[static_cast<std::size_t>(k)] = value;
        }

        m_worklist.clear();
        if (changed.empty()) {
            for (int j = 0; j != m_variables; ++j)
                push(j);
        } else {
            for (const int j : changed)
                push_neighbours(j, costs);
        }

        // Each move strictly decreases the objective, the limit only guards
        // against the rounding errors of the quadratic costs.
        const long int move_limit =
          10 * (static_cast<long int>(m_variables) + length(m_constraints));
        long int moves = 0;
        long int evaluations = 0;

        while (!m_worklist.empty() && moves < move_limit) {
            if ((++evaluations & 1023) == 0 && is_stop_requested())
                break;

            const int j = m_worklist.back();
            m_worklist.pop_back();
            m_listed[static_cast<std::size_t>(j)] = 0;

            const auto delta = flip_delta(x, costs, j);
            if (!is_improving(delta, std::abs(delta)))
                continue;

            flip(x, j);
            if (is_feasible(j)) {
                ++moves;
                push_neighbours(j, costs);
                continue;
            }

            if (const int partner = best_partner(x, costs, j, delta);
                partner >= 0) {
                flip(x, partner);
                ++moves;
                push_neighbours(j, costs);
                push_neighbours(partner, costs);
                continue;
            }

            flip(x, j);
        }

        return moves;
    }

private:
    // Relative tolerance: the costs are float sums and an exchange of two
    // variables of equal costs must not count as an improvement.
    static bool is_improving(const real delta, const real magnitude) noexcept
    {
        return delta < -real{ 1e-6 } * std::max(real{ 1 }, magnitude);
    }

    /** Objective change when @c j flips. */
    template<typename Cost>
    static real flip_delta(bit_array& x, const Cost& costs, const int j)
    {
        if (x[j])
            return -costs(j, x);

        // The cost of a variable includes its x_j * x_j term only while
        // x_j = 1.
        x.set(j);
        const auto delta = costs(j, x);
        x.unset(j);
        return delta;
    }

    void flip(bit_array& x, const int j) noexcept
    {
        const int sign = x[j] ? -1 : 1;
        x.invert(j);

        for (const auto [k, factor] : m_columns(j))
            m_lhs[static_cast<std::size_t>(k)] += sign * factor;
    }

    bool is_satisfied(const int k) const noexcept
    {
        const auto value = m_lhs[static_cast<std::size_t>(k)];
        const auto& cst = m_constraints[static_cast<std::size_t>(k)];

        return cst.min <= value && value <= cst.max;
    }

    bool is_feasible(const int j) const noexcept
    {
        for (const auto [k, factor] : m_columns(j))
            if (!is_satisfied(k))
                return false;

        return true;
    }

    /**
     * With @c j already flipped and some of its constraints violated, finds
     * the variable of the first violated constraint whose flip restores
     * feasibility with the best improvement of the pair.
     *
     * @return the partner, -1 if no flip of a partner improves.
     */
    template<typename Cost>
    int best_partner(bit_array& x,
                     const Cost& costs,
                     const int j,
                     const real delta_j)
    {
        int violated = -1;
        for (const auto [k, factor] : m_columns(j)) {
            if (!is_satisfied(k)) {
                violated = k;
                break;
            }
        }

        int best = -1;
        real best_delta = 0;

        for (const auto& elem :
             m_constraints[static_cast<std::size_t>(violated)].elements) {
            const int i = elem.variable_index;
            if (i == j)
                continue;

            const auto delta_i = flip_delta(x, costs, i);
            const auto delta = delta_j + delta_i;
            if (!is_improving(delta, std::abs(delta_j) + std::abs(delta_i)) ||
                (best >= 0 && delta >= best_delta))
                continue;

            flip(x, i);
            if (is_feasible(i) && is_feasible(j)) {
                best = i;
                best_delta = delta;
            }
            flip(x, i);
        }

        return best;
    }

    void push(const int j)
    {
        if (!m_listed[static_cast<std::size_t>(j)]) {
            m_listed[static_cast<std::size_t>(j)] = 1;
            m_worklist.emplace_back(j);
        }
    }

    /**
     * Lists the variables whose moves may have changed with the flip of
     * @c j: those sharing a constraint with @c j and, for a quadratic
     * objective, those sharing a term with it.
     */
    template<typename Cost>
    void push_neighbours(const int j, const Cost& costs)
    {
        push(j);

        for (const auto [k, factor] : m_columns(j))
            for (const auto& elem :
                 m_constraints[static_cast<std::size_t>(k)].elements)
                push(elem.variable_index);

        if constexpr (requires { costs.quadratic_elements[0].id; }) {
            for (auto q = costs.indices[j], e = costs.indices[j + 1]; q != e;
                 ++q)
                push(costs.quadratic_elements[q].id);
        }
    }

    const std::vector<merged_constraint>& m_constraints;
    constraint_columns m_columns;
    std::vector<long long> m_lhs;
    std::vector<char> m_listed;
    std::vector<int> m_worklist;
    int m_variables;
};

} // namespace baryonyx::itm

#endif
