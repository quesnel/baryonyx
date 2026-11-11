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
#include <utility>
#include <vector>

#include "bit-array.hpp"
#include "itm-common.hpp"
#include "random-distributions.hpp"
#include "utils.hpp"

namespace baryonyx::itm {

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
 * Local search that tries to satisfy the few constraints a solver run left
 * violated. Wedelin oscillates on these last constraints: their multipliers
 * and preferences compensate each other. A weighted flip search (as in
 * Feasibility Jump or weighted WalkSAT) escapes this cycle: each step flips
 * the variable of a violated constraint that most reduces the weighted
 * violation, and the weights of the constraints that stay violated grow at
 * each local minimum. Ties are broken by the cost of the flip.
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
    constexpr long int tabu_tenure = 10;

    const int m = length(constraints);

    // Column view of the constraints: for each variable, its constraints and
    // coefficients.
    std::vector<int> column_start(static_cast<std::size_t>(variables) + 1u,
                                  0);
    for (const auto& cst : constraints)
        for (const auto& elem : cst.elements)
            ++column_start[static_cast<std::size_t>(elem.variable_index) + 1u];

    for (std::size_t j = 1; j < column_start.size(); ++j)
        column_start[j] += column_start[j - 1];

    std::vector<std::pair<int, int>> column(
      static_cast<std::size_t>(column_start.back()));
    {
        auto next = column_start;
        for (int k = 0; k != m; ++k)
            for (const auto& elem : constraints[k].elements)
                column[static_cast<std::size_t>(
                  next[static_cast<std::size_t>(elem.variable_index)]++)] = {
                    k, elem.factor
                };
    }

    std::vector<long long> lhs(static_cast<std::size_t>(m), 0);
    for (int k = 0; k != m; ++k)
        for (const auto& elem : constraints[k].elements)
            if (x[elem.variable_index])
                lhs[static_cast<std::size_t>(k)] += elem.factor;

    const auto violation = [&constraints](int k, long long value) {
        const auto& cst = constraints[k];
        return value < cst.min   ? cst.min - value
               : value > cst.max ? value - cst.max
                                 : 0LL;
    };

    std::vector<int> violated;
    std::vector<int> position(static_cast<std::size_t>(m), -1);
    const auto update_violated = [&](int k) {
        const bool is_violated = violation(k, lhs[k]) > 0;
        auto& pos = position[static_cast<std::size_t>(k)];

        if (is_violated && pos < 0) {
            pos = length(violated);
            violated.emplace_back(k);
        } else if (!is_violated && pos >= 0) {
            const int last = violated.back();
            violated[static_cast<std::size_t>(pos)] = last;
            position[static_cast<std::size_t>(last)] = pos;
            violated.pop_back();
            pos = -1;
        }
    };

    for (int k = 0; k != m; ++k)
        update_violated(k);

    std::vector<long long> weight(static_cast<std::size_t>(m), 1);
    std::vector<long int> tabu_until(static_cast<std::size_t>(variables), 0);

    const auto flip = [&](int j) {
        const int sign = x[j] ? -1 : 1;
        x.invert(j);

        for (int c = column_start[static_cast<std::size_t>(j)],
                 e = column_start[static_cast<std::size_t>(j) + 1u];
             c != e;
             ++c) {
            const auto [k, factor] = column[static_cast<std::size_t>(c)];
            lhs[static_cast<std::size_t>(k)] += sign * factor;
            update_violated(k);
        }
    };

    for (long int step = 1; step <= flip_limit && !violated.empty(); ++step) {
        if ((step & 1023) == 0 && is_stop_requested())
            break;

        const int k = violated[static_cast<std::size_t>(
          uniform_below(rng, static_cast<std::uint64_t>(violated.size())))];
        const auto& elements = constraints[k].elements;

        int best = -1;
        long long best_delta = std::numeric_limits<long long>::max();
        real best_cost = std::numeric_limits<real>::infinity();
        std::uint64_t ties = 0;

        for (const auto& elem : elements) {
            const int j = elem.variable_index;
            if (tabu_until[static_cast<std::size_t>(j)] >= step)
                continue;

            const int sign = x[j] ? -1 : 1;
            long long delta = 0;

            for (int c = column_start[static_cast<std::size_t>(j)],
                     e = column_start[static_cast<std::size_t>(j) + 1u];
                 c != e;
                 ++c) {
                const auto [r, factor] = column[static_cast<std::size_t>(c)];
                const auto value = lhs[static_cast<std::size_t>(r)];
                delta += weight[static_cast<std::size_t>(r)] *
                         (violation(r, value + sign * factor) -
                          violation(r, value));
            }

            const real cost = static_cast<real>(sign) * costs(j, x);

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
            for (const int v : violated)
                ++weight[static_cast<std::size_t>(v)];

            best = elements[static_cast<std::size_t>(uniform_below(
                              rng, static_cast<std::uint64_t>(elements.size())))]
                     .variable_index;
        }

        flip(best);
        tabu_until[static_cast<std::size_t>(best)] = step + tabu_tenure;
    }

    return length(violated);
}

} // namespace baryonyx::itm

#endif
