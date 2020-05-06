/* Copyright (C) 2016-2019 INRA
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

#ifndef ORG_VLEPROJECT_BARYONYX_EXHAUSTIVE_SOLVER_HPP
#define ORG_VLEPROJECT_BARYONYX_EXHAUSTIVE_SOLVER_HPP

#include <limits>
#include <set>
#include <tuple>
#include <vector>

#include "itm-common.hpp"
#include "private.hpp"
#include "utils.hpp"

#include <cstdint>

namespace baryonyx {
namespace itm {

template<typename Mode, typename Float>
struct exhaustive_solver
{
    struct item
    {
        Float r = 0.0;
        int factor = 0;
        int variable = 0;
        int result = 0;

        item() noexcept = default;

        constexpr item(Float r_, int factor_, int variable_) noexcept
          : r(r_)
          , factor(factor_)
          , variable(variable_)
        {}
    };

    struct constraint
    {
        constraint() noexcept = default;

        constraint(const int k_) noexcept
          : k(k_)
        {}

        constraint(const int k_,
                   const int start_,
                   const int solutions_,
                   const int bk_min_,
                   const int bk_max_) noexcept
          : k(k_)
          , start(start_)
          , solutions(solutions_)
          , bk_min(bk_min_)
          , bk_max(bk_max_)
        {}

        int k;
        int start;
        int solutions;
        int bk_min = 0;
        int bk_max = 0;

        bool operator==(const constraint& other) const noexcept
        {
            return k == other.k;
        }

        bool operator<(const constraint& other) const noexcept
        {
            return k < other.k;
        }
    };

    std::vector<item> items;
    std::vector<int> flat_constraints;
    std::vector<std::int8_t> walkers;
    std::set<constraint> constraints;

    void reserve(std::size_t max_variables, std::size_t max_z_constraints)
    {
        bx_assert(max_variables > 0);
        bx_assert(max_z_constraints > 0);

        items.reserve(max_variables);
        flat_constraints.reserve(max_variables * max_z_constraints);
        walkers.reserve(max_variables);
    }

    template<typename C>
    void build_constraints(int k,
                           const C& constraint_elements,
                           int bk_min,
                           int bk_max)
    {
        const auto constraint_size = length(constraint_elements);
        const auto start_it_flat_constraints = length(flat_constraints);
        auto nb_solution = 0;

        bx_ensures(k >= 0);
        bx_ensures(constraint_size > 0);

        walkers.resize(constraint_elements.size(), 0);
        walkers.back() = 1;
        bool end = false;

        do {
            auto i = length(walkers) - 1;

            do {
                int sum = 0;
                for (int j = 0; j != constraint_size; ++j)
                    if (walkers[j])
                        sum += constraint_elements[j].factor;

                if (bk_min <= sum && sum <= bk_max) {
                    ++nb_solution;
                    for (int j = 0; j != constraint_size; ++j) {
                        flat_constraints.emplace_back(
                          walkers[j] ? constraint_elements[j].factor : 0);
                    }
                }

                ++walkers[i];

                if (walkers[i] > 1) {
                    walkers[i] = 0;

                    if (i == 0) {
                        end = true;
                        break;
                    } else {
                        --i;
                    }
                } else {
                    break;
                }
            } while (!end);
        } while (!end);

        bx_expects(start_it_flat_constraints >= 0);
        bx_expects(nb_solution > 0);

        constraints.emplace(
          k, start_it_flat_constraints, nb_solution, bk_min, bk_max);
    }

    static Float init_z() noexcept
    {
        if constexpr (std::is_same_v<Mode, minimize_tag>)
            return +std::numeric_limits<Float>::infinity();
        else
            return -std::numeric_limits<Float>::infinity();
    }

    bool is_best_solution(Float z, Float old_z) const noexcept
    {
        if constexpr (std::is_same_v<Mode, minimize_tag>)
            return z < old_z;
        else
            return z > old_z;
    }

    template<typename R>
    int solve(int k, R& reduced_cost, int r_size)
    {
        const auto it_constraint = constraints.find(k);
        bx_expects(it_constraint != constraints.end());

        items.resize(r_size);
        for (int i = 0; i != r_size; ++i) {
            items[i].r = reduced_cost[i].value;
            items[i].variable = reduced_cost[i].id;
            items[i].factor = reduced_cost[i].f;
            items[i].result = 0;
        }

        Float z_best = 0;
        auto best = 0;
        auto start_solution = it_constraint->start;

        for (int j = 0; j != r_size; ++j)
            if (flat_constraints[start_solution + j])
                z_best += reduced_cost[j].value;

        for (auto i = 1; i != it_constraint->solutions; ++i) {
            start_solution = it_constraint->start + (i * r_size);

            Float z = 0;
            for (int j = 0; j != r_size; ++j)
                if (flat_constraints[start_solution + j])
                    z += reduced_cost[j].value;

            if (is_best_solution(z, z_best)) {
                z_best = z;
                best = i;
            }
        }

        start_solution = it_constraint->start + (best * r_size);
        for (int i = 0; i != r_size; ++i)
            items[i].result = flat_constraints[start_solution + i] ? 1 : 0;

        std::sort(std::begin(items),
                  std::end(items),
                  [](const auto& lhs, const auto& rhs) {
                      if (lhs.result == rhs.result) {
                          if constexpr (std::is_same_v<Mode, minimize_tag>)
                              return lhs.r < rhs.r;
                          else
                              return lhs.r > rhs.r;
                      } else
                          return lhs.result > rhs.result;
                  });

        auto middle =
          std::find_if(std::begin(items),
                       std::end(items),
                       [](const auto& item) { return item.result == 0; });

        for (std::size_t i = 0, e = items.size(); i != e; ++i) {
            reduced_cost[i].value = items[i].r;
            reduced_cost[i].id = items[i].variable;
            reduced_cost[i].f = items[i].factor;
        }

        if (middle == std::end(items))
            return items[0].result == 0 ? -1 : r_size;

        return static_cast<int>(std::distance(std::begin(items), middle) - 1);
    }
};

} // namespace itm
} // namespace baryonyx

#endif // ORG_VLEPROJECT_BARYONYX_EXHAUSTIVE_SOLVER_HPP
