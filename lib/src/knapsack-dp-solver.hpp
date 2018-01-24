/* Copyright (C) 2017 INRA
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublipnse, and/or sell copies of the Software, and to
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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_KNAPSACK_SOLVER_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_KNAPSACK_SOLVER_HPP

#include "fixed_2darray.hpp"
#include "private.hpp"

#include <vector>

#include <cassert>

namespace baryonyx {
namespace details {

template<typename modeT, typename floatingpointT>
struct knapsack_dp_solver
{
    struct item
    {
        floatingpointT r = 0.0;
        int factor = 0;
        int variable = 0;
    };

    fixed_array<item> items;
    fixed_2darray<floatingpointT> best;
    int capacity;

    template<typename iteratorT>
    static inline void sort_r_items(iteratorT begin,
                                    iteratorT end,
                                    maximize_tag)
    {
        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return rhs.r < lhs.r;
        });
    }

    template<typename iteratorT>
    static inline void sort_r_items(iteratorT begin,
                                    iteratorT end,
                                    minimize_tag)
    {
        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return lhs.r < rhs.r;
        });
    }

    knapsack_dp_solver(std::size_t size_, int capacity_)
      : items(size_)
      , best(size_ + 1, capacity_ + 1, 0)
      , capacity(capacity_)
    {
    }

    floatingpointT get_best(floatingpointT lhs,
                            floatingpointT rhs,
                            maximize_tag) noexcept
    {
        return std::max(lhs, rhs);
    }

    floatingpointT get_best(floatingpointT lhs,
                            floatingpointT rhs,
                            minimize_tag) noexcept
    {
        return std::min(lhs, rhs);
    }

    floatingpointT init(maximize_tag)
    {
        return -std::numeric_limits<floatingpointT>::infinity();
    }

    floatingpointT init(minimize_tag)
    {
        return +std::numeric_limits<floatingpointT>::infinity();
    }

    template<typename R>
    int solve(R& reduced_cost)
    {
        const std::size_t n = items.size();
        const std::size_t W = static_cast<size_t>(capacity);
        std::size_t i, j;

        for (j = 1; j <= W; ++j)
            best(0, j) = init(modeT());

        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= W; ++j) {
                if (items[i - 1].factor <= (int)j)
                    best(i, j) =
                      get_best(best(i - 1, j),
                               items[i - 1].r + best(i - 1,

                                                     j - items[i - 1].factor),
                               modeT());
                else
                    best(i, j) = best(i - 1, j);
            }
        }

        i = items.size();
        j = static_cast<std::size_t>(capacity);

        std::vector<int> variables;

        for (; i > 0 and j > 0; --i) {
            if (best(i, j) != best(i - 1, j)) {
                variables.emplace_back(items[i - 1].variable);
                j -= items[i - 1].factor;
            }
        }

        // Finally, according to the items and reduced_cost vectors, we sort
        // reduced_cost and returns the last selected variable.

        if (variables.size() == items.size()) // If all variables are
            return length(variables);         // selected, return size
                                              // of one vector.

        // Otherwise, swap reduced cost elements that appear in best solution
        // with the no selected elements.

        {
            int selected = length(variables);
            int i = 0;
            const int e = length(items);

            for (auto elem : variables) {
                auto it = std::find_if(
                  items.begin(), items.end(), [elem](const auto& r) {
                      return r.variable == elem;
                  });

                std::swap(*it, items[i]);
                ++i;
            }

            sort_r_items(items.begin() + selected, items.end(), modeT());

            int sum = 0;
            for (i = 0; i != selected; ++i)
                sum += items[i].factor;

            if (i == 0 and sum > capacity)
                return -1;

            for (; i != e; ++i) {
                sum += items[i].factor;

                if (sum > capacity or stop_iterating(items[i].r, modeT()))
                    break;
            }

            if (i == e)
                return length(items);

            selected = i;

            // Sort reduced cost according to the items vector order.

            for (i = 0; i != e; ++i) {
                int elem = items[i].variable;

                auto it =
                  std::find_if(reduced_cost.begin(),
                               reduced_cost.begin() + items.size(),
                               [elem](const auto& r) { return r.id == elem; });

                std::swap(reduced_cost[i], *it);
            }

            return selected;
        }
    }
};

} // namespace details

/*
 * @note For @c iteratorT (iterator to factor), we convert to absolute value
 *     negative factor to avoid the use of @c invert_a in the AP matrix.
 */
template<typename modeT,
         typename floatingpointT,
         typename AP,
         typename R,
         typename iteratorT>
int
knapsack_dp_solver(const AP& ap,
                   R& reduced_cost,
                   iteratorT begin,
                   iteratorT end,
                   int bound)
{
    assert(bound >= 0 && "No negative bound");

    details::knapsack_dp_solver<modeT, floatingpointT> slv(
      std::distance(begin, end), bound);

    for (int i = 0; begin != end; ++begin, ++i) {
        slv.items[i].r = reduced_cost[i].value;
        slv.items[i].factor = std::abs(ap.A()[begin->value]);
        slv.items[i].variable = reduced_cost[i].id;
    }

    return slv.solve(reduced_cost);
}

/*
 * @note For @c iteratorT (iterator to factor), we convert to absolute value
 *     negative factor to avoid the use of @c invert_a in the AP matrix.
 */
template<typename modeT,
         typename floatingpointT,
         typename R,
         typename iteratorT>
int
knapsack_dp_solver(R& reduced_cost, iteratorT begin, iteratorT end, int bound)
{
    assert(bound >= 0 && "No negative bound");

    details::knapsack_dp_solver<modeT, floatingpointT> slv(
      std::distance(begin, end), bound);

    for (int i = 0; begin != end; ++begin, ++i) {
        slv.items[i].factor = std::abs(*begin);
        slv.items[i].r = reduced_cost[i].value;
        slv.items[i].variable = reduced_cost[i].id;
    }

    return slv.solve(reduced_cost);
}

} // namespace baryonyx

#endif
