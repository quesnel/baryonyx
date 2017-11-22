/* Copyright (C) 2017 INRA
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

#ifndef ORG_VLEPROJECT_BARYONYX_BRANCH_AND_BOUND_HPP
#define ORG_VLEPROJECT_BARYONYX_BRANCH_AND_BOUND_HPP

#include "fixed_array.hpp"
#include "private.hpp"
#include "utils.hpp"

#include <algorithm>
#include <limits>
#include <vector>

#include <cassert>

namespace baryonyx {
namespace details {

template<typename modeT, typename floatingpointT>
struct branch_and_bound_solver
{
    struct item
    {
        floatingpointT r = 0.0;
        int factor = 0;
        int variable = 0;
    };

    struct node
    {
        std::vector<int> variables;
        floatingpointT z = 0.0;
        floatingpointT sumr = 0.0;
        int sumfactor = 0;
        int level = 0;
    };

    fixed_array<item> items;
    int bound;

    branch_and_bound_solver(std::size_t size_, int bound_)
      : items(size_)
      , bound(bound_)
    {
    }

    template<typename iteratorT>
    static inline void sort(iteratorT begin, iteratorT end, maximize_tag)
    {
        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return lhs.value < rhs.value;
        });
    }

    template<typename iteratorT>
    static inline void sort(iteratorT begin, iteratorT end, minimize_tag)
    {
        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return rhs.value < lhs.value;
        });
    }

    static inline void sort_items(fixed_array<item>& c, maximize_tag) noexcept
    {
        std::sort(c.begin(), c.end(), [](const auto& lhs, const auto& rhs) {
            return (lhs.r / lhs.factor) < (rhs.r / rhs.factor);
        });
    }

    static inline void sort_items(fixed_array<item>& c, minimize_tag) noexcept
    {
        std::sort(c.begin(), c.end(), [](const auto& lhs, const auto& rhs) {
            return (rhs.r / rhs.factor) < (lhs.r / lhs.factor);
        });
    }

    static inline bool is_best(floatingpointT lhs,
                               floatingpointT rhs,
                               maximize_tag) noexcept
    {
        return lhs > rhs;
    }

    static inline bool is_best(floatingpointT lhs,
                               floatingpointT rhs,
                               minimize_tag) noexcept
    {
        return lhs < rhs;
    }

    floatingpointT bound_node(const node& u,
                              int W,
                              const fixed_array<item>& items) noexcept
    {
        if (u.sumfactor > W)
            return 0.0;

        auto n = length(items);
        auto bound = u.sumr;
        auto j = u.level + 1;
        auto factor = u.sumfactor;

        while ((j < n) and (factor + items[j].factor <= W)) {
            factor += items[j].factor;
            bound += items[j].r;
            ++j;
        }

        if (j < n)
            bound += (W - factor) * items[j].r / items[j].factor;

        return bound;
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
        assert(std::find_if(items.cbegin(),
                            items.cend(),
                            [](auto elem) { return elem.factor <= 0; }) ==
                 items.cend() &&
               "No negative factor.");

        node best;
        best.sumr = init(modeT());

        sort_items(items, modeT());

        std::vector<node> queue;
        queue.emplace_back();
        queue.back().z = bound_node(best, bound, items);

        while (not queue.empty()) {
            node u = std::move(queue.back());
            queue.pop_back();

            if ((std::size_t)u.level == items.size())
                continue;

            node next_added = u;
            next_added.variables.emplace_back(items[u.level].variable);
            next_added.level += 1;
            next_added.sumfactor += items[u.level].factor;
            next_added.sumr += items[u.level].r;

            if (is_best(u.z, best.sumr, modeT())) {
                if (next_added.sumfactor <= bound) {
                    if (is_best(next_added.sumr, best.sumr, modeT()))
                        best = next_added;

                    if (is_best(next_added.z, best.sumr, modeT()))
                        queue.emplace_back(next_added);
                }
            }

            node next_not_added = u;
            next_not_added.level += 1;
            next_not_added.z = bound_node(next_not_added, bound, items);

            if (is_best(next_not_added.z, best.sumr, modeT()))
                queue.emplace_back(next_not_added);
        }

        // Finally, according to the items and reduced_cost vectors, we sort
        // reduced_cost and returns the last selected variable.

        if (best.variables.size() == items.size()) // If all variables are
            return length(best.variables);         // selected, return size
                                                   // of one vector.

        // Otherwise, swap reduced cost elements that appear in best solution
        // with the no selected elements.

        int nb = length(best.variables);
        int i = 0;

        for (auto elem : best.variables) {
            auto it =
              std::find_if(reduced_cost.begin(),
                           reduced_cost.begin() + items.size(),
                           [elem](const auto& r) { return r.id == elem; });

            std::swap(*it, reduced_cost[i]);
            ++i;
        }

        branch_and_bound_solver::sort(reduced_cost.begin() + nb,
                                      reduced_cost.begin() + items.size(),
                                      modeT());

        return nb;
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
branch_and_bound_solver(const AP& ap,
                        R& reduced_cost,
                        iteratorT begin,
                        iteratorT end,
                        int bound)
{
    assert(bound >= 0 && "No negative bound");

    details::branch_and_bound_solver<modeT, floatingpointT> slv(
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
branch_and_bound_solver(R& reduced_cost,
                        iteratorT begin,
                        iteratorT end,
                        int bound)
{
    assert(bound >= 0 && "No negative bound");

    details::branch_and_bound_solver<modeT, floatingpointT> slv(
      std::distance(begin, end), bound);

    for (int i = 0; begin != end; ++begin, ++i) {
        slv.items[i].r = reduced_cost[i].value;
        slv.items[i].factor = std::abs(*begin);
        slv.items[i].variable = reduced_cost[i].id;
    }

    return slv.solve(reduced_cost);
}

} // namespace baryonyx

#endif // ORG_VLEPROJECT_BARYONYX_BRANCH_AND_BOUND_HPP
