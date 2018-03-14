/* Copyright (C) 2016-2018 INRA
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
#include "itm-common.hpp"
#include "private.hpp"
#include "utils.hpp"

#include <algorithm>
#include <limits>
#include <queue>

#include <cassert>

namespace baryonyx {
namespace itm {
namespace details {

template<typename floatingpointT>
inline bool
stop_iterating(floatingpointT value, minimize_tag) noexcept
{
    return value > 0;
}

template<typename floatingpointT>
inline bool
stop_iterating(floatingpointT value, maximize_tag) noexcept
{
    return value < 0;
}

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
        node(floatingpointT _z,
             floatingpointT _sumr,
             int _sumfactor,
             int _level)
          : z(_z)
          , sumr(_sumr)
          , sumfactor(_sumfactor)
          , level(_level)
        {}

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
    {}

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

    static inline void sort_items(fixed_array<item>& c, maximize_tag) noexcept
    {
        std::sort(c.begin(), c.end(), [](const auto& lhs, const auto& rhs) {
            return (lhs.r / static_cast<floatingpointT>(lhs.factor)) >
                   (rhs.r / static_cast<floatingpointT>(rhs.factor));
        });
    }

    static inline void sort_items(fixed_array<item>& c, minimize_tag) noexcept
    {
        std::sort(c.begin(), c.end(), [](const auto& lhs, const auto& rhs) {
            return (rhs.r / static_cast<floatingpointT>(rhs.factor)) >
                   (lhs.r / static_cast<floatingpointT>(lhs.factor));
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

    floatingpointT bound_node(const node& u, int W) noexcept
    {
        if (u.sumfactor >= W)
            return init(modeT());

        auto bound_ret = u.sumr;
        auto sumf = u.sumfactor;

        int j = u.level + 1;
        int n = length(items);

        for (; j != n; ++j) {
            if (sumf + items[j].factor > W)
                break;

            bound_ret += items[j].r;
            sumf += items[j].factor;
        }

        if (j < n)
            bound_ret +=
              static_cast<floatingpointT>(W - sumf) *
              (items[j].r / static_cast<floatingpointT>(items[j].factor));

        return bound_ret;
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
        std::queue<node> queue;

        node best(0, 0, 0, -1);

        sort_items(items, modeT());

        queue.emplace(bound_node(best, bound), floatingpointT(0), 0, -1);
        best.sumr = init(modeT());

        while (not queue.empty()) {
            node u = queue.front();
            queue.pop();

            int next_level = u.level + 1;
            if (u.level == -1)
                next_level = 0;

            if (u.level == length(items) - 1)
                continue;

            if (is_best(u.z, best.sumr, modeT())) {
                node next_added = u;
                next_added.variables.emplace_back(items[next_level].variable);
                next_added.level = next_level;
                next_added.sumfactor += items[next_level].factor;
                next_added.sumr += items[next_level].r;

                if (next_added.sumfactor <= bound) {
                    if (is_best(next_added.sumr, best.sumr, modeT()))
                        best = next_added;

                    if (is_best(next_added.z, best.sumr, modeT()))
                        queue.push(next_added);
                }

                node next_not_added = u;
                next_not_added.level += 1;
                next_not_added.z = bound_node(next_not_added, bound);

                if (is_best(next_not_added.z, best.sumr, modeT()))
                    queue.push(next_not_added);
            }
        }

        // Finally, according to the items and reduced_cost vectors, we sort
        // reduced_cost and returns the last selected variable.

        if (best.variables.size() == items.size()) // If all variables are
            return length(best.variables);         // selected, return size
                                                   // of one vector.

        // Otherwise, swap reduced cost elements that appear in best solution
        // with the no selected elements.

        int selected = length(best.variables);
        int i = 0;
        const int e = length(items);

        for (auto elem : best.variables) {
            auto it =
              std::find_if(items.begin(), items.end(), [elem](const auto& r) {
                  return r.variable == elem;
              });

            std::swap(*it, items[i]);
            ++i;
        }

        sort_r_items(items.begin() + selected, items.end(), modeT());

        int sum = 0;
        for (i = 0; i != selected; ++i)
            sum += items[i].factor;

        if (i == 0 and sum > bound)
            return -1;

        for (; i != e; ++i) {
            sum += items[i].factor;

            if (sum > bound or stop_iterating(items[i].r, modeT()))
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
};

} // namespace details

/*
 * @note For @c iteratorT (iterator to factor), we convert to absolute value
 *     negative factor to avoid the use of @c invert_a in the AP matrix.
 */
template<typename modeT,
         typename floatingpointT,
         typename A,
         typename R,
         typename iteratorT>
int
branch_and_bound_solver(const A& a,
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
        slv.items[i].factor = std::abs(a[begin->value]);
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

} // namespace itm
} // namespace baryonyx

#endif // ORG_VLEPROJECT_BARYONYX_BRANCH_AND_BOUND_HPP
