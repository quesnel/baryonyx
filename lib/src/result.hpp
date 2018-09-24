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

#ifndef ORG_VLEPROJECT_BARYONYX_LIB_PRIVATE_RESULT_HPP
#define ORG_VLEPROJECT_BARYONYX_LIB_PRIVATE_RESULT_HPP

#include <baryonyx/core>

#include "itm-common.hpp"

#include <iosfwd>

namespace baryonyx {
namespace detail {

template<typename Mode>
inline bool
store_one_solution(const context_ptr& ctx,
                   result& res,
                   const std::vector<bool>& solution,
                   double value)
{
    if (itm::is_better_solution(value, res.solutions.back().value, Mode())) {
        res.solutions.back() = { solution, value };
        return true;
    }

    return false;
}

template<typename Mode>
inline bool
store_bound_solutions(const context_ptr& ctx,
                      result& res,
                      const std::vector<bool>& solution,
                      double value)
{
    if (res.solutions.size() == static_cast<size_t>(1)) {
        if (itm::is_better_solution(
              value, res.solutions.back().value, Mode())) {
            res.solutions.emplace_back(solution, value);
        } else {
            res.solutions.emplace(res.solutions.begin(), solution, value);
        }
    }

    if (itm::is_better_solution(value, res.solutions.back().value, Mode())) {
        res.solutions.back() = { solution, value };
    } else if (itm::is_better_solution(
                 res.solutions.front().value, value, Mode())) {
        res.solutions.front() = { solution, value };
    }

    return res.solutions.back().value == value;
}

template<typename Mode>
inline bool
store_five_solutions(const context_ptr& ctx,
                     result& res,
                     const std::vector<bool>& solution,
                     double value)
{
    auto it = res.solutions.rbegin();
    auto et = res.solutions.rend();

    if (res.solutions.size() == static_cast<size_t>(5)) {
        for (; it != et; ++it) {

            // If the current element is better, we need to shift all
            // previous solutions, drop the worst and replace the
            // solution.

            if (itm::is_better_solution(value, it->value, Mode())) {
                auto found = it.base();
                auto first = res.solutions.begin() + 1;

                for (; first != found; ++first)
                    std::iter_swap(first - 1, first);

                *found = { solution, value };
            }
        }
    } else {
        for (; it != et; ++it)
            if (itm::is_better_solution(value, it->value, Mode()))
                res.solutions.emplace(it.base(), solution, value);
    }

    return res.solutions.back().value == value;
}

} // namespace detail

template<typename Mode>
bool
store_solution(const context_ptr& ctx,
               result& res,
               const std::vector<bool>& solution,
               double value)
{
    // If the result solutions vector is empty, the solution and value
    // are pushed into the solution vector. All detail functions are
    // sure to have at lest one solution is the result solutions
    // vector.

    if (res.solutions.empty()) {
        res.solutions.emplace_back(solution, value);
        res.remaining_constraints = 0;
        res.status = result_status::success;
        return true;
    } else {
        switch (ctx->parameters.storage) {
        case solver_parameters::storage_type::one:
            return detail::store_one_solution<Mode>(ctx, res, solution, value);
        case solver_parameters::storage_type::bound:
            return detail::store_bound_solutions<Mode>(
              ctx, res, solution, value);
        case solver_parameters::storage_type::five:
            return detail::store_five_solutions<Mode>(
              ctx, res, solution, value);
        }
    }

    return detail::store_one_solution<Mode>(ctx, res, solution, value);
}

inline bool
store_advance(result& res, int constraint_remaining)
{
    if (res.remaining_constraints > constraint_remaining) {
        res.remaining_constraints = constraint_remaining;
        return true;
    }

    return false;
}

std::istream&
operator>>(std::istream& is, result& r);

struct best_solution_writer
{
    const result& res;

    best_solution_writer(const result& res_)
      : res(res_)
    {}
};

std::ostream&
operator<<(std::ostream& os, const best_solution_writer& writer);

} // namespace baryonyx

#endif
