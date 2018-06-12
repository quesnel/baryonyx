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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_COMMON_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_COMMON_HPP

#include <baryonyx/core>

#include <tuple>

#include <cassert>

namespace baryonyx {
namespace itm {

/**
 * x_type is a std::vector<bool> instead of a baryonyx::fixed_array<bool> to
 * use the specialized version of vector, which is used for elements of type
 * bool and optimizes for space.
 */
using x_type = std::vector<bool>;

struct maximize_tag
{};

struct minimize_tag
{};

inline double
best_solution_value(const baryonyx::result& res) noexcept
{
    assert(res.status == baryonyx::result_status::success);
    assert(not res.solutions.empty());

    return res.solutions.back().value;
}

inline bool
is_better_solution(const baryonyx::result& lhs,
                   const baryonyx::result& rhs,
                   maximize_tag) noexcept
{
    return best_solution_value(lhs) > best_solution_value(rhs);
}

inline bool
is_better_solution(const baryonyx::result& lhs,
                   const baryonyx::result& rhs,
                   minimize_tag) noexcept
{
    return best_solution_value(lhs) < best_solution_value(rhs);
}

struct merged_constraint
{
    merged_constraint(std::vector<function_element> elements_,
                      int min_,
                      int max_,
                      int id_)
      : elements(std::move(elements_))
      , min(min_)
      , max(max_)
      , id(id_)
    {}

    std::vector<function_element> elements;
    int min;
    int max;
    int id;
};

std::vector<merged_constraint>
make_merged_constraints(const context_ptr& ctx, const problem& pb);

baryonyx::result
solve(const baryonyx::context_ptr& ctx, const baryonyx::problem& pb);

baryonyx::result
optimize(const baryonyx::context_ptr& ctx,
         const baryonyx::problem& pb,
         int thread);

/**
 * @brief Auto-tune baryonyx solver parameters to optimize the problem.
 *
 * @details The @c nlopt library or a manual test is used to found the best
 *     parameters @c baryonyx::optimize function for the specified problem.
 *
 * @param ctx A context.
 * @param pb Problem definition.
 * @param thread Number of thread available for running optimization.
 * @return A representation of the result.
 * @throw @c baryonyx::solver_error
 */
baryonyx::result
automatic_optimizer(const baryonyx::context_ptr& ctx,
                    const baryonyx::problem& pb,
                    int thread);

} // namespace itm
} // namespace baryonyx

#endif
