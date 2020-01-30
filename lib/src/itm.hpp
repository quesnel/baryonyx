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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_HPP

#include "private.hpp"
#include "problem.hpp"

namespace baryonyx {
namespace itm {

result
solve_equalities_01(const context& ctx, const problem& pb);

result
optimize_equalities_01(const context& ctx, const problem& pb);

result
solve_equalities_101(const context& ctx, const problem& pb);

result
optimize_equalities_101(const context& ctx, const problem& pb);

result
solve_inequalities_01(const context& ctx, const problem& pb);

result
optimize_inequalities_01(const context& ctx, const problem& pb);

result
solve_inequalities_101(const context& ctx, const problem& pb);

result
optimize_inequalities_101(const context& ctx, const problem& pb);

result
solve_inequalities_Z(const context& ctx, const problem& pb);

result
optimize_inequalities_Z(const context& ctx, const problem& pb);

result
solve_random_inequalities_101(const context& ctx, const problem& pb);

result
optimize_random_inequalities_101(const context& ctx, const problem& pb);

result
solve_random_equalities_101(const context& ctx, const problem& pb);

result
optimize_random_equalities_101(const context& ctx, const problem& pb);

result
solve_random_inequalities_01(const context& ctx, const problem& pb);

result
optimize_random_inequalities_01(const context& ctx, const problem& pb);

result
solve_random_equalities_01(const context& ctx, const problem& pb);

result
optimize_random_equalities_01(const context& ctx, const problem& pb);

/**
 * @brief Select the correct solver and solve the problem.
 *
 * @param ctx Current context.
 * @param pb Preprocessed problem.
 * @return A representation of the result.
 * @throw @c baryonyx::solver_error
 */
inline result
solve(const context& ctx, const problem& pb)
{
    if (ctx.parameters.solver == solver_parameters::solver_type::bastert) {
        if (ctx.method == "buffered")
            return solve_equalities_01(ctx, pb);

        switch (pb.problem_type) {
        case problem_solver_type::equalities_01:
            return solve_equalities_01(ctx, pb);

        case problem_solver_type::equalities_101:
            return solve_equalities_101(ctx, pb);

        case problem_solver_type::equalities_Z:
            return solve_inequalities_Z(ctx, pb);

        case problem_solver_type::inequalities_01:
            return solve_inequalities_01(ctx, pb);

        case problem_solver_type::inequalities_101:
            return solve_inequalities_101(ctx, pb);

        case problem_solver_type::inequalities_Z:
            return solve_inequalities_Z(ctx, pb);
        }
    } else {
        switch (pb.problem_type) {
        case problem_solver_type::equalities_01:
            return optimize_random_equalities_01(ctx, pb);

        case problem_solver_type::equalities_101:
            return optimize_random_equalities_101(ctx, pb);

        case problem_solver_type::equalities_Z:
            return result(result_status::internal_error);

        case problem_solver_type::inequalities_01:
            return optimize_random_inequalities_01(ctx, pb);

        case problem_solver_type::inequalities_101:
            return optimize_random_inequalities_101(ctx, pb);

        case problem_solver_type::inequalities_Z:
            return result(result_status::internal_error);
        }

        return result(result_status::internal_error);
    }

    return result(result_status::internal_error);
}

/**
 * @brief Select the correct optimizez and solve the problem.
 *
 * @param ctx Current context.
 * @param pb Preprocessed problem.
 * @return A representation of the result.
 * @throw @c baryonyx::solver_error
 */
inline result
optimize(const context& ctx, const problem& pb)
{
    if (ctx.parameters.solver == solver_parameters::solver_type::bastert) {
        if (ctx.method == "buffered")
            return optimize_equalities_01(ctx, pb);

        switch (pb.problem_type) {
        case problem_solver_type::equalities_01:
            return optimize_equalities_01(ctx, pb);

        case problem_solver_type::equalities_101:
            return optimize_equalities_101(ctx, pb);

        case problem_solver_type::equalities_Z:
            return optimize_inequalities_Z(ctx, pb);

        case problem_solver_type::inequalities_01:
            return optimize_inequalities_01(ctx, pb);

        case problem_solver_type::inequalities_101:
            return optimize_inequalities_101(ctx, pb);

        case problem_solver_type::inequalities_Z:
            return optimize_inequalities_Z(ctx, pb);
        }
    } else {
        switch (pb.problem_type) {
        case problem_solver_type::equalities_01:
            return optimize_random_equalities_01(ctx, pb);

        case problem_solver_type::equalities_101:
            return optimize_random_equalities_101(ctx, pb);

        case problem_solver_type::equalities_Z:
            return result(result_status::internal_error);

        case problem_solver_type::inequalities_01:
            return optimize_random_inequalities_01(ctx, pb);

        case problem_solver_type::inequalities_101:
            return optimize_random_inequalities_101(ctx, pb);

        case problem_solver_type::inequalities_Z:
            return result(result_status::internal_error);
        }
    }
    return result(result_status::internal_error);
}

/**
 * @brief auto-tune solver parameters to optimize the problem.
 *
 * @details The solver parameters are manually tuned to perform optimization
 *     and to found the best parameters to found the best solution.
 *
 * @param ctx A context.
 * @param pb Problem definition.
 * @return A representation of the result.
 *
 * @throw @c baryonyx::solver_error
 */
result
manual_optimize(const context& ctx, const problem& pb);

/**
 * @brief auto-tune solver parameters to optimize the problem using nlopt.
 *
 * @details The solver parameters are tuned to perform optimization and to
 *     found the best parameters to found the best solution. This function use
 *     the @c nlopt library to search best parameters.
 *
 * @note Available if baryonyx is build with nlopt support otherwise, this
 *     function rely to @c manual_optimize.
 *
 * @param ctx A context.
 * @param pb Problem definition.
 * @return A representation of the result.
 *
 * @throw @c baryonyx::solver_error
 */
result
nlopt_optimize(const context& ctx, const problem& pb);

/**
 * @brief Split the problem with affected/no-affected variable.
 *
 * @details This optimization function tries to found variable with many
 *     changes during an @c update_row and split the problem in two distinct
 *     problem: variable affected to 0, variable affected to 1. Internally, the
 *     problem can use @c optimize, @c manual_optimize or @c nlopt_optimize.
 *
 * @param ctx A context.
 * @param pb Problem definition.
 * @return A representation of the result.
 *
 * @throw @c baryonyx::solver_error
 */
result
branch_optimize(const context& ctx, const problem& pb);

} // namespace itm
} // namespace baryonyx

#endif
