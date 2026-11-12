/* Copyright (C) 2016-2026 INRAE
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
solve_inequalities_101_buffered(const context& ctx, const problem& pb);

result
optimize_inequalities_101_buffered(const context& ctx, const problem& pb);

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
    if (ctx.parameters.use_buffer_solver) {
        switch (pb.problem_type) {
        case problem_solver_type::equalities_01:
        case problem_solver_type::equalities_101:
        case problem_solver_type::inequalities_01:
        case problem_solver_type::inequalities_101:
            return solve_inequalities_101_buffered(ctx, pb);

        default:
            break;
        }
    }

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
/**
 * Optimizes with the genetic and restart optimizers tuned by successive
 * halving (`optimizer=tuned`): half of the time limit compares candidate
 * settings, the rest runs the best one.
 */
result
tuned_optimize(const context& ctx, const problem& pb);

/**
 * Optimizes by large neighbourhood search (`optimizer=lns`): from a first
 * restart solution, frees small neighbourhoods of the best solution and
 * solves the reduced problems.
 */
result
lns_optimize(const context& ctx, const problem& pb);

inline result
optimize(const context& ctx, const problem& pb)
{
    if (ctx.parameters.optimizer == solver_parameters::optimizer_type::tuned)
        return tuned_optimize(ctx, pb);

    if (ctx.parameters.optimizer == solver_parameters::optimizer_type::lns)
        return lns_optimize(ctx, pb);

    if (ctx.parameters.use_buffer_solver) {
        switch (pb.problem_type) {
        case problem_solver_type::equalities_01:
        case problem_solver_type::equalities_101:
        case problem_solver_type::inequalities_01:
        case problem_solver_type::inequalities_101:
            return optimize_inequalities_101_buffered(ctx, pb);

        default:
            break;
        }
    }

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

    return result(result_status::internal_error);
}

/**
 * @brief Split the problem with affected/no-affected variable.
 *
 * @details This optimization function tries to found variable with many
 *     changes during an @c update_row and split the problem in two distinct
 *     problem: variable affected to 0, variable affected to 1. Internally,
 * the problem uses @c optimize.
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
