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

#include "itm-solver-common.hpp"

namespace baryonyx {
namespace itm {

result
solve(const context_ptr& ctx, const problem& pb)
{
    if (ctx->method == "buffered") {
        warning(ctx, "buffered method is experimental\n");
        return itm::solve_inequalities_101coeff_buffered(ctx, pb);
    } else if (!ctx->method.empty()) {
        warning(ctx, "undefined method {}");
    }

    switch (pb.problem_type) {
    case problem_solver_type::equalities_01:
        return itm::solve_equalities_01coeff(ctx, pb);

    case problem_solver_type::equalities_101:
        return itm::solve_equalities_101coeff(ctx, pb);

    case problem_solver_type::equalities_Z:
        warning(ctx, "Z coefficient is experimental\n");
        return itm::solve_inequalities_Zcoeff(ctx, pb);

    case problem_solver_type::inequalities_01:
        return itm::solve_inequalities_01coeff(ctx, pb);

    case problem_solver_type::inequalities_101:
        return itm::solve_inequalities_101coeff(ctx, pb);

    case problem_solver_type::inequalities_Z:
        warning(ctx, "Z coefficient is experimental\n");
        return itm::solve_inequalities_Zcoeff(ctx, pb);
    }

    return result(result_status::internal_error);
}

result
optimize(const context_ptr& ctx, const problem& pb)
{
    if (ctx->method == "buffered") {
        warning(ctx, "buffered method is experimental\n");
        return itm::optimize_inequalities_101coeff_buffered(ctx, pb);
    } else if (!ctx->method.empty()) {
        warning(ctx, "undefined method {}");
    }

    switch (pb.problem_type) {
    case problem_solver_type::equalities_01:
        return itm::optimize_equalities_01coeff(ctx, pb);

    case problem_solver_type::equalities_101:
        return itm::optimize_equalities_101coeff(ctx, pb);

    case problem_solver_type::equalities_Z:
        warning(ctx, "Z coefficient is experimental\n");
        return itm::optimize_inequalities_Zcoeff(ctx, pb);

    case problem_solver_type::inequalities_01:
        return itm::optimize_inequalities_01coeff(ctx, pb);

    case problem_solver_type::inequalities_101:
        return itm::optimize_inequalities_101coeff(ctx, pb);

    case problem_solver_type::inequalities_Z:
        baryonyx::warning(ctx, "Z coefficient is experimental\n");
        return itm::optimize_inequalities_Zcoeff(ctx, pb);
    }

    return result(result_status::internal_error);
}

} // namespace itm
} // namespace baryonyx
