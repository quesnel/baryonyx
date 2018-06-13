/* Copyright (C) 2018 INRA
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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_OPTIMIZER_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_OPTIMIZER_HPP

#include <baryonyx/core>

namespace baryonyx {
namespace itm {

template<typename Solver, typename realT, typename modeT, typename randomT>
inline result
dispatch_optimizer_constraint_order_parameter(const context_ptr& ctx, const problem& pb)
{
    switch (ctx->parameters.order) {
    case solver_parameters::constraint_order::none:
        return optimize_problem<
          Solver<R, modeT, randomT>,
          floatingpointT,
          modeT,
          compute_none<realT, randomT>,
          randomT>(ctx, pb);

    case solver_parameters::constraint_order::reversing:
        return optimize_problem<
          Solver<R, modeT, randomT>,
          floatingpointT,
          modeT,
          compute_reversing<realT, randomT>,
          randomT>(ctx, pb);

    case solver_parameters::constraint_order::random_sorting:
        return optimize_problem<
          Solver<R, modeT, randomT>,
          floatingpointT,
          modeT,
          compute_random<realT, randomT>,
          randomT>(ctx, pb);

    case solver_parameters::constraint_order::infeasibility_decr:
        return optimize_problem<
          Solver<R, modeT, randomT>,
          floatingpointT,
          modeT,
          compute_infeasibility<realT, randomT, compute_infeasibility_decr>,
          randomT>(ctx, pb);

    case solver_parameters::constraint_order::infeasibility_incr:
        return optimize_problem<
          Solver<R, modeT, randomT>,
          floatingpointT,
          modeT,
          compute_infeasibility<realT, randomT, compute_infeasibility_incr>,
          randomT>(ctx, pb);
    }

    return result(result_status::internal_error);
}

template<typename Solver>
inline result
dispatch_optimizer_parameters(const baryonyx::context_ptr& ctx, const problem& pb)
{
    using random_type = std::default_random_engine;

    if (pb.type == baryonyx::objective_function_type::maximize) {
        switch (ctx->parameters.float_type) {
        case solver_parameters::floating_point_type::float_type:
            return dispatch_optimizer_constraint_order_parameter<Solver,
                                              float,
                                              maximize_tag,
                                              random_type>(ctx, pb);
        case solver_parameters::floating_point_type::double_type:
            return dispatch_optimizer_constraint_order_parameter<Solver,
                                              double,
                                              maximize_tag,
                                              random_type>(ctx, pb);
        case solver_parameters::floating_point_type::longdouble_type:
            return dispatch_optimizer_constraint_order_parameter<Solver,
                                              long double,
                                              maximize_tag,
                                              random_type>(ctx, pb);
        }
    } else {
        switch (ctx->parameters.float_type) {
        case solver_parameters::floating_point_type::float_type:
            return dispatch_optimizer_constraint_order_parameter<Solver,
                                              float,
                                              minimize_tag,
                                              random_type>(ctx, pb);
        case solver_parameters::floating_point_type::double_type:
            return dispatch_optimizer_constraint_order_parameter<Solver,
                                              double,
                                              minimize_tag,
                                              random_type>(ctx, pb);
        case solver_parameters::floating_point_type::longdouble_type:
            return dispatch_optimizer_constraint_order_parameter<Solver,
                                              long double,
                                              minimize_tag,
                                              random_type>(ctx, pb);
        }
    }

    return result(result_status::internal_error);
}

result
optimize_equalities_01coeff(const context_ptr& ctx, const problem& pb);

result
optimize_equalities_101coeff(const context_ptr& ctx, const problem& pb);

result
optimize_inequalities_01coeff(const context_ptr& ctx, const problem& pb);

result
optimize_inequalities_101coeff(const context_ptr& ctx, const problem& pb);

result
optimize_inequalities_101coeff_buffered(const context_ptr& ctx, const problem& pb);

result
optimize_inequalities_Zcoeff(const context_ptr& ctx, const problem& pb);

} // namespace itm

inline result
optimize(const context_ptr& ctx, const problem& pb)
{
    if (ctx->method == "buffered") {
        baryonyx::warning(ctx, "buffered method is experimental\n");
        return itm::optimize_inequalities_101coeff_buffered(ctx, pb);
    } else if (not ctx->method.empty()) {
        baryonyx::warning(ctx, "undefined method {}");
    }

    switch (pb.problem_type) {
    case problem_solver_type::equalities_01:
        return itm::optimize_equalities_01coeff(ctx, pb);

    case problem_solver_type::equalities_101:
        return itm::optimize_equalities_101coeff(ctx, pb);

    case problem_solver_type::equalities_Z:
        baryonyx::warning(ctx, "Z coefficient is experimental\n");
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

} // namespace baryonyx

#endif