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

#include "itm-common.hpp"

#include "itm-solver-equalities-01.hpp"
#include "itm-solver-equalities-101.hpp"
#include "itm-solver-inequalities-01.hpp"
#include "itm-solver-inequalities-101-buffered.hpp"
#include "itm-solver-inequalities-101.hpp"
#include "itm-solver-inequalities-Z.hpp"

#include <set>

#include <cassert>

namespace baryonyx {
namespace itm {

template<typename floatingpointT, typename modeT, typename randomT>
struct solver_inequalities_Zcoeff;

template<typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT>
inline result
dispatch_solver(const context_ptr& ctx, const problem& pb)
{
    if (ctx->method == "buffered") {
        baryonyx::warning(ctx, "buffered method is experimental\n");

        return solve_problem<
          solver_inequalities_101coeff_buffered<floatingpointT,
                                                modeT,
                                                randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb);
    }

    switch (pb.problem_type) {
    case problem_solver_type::equalities_01:
        return solve_problem<
          solver_equalities_01coeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb);
    case problem_solver_type::equalities_101:
        return solve_problem<
          solver_equalities_101coeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb);
    case problem_solver_type::equalities_Z:
        baryonyx::warning(ctx, " \u26A0 Z coefficient is experimental\n");
        return solve_problem<
          solver_inequalities_Zcoeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb);
    case problem_solver_type::inequalities_01:
        return solve_problem<
          solver_inequalities_01coeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb);
    case problem_solver_type::inequalities_101:
        return solve_problem<
          solver_inequalities_101coeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb);
    case problem_solver_type::inequalities_Z:
        baryonyx::warning(ctx, " \u26A0 Z coefficient is experimental\n");
        return solve_problem<
          solver_inequalities_Zcoeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb);
    }

    return result(result_status::internal_error);
}

template<typename realT, typename modeT, typename randomT>
inline result
dispatch_solver_parameters(const context_ptr& ctx, const problem& pb)
{
    switch (ctx->parameters.order) {
    case solver_parameters::constraint_order::none:
        return dispatch_solver<realT,
                               modeT,
                               compute_none<realT, randomT>,
                               randomT>(ctx, pb);
    case solver_parameters::constraint_order::reversing:
        return dispatch_solver<realT,
                               modeT,
                               compute_reversing<realT, randomT>,
                               randomT>(ctx, pb);
    case solver_parameters::constraint_order::random_sorting:
        return dispatch_solver<realT,
                               modeT,
                               compute_random<realT, randomT>,
                               randomT>(ctx, pb);
    case solver_parameters::constraint_order::infeasibility_decr:
        return dispatch_solver<
          realT,
          modeT,
          compute_infeasibility<realT, randomT, compute_infeasibility_decr>,
          randomT>(ctx, pb);
    case solver_parameters::constraint_order::infeasibility_incr:
        return dispatch_solver<
          realT,
          modeT,
          compute_infeasibility<realT, randomT, compute_infeasibility_incr>,
          randomT>(ctx, pb);
    }

    return result(result_status::internal_error);
}

result
solve(const baryonyx::context_ptr& ctx, const problem& pb)
{
    info(ctx, "Solve mode\n");
    baryonyx::print(ctx);

    using random_type = std::default_random_engine;

    if (pb.type == baryonyx::objective_function_type::maximize) {
        switch (ctx->parameters.float_type) {
        case solver_parameters::floating_point_type::float_type:
            return dispatch_solver_parameters<float,
                                              maximize_tag,
                                              random_type>(ctx, pb);
        case solver_parameters::floating_point_type::double_type:
            return dispatch_solver_parameters<double,
                                              maximize_tag,
                                              random_type>(ctx, pb);
        case solver_parameters::floating_point_type::longdouble_type:
            return dispatch_solver_parameters<long double,
                                              maximize_tag,
                                              random_type>(ctx, pb);
        }
    } else {
        switch (ctx->parameters.float_type) {
        case solver_parameters::floating_point_type::float_type:
            return dispatch_solver_parameters<float,
                                              minimize_tag,
                                              random_type>(ctx, pb);
        case solver_parameters::floating_point_type::double_type:
            return dispatch_solver_parameters<double,
                                              minimize_tag,
                                              random_type>(ctx, pb);
        case solver_parameters::floating_point_type::longdouble_type:
            return dispatch_solver_parameters<long double,
                                              minimize_tag,
                                              random_type>(ctx, pb);
        }
    }

    return result(result_status::internal_error);
}

} // namespace itm
} // namespace baryonyx
