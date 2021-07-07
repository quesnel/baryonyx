/* Copyright (C) 2018-2021 INRAE
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

#include "debug.hpp"
#include "itm.hpp"
#include "private.hpp"
#include "problem.hpp"
#include "utils.hpp"

#ifdef BARYONYX_HAVE_NLOPT
#include <nlopt.h>
#include <utility>

enum param
{
    param_theta = 0,
    param_delta,
    param_kappa_min,
    param_kappa_step,
    param_init_random,
    param_init_policy_random
};

struct nlopt_data
{
    nlopt_data(baryonyx::context& ctx_, baryonyx::problem pb_)
      : ctx(ctx_)
      , pb(std::move(pb_))
    {}

    baryonyx::context& ctx;
    const baryonyx::problem pb;
};

static double
nlopt_optimize_fun(unsigned /*n*/,
                   const double* x,
                   double* /*grad*/,
                   void* data_orig) noexcept
{
    try {
        auto* data = reinterpret_cast<nlopt_data*>(data_orig);

        data->ctx.parameters.theta = x[static_cast<int>(param_theta)];
        data->ctx.parameters.delta = x[static_cast<int>(param_delta)];
        data->ctx.parameters.kappa_min = x[static_cast<int>(param_kappa_min)];
        data->ctx.parameters.kappa_step =
          x[static_cast<int>(param_kappa_step)];
        data->ctx.parameters.init_policy_random =
          x[static_cast<int>(param_init_policy_random)];

        auto ret = baryonyx::itm::optimize(data->ctx, data->pb);
        if (!(ret))
            return HUGE_VAL;

        baryonyx::notice(data->ctx,
                         "theta: {} delta: {} kappa_min: {} kappa_step: {} "
                         "init_policy_random: {}: {:f}\n",
                         data->ctx.parameters.theta,
                         data->ctx.parameters.delta,
                         data->ctx.parameters.kappa_min,
                         data->ctx.parameters.kappa_step,
                         data->ctx.parameters.init_policy_random,
                         ret.solutions.back().value);

        return ret.solutions.back().value;
    } catch (const std::exception& e) {
        fmt::print("Exception in nlopt_optimize_fun: {}\n", e.what());
    }
    return HUGE_VAL;
}

static baryonyx::result
optimize(const baryonyx::context& ctx, const baryonyx::problem& pb)
{
    baryonyx::context internal(ctx);
    internal.log_priority = baryonyx::context::message_type::warning;

    nlopt_data data(internal, pb);

    double low[5] = { 0.0, 0.0001, 0.0, 1e-7, 0.0 };
    double up[5] = { 1.0, 0.1, 0.5, 0.01, 1.0 };
    double x[5] = { 0.5, 0.001, 0.1, 0.001, 0.5 };

    nlopt_opt opt = nlopt_create(NLOPT_LN_NELDERMEAD, 5u);

    nlopt_set_maxtime(opt, 3600.0);
    nlopt_set_vector_storage(opt, 100u);
    nlopt_set_lower_bounds(opt, low);
    nlopt_set_upper_bounds(opt, up);
    nlopt_set_min_objective(
      opt, nlopt_optimize_fun, reinterpret_cast<void*>(&data));

    // Seems interesting to loop over:
    // - more time-limit
    // - more loop-limit
    // - more optimization algorithms.

    double value;
    auto result = nlopt_optimize(opt, x, &value);

    nlopt_destroy(opt);

    if (result >= 1 || result == -4) {
        baryonyx::notice(
          ctx,
          "  - nlopt optimization found solution {:f}: with theta:{} "
          "delta:{} kappa-min:{} kappa-step:{} init-random:{}\n",
          value,
          x[param_theta],
          x[param_delta],
          x[param_kappa_min],
          x[param_kappa_step],
          x[param_init_policy_random]);

        internal.parameters.theta = x[param_theta];
        internal.parameters.delta = x[param_delta];
        internal.parameters.kappa_min = x[param_kappa_min];
        internal.parameters.kappa_step = x[param_kappa_step];
        internal.parameters.init_policy_random = x[param_init_policy_random];

        return baryonyx::itm::optimize(internal, pb);
    } else {
        baryonyx::notice(ctx,
                         "  - nlopt optimization fail. Toggle to manual "
                         "parameters optimization.\n");

        return baryonyx::itm::manual_optimize(ctx, pb);
    }
}
#endif

namespace baryonyx {
namespace itm {

#ifdef BARYONYX_HAVE_NLOPT
result
nlopt_optimize(const context& ctx, const problem& pb)
{
    notice(ctx, "- auto-tune parameters (nlopt) starts\n");

    return ::optimize(ctx, pb);
}
#else
result
nlopt_optimize(const context& ctx, const problem& pb)
{
    notice(ctx, "- auto-tune parameters (nlopt) starts\n");

    error(ctx,
          "  Baryonyx does not have nlopt. Toggle to manual parameters "
          "optimization.\n");

    return manual_optimize(ctx, pb);
}
#endif

} // namespace itm
} // namespace baryonyx
