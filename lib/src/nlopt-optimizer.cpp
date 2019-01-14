/* Copyright (C) 2018-2019 INRA
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
#include <nlopt.hpp>
#include <utility>
#endif

#ifdef BARYONYX_HAVE_NLOPT
enum param
{
    param_theta = 0,
    param_delta,
    param_kappa_min,
    param_kappa_step,
    param_init_random
};

struct nlopt_data
{
    nlopt_data(const baryonyx::context_ptr& ctx_, baryonyx::problem pb_)
      : ctx(ctx_)
      , pb(std::move(pb_))
    {}

    const baryonyx::context_ptr& ctx;
    const baryonyx::problem pb;
};

static double
nlopt_optimize_fun(const std::vector<double>& x,
                   std::vector<double>& /*grad*/,
                   void* data_orig)
{
    try {
        auto* data = reinterpret_cast<nlopt_data*>(data_orig);

        data->ctx->parameters.theta = x[static_cast<int>(param_theta)];
        data->ctx->parameters.delta = x[static_cast<int>(param_delta)];
        data->ctx->parameters.kappa_min = x[static_cast<int>(param_kappa_min)];
        data->ctx->parameters.kappa_step =
          x[static_cast<int>(param_kappa_step)];
        data->ctx->parameters.init_random =
          x[static_cast<int>(param_init_random)];

        auto ret = baryonyx::itm::optimize(data->ctx, data->pb);
        if (!(ret))
            return HUGE_VAL;

        baryonyx::notice(data->ctx,
                         "theta: {} delta: {} kappa_min: {} kappa_step: {} "
                         "init_random: {}: {:f}\n",
                         data->ctx->parameters.theta,
                         data->ctx->parameters.delta,
                         data->ctx->parameters.kappa_min,
                         data->ctx->parameters.kappa_step,
                         data->ctx->parameters.init_random,
                         ret.solutions.back().value);

        return ret.solutions.back().value;
    } catch (const std::exception& e) {
        fmt::print("Exception in nlopt_optimize_fun: {}\n", e.what());
    }
    return HUGE_VAL;
}

static baryonyx::result
optimize(const baryonyx::context_ptr& ctx, const baryonyx::problem& pb)
{
    auto old_log_priority = ctx->log_priority;
    ctx->log_priority = baryonyx::context::message_type::notice;

    nlopt_data data(ctx, pb);

    const std::vector<double> low{ 0, 0.0001, 0.0, 1e-7, 0 };
    const std::vector<double> up{ 1, 0.1, 0.5, 0.01, 1 };
    std::vector<double> x{ 0.5, 0.001, 0.1, 0.001, 0.5 };

    nlopt::opt opt(nlopt::LN_NELDERMEAD, 5);
    opt.set_maxtime(3600);
    opt.set_vector_storage(100);
    opt.set_lower_bounds(low);
    opt.set_upper_bounds(up);
    opt.set_min_objective(nlopt_optimize_fun, reinterpret_cast<void*>(&data));

    //
    // Seems interesting to loop over:
    // - more time-limit
    // - more loop-limit
    // - more optimization algorithms.
    //

    double value;
    auto result = opt.optimize(x, value);

    ctx->log_priority = old_log_priority;

    if (result >= 1 || result == -4) {
        baryonyx::notice(
          ctx,
          "  - nlopt optimization found solution {:f}: with theta:xxx "
          "delta:{} kappa-min:{} kappa-step:{} init-random:{}\n",
          value,
          x[param_theta],
          x[param_delta],
          x[param_kappa_min],
          x[param_kappa_step],
          x[param_init_random]);

        ctx->parameters.theta = x[param_theta];
        ctx->parameters.delta = x[param_delta];
        ctx->parameters.kappa_min = x[param_kappa_min];
        ctx->parameters.kappa_step = x[param_kappa_step];
        ctx->parameters.init_random = x[param_init_random];

        return baryonyx::itm::optimize(ctx, pb);
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

result
nlopt_optimize(const context_ptr& ctx, const problem& pb)
{
    notice(ctx, "- auto-tune parameters (nlopt) starts\n");

#ifdef BARYONYX_HAVE_NLOPT
    return ::optimize(ctx, pb);
#else
    error(ctx,
          "  Baryonyx does not have nlopt. Toggle to manual parameters "
          "optimization.\n");

    return manual_optimize(ctx, pb);
#endif
}

} // namespace itm
} // namespace baryonyx
