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

#include <array>

#include "private.hpp"
#include "utils.hpp"

#ifdef BARYONYX_HAVE_NLOPT
#include <nlopt.hpp>
#include <utility>
#endif

struct manual_course
{
    std::array<double, 5> theta = { { 0.0, 0.3, 0.5, 0.7, 1.0 } };
    std::array<double, 5> delta = { { -1, 1e-2, 0.05, 0.001, 0.0003 } };
    std::array<double, 5> kappa_min = { { 0, 1e-2, 0.05, 0.1, 0.3 } };
    std::array<double, 5> kappa_step = { { 1e-7, 1e-5, 1e-3, 1e-2, 1e-1 } };
    std::array<double, 5> init_random = { { 0.0, 0.3, 0.5, 0.7, 1.0 } };
    std::array<int, 5> it = { { 0, 0, 0, 0, 0 } };

    void reset()
    {
        std::fill(it.begin(), it.end(), 0);
    }

    bool next()
    {
        auto size = static_cast<int>(it.size() - 1);

        for (int i = size; i >= 0; --i) {
            if (it[i] + 1 > 4) {
                it[i] = 0;
            } else {
                ++it[i];
                return true;
            }
        }

        return false;
    }
};

static baryonyx::result
manual_optimize(const baryonyx::context_ptr& ctx, baryonyx::problem& pb)
{
    ctx->parameters.auto_tune =
      baryonyx::solver_parameters::auto_tune_parameters::disabled;

    manual_course array;

    std::array<int, 5> best_params;
    double best = +HUGE_VAL;

    do {
        ctx->parameters.theta = array.theta[array.it[0]];
        ctx->parameters.delta = array.delta[array.it[1]];
        ctx->parameters.kappa_min = array.kappa_min[array.it[2]];
        ctx->parameters.kappa_step = array.kappa_step[array.it[3]];
        ctx->parameters.init_random = array.init_random[array.it[4]];

        auto copy_pb = pb;
        auto ret = baryonyx::optimize(ctx, copy_pb);
        if (ret) {
            if (best > ret.solutions.back().value) {
                best = ret.solutions.back().value;
                best_params = array.it;
            }
        }
    } while (array.next());

    baryonyx::notice(
      ctx,
      "  - manual optimization found solution {}: with theta:{} "
      "delta:{} kappa-min:{} kappa-step:{} init-random:{}\n",
      best,
      array.theta[array.it[0]],
      array.delta[array.it[1]],
      array.kappa_min[array.it[2]],
      array.kappa_step[array.it[3]],
      array.init_random[array.it[4]]);

    return baryonyx::optimize(ctx, pb);
}

#ifdef BARYONYX_HAVE_NLOPT
enum param
{
    param_theta = 0,
    // param_delta,
    param_kappa_min,
    param_kappa_step,
    param_init_random
};

struct nlopt_data
{
    nlopt_data(const baryonyx::context_ptr& ctx_, baryonyx::problem  pb_)
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
        // data->ctx->parameters.delta = x[static_cast<int>(param_delta)];
        data->ctx->parameters.kappa_min = x[static_cast<int>(param_kappa_min)];
        data->ctx->parameters.kappa_step =
          x[static_cast<int>(param_kappa_step)];
        data->ctx->parameters.init_random =
          x[static_cast<int>(param_init_random)];

        auto copy_pb(data->pb);
        auto ret = baryonyx::optimize(data->ctx, copy_pb);
        if (not(ret))
            return HUGE_VAL;

        // baryonyx::notice(data->ctx,
        //                  "theta: {} delta: {} kappa_min: {} kappa_step: {} "
        //                  "init_random: {}: {}\n",
        //                  data->ctx->parameters.theta,
        //                  data->ctx->parameters.delta,
        //                  data->ctx->parameters.kappa_min,
        //                  data->ctx->parameters.kappa_step,
        //                  data->ctx->parameters.init_random,
        //                  ret.solutions.back().value);

        baryonyx::notice(data->ctx,
                         "theta: {} kappa_min: {} kappa_step: {} "
                         "init_random: {}: {}\n",
                         data->ctx->parameters.theta,
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
nlopt_optimize(const baryonyx::context_ptr& ctx, baryonyx::problem& pb)
{
    auto old_log_priority = ctx->log_priority;
    ctx->log_priority = baryonyx::context::message_type::notice;

    ctx->parameters.auto_tune =
      baryonyx::solver_parameters::auto_tune_parameters::disabled;
    nlopt_data data(ctx, pb);

    const std::vector<double> low{ 0, 0.0, 1e-7, 0 };
    const std::vector<double> up{ 1, 0.5, 0.01, 1 };
    std::vector<double> x{ 0.5, 0, 0.001, 0.5 };

    nlopt::opt opt(nlopt::LN_NELDERMEAD, 4);
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

    if (result >= 1 or result == -4) {
        baryonyx::notice(
          ctx,
          "  - nlopt optimization found solution {}: with theta:xxx "
          "delta:{} kappa-min:{} kappa-step:{} init-random:{}\n",
          value,
          x[param_theta],
          // x[param_delta],
          x[param_kappa_min],
          x[param_kappa_step],
          x[param_init_random]);

        ctx->parameters.theta = x[param_theta];
        // ctx->parameters.delta = x[param_delta];
        ctx->parameters.kappa_min = x[param_kappa_min];
        ctx->parameters.kappa_step = x[param_kappa_step];
        ctx->parameters.init_random = x[param_init_random];

        return baryonyx::optimize(ctx, pb);
    } else {
        baryonyx::notice(ctx,
                         "  - nlopt optimization fail. Toggle to manual "
                         "parameters optimization.\n");

        return manual_optimize(ctx, pb);
    }
}
#endif

namespace baryonyx {
namespace itm {

baryonyx::result
automatic_optimizer(const baryonyx::context_ptr& ctx,
                    baryonyx::problem& pb,
                    int thread)
{
    baryonyx::notice(ctx, "- Automatic optimization starts\n");

    assert(ctx->parameters.auto_tune !=
           baryonyx::solver_parameters::auto_tune_parameters::disabled);

    (void)thread;

#ifdef BARYONYX_HAVE_NLOPT
    if (ctx->parameters.auto_tune ==
        baryonyx::solver_parameters::auto_tune_parameters::manual)
        return ::manual_optimize(ctx, pb);
    else
        return ::nlopt_optimize(ctx, pb);
#else
    if (ctx->auto_tune == context::auto_tune_parameters::nlop)
        baryonyx::error(ctx,
                        "Baryonyx does not have nlopt. "
                        "Toggle to manual parameters optimization.\n");

    return ::manual_optimize(ctx, pb);
#endif
}
} // namespace itm
} // namespace baryonyx
