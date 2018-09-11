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

#include "debug.hpp"
#include "itm.hpp"
#include "private.hpp"
#include "problem.hpp"
#include "utils.hpp"

struct manual_course
{
    std::unique_ptr<double[]> theta;
    std::unique_ptr<double[]> delta;
    std::unique_ptr<double[]> kappa_min;
    std::unique_ptr<double[]> kappa_step;
    std::unique_ptr<double[]> init_random;
    std::unique_ptr<int[]> iterators;

    const int length;

    manual_course(double theta_,
                  double delta_,
                  double kappa_min_,
                  double kappa_step_,
                  double init_random_,
                  int length_)
      : theta(std::make_unique<double[]>(length_))
      , delta(std::make_unique<double[]>(length_))
      , kappa_min(std::make_unique<double[]>(length_))
      , kappa_step(std::make_unique<double[]>(length_))
      , init_random(std::make_unique<double[]>(length_))
      , iterators(std::make_unique<int[]>(length_))
      , length(length_)
    {
        bx_ensures(length > 2);

        theta[0] = { 0.5 };
        theta[1] = 1.0 / static_cast<double>(length);

        delta[0] = { -1 };
        delta[1] = 0.1 / static_cast<double>(length);

        kappa_min[0] = { 0 };
        kappa_min[1] = 1e-2 / static_cast<double>(length);

        kappa_step[0] = { 1e-3 };
        kappa_step[1] = 1e-7 / static_cast<double>(length);

        init_random[0] = { 0.5 };
        init_random[1] = 0.9 / static_cast<double>(length);

        for (int i = 2; i != length; ++i)
            theta[i] = theta[i - 1] + 1.0 / static_cast<double>(length);

        for (int i = 2; i != length; ++i)
            delta[i] =
              theta[i - 1] + 0.1 / (5.0 * static_cast<double>(length));

        for (int i = 2; i != length; ++i)
            kappa_min[i] = theta[i - 1] + 1e-2 / static_cast<double>(length);

        for (int i = 2; i != length; ++i)
            kappa_step[i] =
              theta[i - 1] + 1e-7 / (5.0 * static_cast<double>(length));

        for (int i = 2; i != length; ++i)
            init_random[i] = theta[i - 1] + 0.9 / static_cast<double>(length);

        reset();
    }

    void reset()
    {
        std::fill(iterators.get(), iterators.get() + length, 0);
    }

    bool next()
    {
        for (int i = length - 1; i >= 0; --i) {
            if (iterators[i] + 1 >= length) {
                iterators[i] = 0;
            } else {
                ++iterators[i];
                return true;
            }
        }

        return false;
    }
};

static baryonyx::result
optimize(const baryonyx::context_ptr& ctx, const baryonyx::problem& pb)
{
    manual_course array(ctx->parameters.theta,
                        ctx->parameters.delta,
                        ctx->parameters.kappa_min,
                        ctx->parameters.kappa_step,
                        ctx->parameters.init_random,
                        5);

    auto old_log_priority = ctx->log_priority;
    ctx->log_priority = baryonyx::context::message_type::notice;

    auto best_params = std::make_unique<int[]>(array.length);
    double best = +HUGE_VAL;

    do {
        ctx->parameters.theta = array.theta[array.iterators[0]];
        ctx->parameters.delta = array.delta[array.iterators[1]];
        ctx->parameters.kappa_min = array.kappa_min[array.iterators[2]];
        ctx->parameters.kappa_step = array.kappa_step[array.iterators[3]];
        ctx->parameters.init_random = array.init_random[array.iterators[4]];

        baryonyx::notice(ctx,
                         "  - optimization with theta:{} delta:{} "
                         "kappa-min:{} kappa-step:{} init-random:{} ",
                         ctx->parameters.theta,
                         ctx->parameters.delta,
                         ctx->parameters.kappa_min,
                         ctx->parameters.kappa_step,
                         ctx->parameters.init_random);

        auto ret = baryonyx::itm::optimize(ctx, pb);
        if (ret) {
            baryonyx::notice(ctx, "{:f}\n", ret.solutions.back().value);
            if (best > ret.solutions.back().value) {
                best = ret.solutions.back().value;

                std::copy(array.iterators.get(),
                          array.iterators.get() + array.length,
                          best_params.get());
            }
        } else {
            baryonyx::notice(ctx, "no solution\n");
        }
    } while (array.next());

    baryonyx::notice(
      ctx,
      "  - manual optimization found solution {:f}: with theta:{} "
      "delta:{} kappa-min:{} kappa-step:{} init-random:{}\n",
      best,
      array.theta[array.iterators[0]],
      array.delta[array.iterators[1]],
      array.kappa_min[array.iterators[2]],
      array.kappa_step[array.iterators[3]],
      array.init_random[array.iterators[4]]);

    ctx->log_priority = old_log_priority;

    return baryonyx::itm::optimize(ctx, pb);
}

namespace baryonyx {
namespace itm {

result
manual_optimize(const baryonyx::context_ptr& ctx, const baryonyx::problem& pb)
{
    baryonyx::notice(ctx, "- auto-tune parameters (manual) starts\n");

    return ::optimize(ctx, pb);
}

} // namespace itm
} // namespace baryonyx
