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
optimize(const baryonyx::context_ptr& ctx, const baryonyx::problem& pb)
{
    manual_course array;

    auto old_log_priority = ctx->log_priority;
    ctx->log_priority = baryonyx::context::message_type::notice;

    std::array<int, 5> best_params;
    double best = +HUGE_VAL;

    do {
        ctx->parameters.theta = array.theta[array.it[0]];
        ctx->parameters.delta = array.delta[array.it[1]];
        ctx->parameters.kappa_min = array.kappa_min[array.it[2]];
        ctx->parameters.kappa_step = array.kappa_step[array.it[3]];
        ctx->parameters.init_random = array.init_random[array.it[4]];

        auto ret = baryonyx::itm::optimize(ctx, pb);
        if (ret) {
            if (best > ret.solutions.back().value) {
                best = ret.solutions.back().value;
                best_params = array.it;
            }
        }
    } while (array.next());

    baryonyx::notice(
      ctx,
      "  - manual optimization found solution {:f}: with theta:{} "
      "delta:{} kappa-min:{} kappa-step:{} init-random:{}\n",
      best,
      array.theta[array.it[0]],
      array.delta[array.it[1]],
      array.kappa_min[array.it[2]],
      array.kappa_step[array.it[3]],
      array.init_random[array.it[4]]);

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
