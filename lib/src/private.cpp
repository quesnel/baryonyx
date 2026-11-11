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

#include "private.hpp"
#include "utils.hpp"

namespace baryonyx {

constexpr fmt::text_style context::message_style[] = {
    (fmt::emphasis::bold | fmt::fg(fmt::terminal_color::red)), // emerg
    (fmt::emphasis::bold | fmt::fg(fmt::terminal_color::red)), // alert
    (fmt::emphasis::bold | fmt::fg(fmt::terminal_color::red)), // crit
    (fmt::fg(fmt::terminal_color::red)),                       // err
    (fmt::fg(fmt::terminal_color::yellow)),                    // warning
    (fmt::fg(fmt::terminal_color::white)),                     // notice
    (fmt::fg(fmt::terminal_color::white)),                     // info
    (fmt::fg(fmt::terminal_color::blue))                       // debug
};

void
context_set_solver_parameters(const context_ptr& ctx,
                              const solver_parameters& params)
{
    ctx->parameters.time_limit = params.time_limit <= 0
                                   ? std::numeric_limits<double>::infinity()
                                   : params.time_limit;

    if (params.theta >= 0 && params.theta <= 1)
        ctx->parameters.theta = params.theta;

    if (params.delta == -1 || params.delta >= 0)
        ctx->parameters.delta = params.delta;

    if (params.kappa_min >= 0 &&
        params.kappa_min < solver_parameters::kappa_max &&
        params.kappa_step >= 0 && params.kappa_step < 1) {
        ctx->parameters.kappa_min = params.kappa_min;
        ctx->parameters.kappa_step = params.kappa_step;
    }

    if (params.alpha >= 0)
        ctx->parameters.alpha = params.alpha;

    ctx->parameters.seed = params.seed;
    ctx->parameters.thread = params.thread;

    ctx->parameters.limit = (params.limit <= 0)
                              ? std::numeric_limits<long int>::max()
                              : params.limit;

    ctx->parameters.order = params.order;
    ctx->parameters.cost_norm = params.cost_norm;
    ctx->parameters.mode = params.mode;
    ctx->parameters.preprocessor = params.preprocessor;
    ctx->parameters.observer = params.observer;
    ctx->parameters.use_buffer_solver = params.use_buffer_solver;
    ctx->parameters.optimizer = params.optimizer;
    ctx->parameters.debug = params.debug;
}

solver_parameters
context_get_solver_parameters(const context_ptr& ctx)
{
    return ctx->parameters;
}

} // namespace baryonyx
