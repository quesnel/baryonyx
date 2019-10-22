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

#include "private.hpp"
#include "utils.hpp"

namespace baryonyx {

const fmt::text_style context::message_style[] = {
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
context_set_parameters(const context_ptr& ctx,
                       const std::string& name,
                       std::string value)
{
    if (name == "method")
        ctx->method = std::move(value);
    else
        warning(ctx, "context: unknown variable {}.\n", name);
}

void
context_set_solver_parameters(const context_ptr& ctx,
                              const solver_parameters& params)
{
    if (params.time_limit > 0)
        ctx->parameters.time_limit = params.time_limit;
    else
        ctx->parameters.time_limit = -1;

    if (params.theta >= 0 && params.theta <= 1)
        ctx->parameters.theta = params.theta;

    if (params.delta == -1 || params.delta >= 0)
        ctx->parameters.delta = params.delta;

    if (params.kappa_min < params.kappa_max && params.kappa_min >= 0 &&
        params.kappa_min < 1 && params.kappa_max <= 1 &&
        params.kappa_step >= 0 && params.kappa_step < 1) {
        ctx->parameters.kappa_min = params.kappa_min;
        ctx->parameters.kappa_step = params.kappa_step;
        ctx->parameters.kappa_max = params.kappa_max;
    }

    if (params.alpha >= 0)
        ctx->parameters.alpha = params.alpha;

    if (std::isnormal(params.pushing_k_factor))
        ctx->parameters.pushing_k_factor = params.pushing_k_factor;

    if (std::isnormal(params.pushing_objective_amplifier))
        ctx->parameters.pushing_objective_amplifier =
          params.pushing_objective_amplifier;

    if (params.init_random >= 0 && params.init_random <= 1)
        ctx->parameters.init_random = params.init_random;

    if (params.init_policy_random >= 0 && params.init_policy_random <= 1)
        ctx->parameters.init_policy_random = params.init_policy_random;

    ctx->parameters.seed = params.seed;
    ctx->parameters.thread = params.thread;

    if (params.limit >= -1)
        ctx->parameters.limit = params.limit;
    else
        ctx->parameters.limit = -1;

    if (params.print_level >= 0)
        ctx->parameters.print_level = params.print_level;

    // Convert the value [0..1] into a percentage of loop or use directly the
    // number if the value.

    if (params.w < 0) {
        ctx->parameters.w = 0;
    } else {
        if (params.limit > 0) {
            auto limit = static_cast<double>(params.limit);

            if (params.w < 1)
                ctx->parameters.w = limit * params.w;
            else
                ctx->parameters.w = std::min(limit, params.w);
        } else {
            ctx->parameters.w = std::floor(params.w);
        }
    }

    if (params.pushes_limit >= 0)
        ctx->parameters.pushes_limit = params.pushes_limit;

    if (params.pushing_iteration_limit >= 0)
        ctx->parameters.pushing_iteration_limit =
          params.pushing_iteration_limit;

    ctx->parameters.pre_order = params.pre_order;
    ctx->parameters.order = params.order;
    ctx->parameters.float_type = params.float_type;
    ctx->parameters.init_policy = params.init_policy;
    ctx->parameters.cost_norm = params.cost_norm;
    ctx->parameters.mode = params.mode;
    ctx->parameters.preprocessor = params.preprocessor;
    ctx->parameters.observer = params.observer;
    ctx->parameters.storage = params.storage;
    ctx->parameters.solver = params.solver;
}

solver_parameters
context_get_solver_parameters(const context_ptr& ctx)
{
    return ctx->parameters;
}

} // namespace baryonyx
