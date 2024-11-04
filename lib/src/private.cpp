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

    if (params.kappa_min < params.kappa_max && params.kappa_min >= 0 &&
        params.kappa_min < 1 && params.kappa_max <= 1 &&
        params.kappa_step >= 0 && params.kappa_step < 1) {
        ctx->parameters.kappa_min = params.kappa_min;
        ctx->parameters.kappa_step = params.kappa_step;
        ctx->parameters.kappa_max = params.kappa_max;
    }

    if (params.alpha >= 0)
        ctx->parameters.alpha = params.alpha;

    if (std::isfinite(params.pushing_k_factor))
        ctx->parameters.pushing_k_factor = params.pushing_k_factor;

    if (std::isfinite(params.pushing_objective_amplifier))
        ctx->parameters.pushing_objective_amplifier =
          params.pushing_objective_amplifier;

    if (std::isfinite(params.init_crossover_bastert_insertion))
        ctx->parameters.init_crossover_bastert_insertion =
          std::clamp(params.init_crossover_bastert_insertion, 0.0, 1.0);

    if (std::isfinite(params.init_crossover_solution_selection_mean))
        ctx->parameters.init_crossover_solution_selection_mean =
          std::clamp(params.init_crossover_solution_selection_mean, 0.0, 1.0);

    if (std::isfinite(params.init_crossover_solution_selection_stddev))
        ctx->parameters.init_crossover_solution_selection_stddev = std::clamp(
          params.init_crossover_solution_selection_stddev, 0.0, 1.0);

    if (std::isfinite(params.init_mutation_variable_mean))
        ctx->parameters.init_mutation_variable_mean =
          std::clamp(params.init_mutation_variable_mean, 0.0, 1.0);

    if (std::isfinite(params.init_mutation_variable_stddev))
        ctx->parameters.init_mutation_variable_stddev =
          std::clamp(params.init_mutation_variable_stddev, 0.0, 1.0);

    if (std::isfinite(params.init_mutation_value_mean))
        ctx->parameters.init_mutation_value_mean =
          std::clamp(params.init_mutation_value_mean, 0.0, 1.0);

    if (std::isfinite(params.init_mutation_value_stddev))
        ctx->parameters.init_mutation_value_stddev =
          std::clamp(params.init_mutation_value_stddev, 0.0, 1.0);

    if (params.init_policy_random >= 0 && params.init_policy_random <= 1)
        ctx->parameters.init_policy_random = params.init_policy_random;

    ctx->parameters.seed = params.seed;
    ctx->parameters.thread = params.thread;

    ctx->parameters.limit = (params.limit <= 0)
                              ? std::numeric_limits<long int>::max()
                              : params.limit;

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

    ctx->parameters.pushes_limit =
      params.pushes_limit <= 0 ? 0 : params.pushes_limit;

    ctx->parameters.pushing_iteration_limit =
      params.pushing_iteration_limit <= 0 ? 0 : params.pushing_iteration_limit;

    if (params.init_population_size >= 5)
        ctx->parameters.init_population_size = params.init_population_size;

    ctx->parameters.init_kappa_improve_start =
      std::clamp(params.init_kappa_improve_start, 0.0, 1.0);

    ctx->parameters.init_kappa_improve_increase =
      std::clamp(params.init_kappa_improve_increase, 0.0, 1.0);

    ctx->parameters.init_kappa_improve_stop =
      std::clamp(params.init_kappa_improve_stop,
                 ctx->parameters.init_kappa_improve_start,
                 1.0);

    ctx->parameters.pre_order = params.pre_order;
    ctx->parameters.order = params.order;
    ctx->parameters.init_policy = params.init_policy;
    ctx->parameters.cost_norm = params.cost_norm;
    ctx->parameters.mode = params.mode;
    ctx->parameters.preprocessor = params.preprocessor;
    ctx->parameters.observer = params.observer;
    ctx->parameters.storage = params.storage;
    ctx->parameters.solver = params.solver;
    ctx->parameters.use_buffer_solver = params.use_buffer_solver;
    ctx->parameters.debug = params.debug;
}

solver_parameters
context_get_solver_parameters(const context_ptr& ctx)
{
    return ctx->parameters;
}

} // namespace baryonyx
