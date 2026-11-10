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

#include <cmath>
#include <functional>
#include <limits>
#include <string_view>
#include <vector>

#include "debug.hpp"
#include "itm.hpp"
#include "private.hpp"
#include "problem.hpp"
#include "utils.hpp"

namespace {

namespace bx = baryonyx;

/** A tuned parameter: its name, how to set it and the values to try. */
struct tuned_parameter
{
    std::string_view name;
    std::function<void(bx::solver_parameters&, double)> set;
    std::vector<double> values;
};

std::vector<tuned_parameter>
tuned_parameters()
{
    using optimizer_type = bx::solver_parameters::optimizer_type;

    return {
        { "optimizer",
          [](auto& p, double v) {
              p.optimizer =
                v == 0.0 ? optimizer_type::genetic : optimizer_type::restart;
          },
          { 0.0, 1.0 } },
        { "theta",
          [](auto& p, double v) { p.theta = v; },
          { 0.2, 0.4, 0.6, 0.8 } },
        { "delta",
          [](auto& p, double v) { p.delta = v; },
          { -1.0, 1e-3, 1e-2, 1e-1 } },
        { "kappa-min",
          [](auto& p, double v) { p.kappa_min = v; },
          { 0.0, 0.05, 0.1, 0.2 } },
        { "kappa-step",
          [](auto& p, double v) { p.kappa_step = v; },
          { 1e-4, 3e-4, 1e-3, 3e-3, 1e-2 } },
        { "init-policy-random",
          [](auto& p, double v) { p.init_policy_random = v; },
          { 0.0, 0.25, 0.5, 0.75 } },
    };
}

// A feasible solution beats an infeasible one, then the objective decides,
// or the number of violated constraints between infeasible ones.
bool
is_better(const bx::result& lhs,
          const bx::result& rhs,
          bx::objective_function_type type)
{
    if (rhs.solutions.empty())
        return !lhs.solutions.empty();

    if (lhs.solutions.empty())
        return false;

    const bool lhs_ok = lhs.status == bx::result_status::success;
    const bool rhs_ok = rhs.status == bx::result_status::success;

    if (lhs_ok != rhs_ok)
        return lhs_ok;

    if (!lhs_ok)
        return lhs.remaining_constraints < rhs.remaining_constraints;

    const auto l = lhs.solutions.back().value;
    const auto r = rhs.solutions.back().value;

    return type == bx::objective_function_type::maximize ? l > r : l < r;
}

std::string
describe(const bx::result& r)
{
    if (r.solutions.empty())
        return "no result";

    if (r.status != bx::result_status::success)
        return fmt::format("{} constraints remaining",
                           r.remaining_constraints);

    return fmt::format("{:g}", r.solutions.back().value);
}

// Coordinate search over the solver parameters: each parameter is tried in
// turn with the others fixed and keeps its best value. Half of the time
// budget goes to these trials, the other half to a final run with the best
// parameters.
bx::result
optimize(const bx::context& ctx, const bx::problem& pb)
{
    constexpr double default_budget = 10.0;

    double budget = ctx.parameters.time_limit;
    if (!std::isfinite(budget) || budget <= 0.0) {
        bx::warning(ctx,
                    "  - auto-tune (manual) needs a time limit: {}s used.\n",
                    default_budget);
        budget = default_budget;
    }

    const auto parameters = tuned_parameters();

    int trials = 1;
    for (const auto& parameter : parameters)
        trials += bx::length(parameter.values);

    bx::context internal(ctx);
    internal.start = nullptr;
    internal.update = nullptr;
    internal.finish = nullptr;
    internal.log_priority = bx::context::message_type::warning;
    internal.parameters.time_limit = 0.5 * budget / trials;

    if (ctx.start)
        ctx.start(ctx.parameters);

    auto best_parameters = internal.parameters;
    auto best = bx::itm::optimize(internal, pb);
    auto best_tuned = best;

    bx::notice(ctx,
               "  - auto-tune (manual) baseline: {}\n",
               describe(best_tuned));

    for (const auto& parameter : parameters) {
        for (const double value : parameter.values) {
            internal.parameters = best_parameters;
            parameter.set(internal.parameters, value);

            auto ret = bx::itm::optimize(internal, pb);

            bx::notice(ctx,
                       "  - auto-tune (manual) {}={:g}: {}\n",
                       parameter.name,
                       value,
                       describe(ret));

            if (is_better(ret, best_tuned, pb.type)) {
                best_tuned = ret;
                best_parameters = internal.parameters;
            }
        }
    }

    bx::notice(ctx,
               "  - auto-tune (manual) best: theta={:g} delta={:g} "
               "kappa-min={:g} kappa-step={:g} init-policy-random={:g} "
               "optimizer={}\n",
               best_parameters.theta,
               best_parameters.delta,
               best_parameters.kappa_min,
               best_parameters.kappa_step,
               best_parameters.init_policy_random,
               best_parameters.optimizer ==
                   bx::solver_parameters::optimizer_type::restart
                 ? "restart"
                 : "genetic");

    internal.parameters = best_parameters;
    internal.parameters.time_limit = 0.5 * budget;
    best = bx::itm::optimize(internal, pb);

    if (is_better(best_tuned, best, pb.type))
        best = std::move(best_tuned);

    if (ctx.finish)
        ctx.finish(best);

    return best;
}

} // namespace

namespace baryonyx {
namespace itm {

result
manual_optimize(const baryonyx::context& ctx, const baryonyx::problem& pb)
{
    baryonyx::notice(ctx, "- auto-tune parameters (manual) starts\n");
    return ::optimize(ctx, pb);
}

} // namespace itm
} // namespace baryonyx
