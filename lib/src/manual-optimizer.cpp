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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>

#include "debug.hpp"
#include "itm-common.hpp"
#include "itm.hpp"
#include "private.hpp"
#include "problem.hpp"
#include "random-distributions.hpp"
#include "random-engine.hpp"
#include "utils.hpp"

namespace {

namespace bx = baryonyx;

struct candidate
{
    bx::solver_parameters parameters;
    bx::result result;
};

template<typename T, std::size_t N>
T
pick(bx::xoshiro256pp& rng, const T (&values)[N])
{
    return values[bx::uniform_below(rng, N)];
}

// Random parameters in the ranges the solver is usually tuned in.
bx::solver_parameters
random_parameters(bx::solver_parameters p, bx::xoshiro256pp& rng)
{
    using optimizer_type = bx::solver_parameters::optimizer_type;

    constexpr optimizer_type optimizers[] = { optimizer_type::genetic,
                                              optimizer_type::restart };
    constexpr double thetas[] = { 0.2, 0.35, 0.5, 0.65, 0.8 };
    constexpr double deltas[] = { -1.0, 1e-3, 1e-2, 1e-1 };
    constexpr double kappa_mins[] = { 0.0, 0.05, 0.1, 0.2 };
    constexpr double kappa_steps[] = { 1e-4, 3e-4, 1e-3, 3e-3, 1e-2 };
    constexpr double init_randoms[] = { 0.0, 0.25, 0.5, 0.75 };

    p.optimizer = pick(rng, optimizers);
    p.theta = pick(rng, thetas);
    p.delta = pick(rng, deltas);
    p.kappa_min = pick(rng, kappa_mins);
    p.kappa_step = pick(rng, kappa_steps);
    p.init_policy_random = pick(rng, init_randoms);

    return p;
}

std::string
describe(const bx::solver_parameters& p)
{
    return fmt::format(
      "optimizer={} theta={:g} delta={:g} kappa-min={:g} kappa-step={:g} "
      "init-policy-random={:g}",
      p.optimizer == bx::solver_parameters::optimizer_type::restart
        ? "restart"
        : "genetic",
      p.theta,
      p.delta,
      p.kappa_min,
      p.kappa_step,
      p.init_policy_random);
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

// Successive halving: every remaining candidate runs for the same time, the
// better half survives and the time per run doubles, so each round costs the
// same and the last survivors are compared on long, less noisy runs. Half of
// the budget goes to these rounds, the rest to a final run of the winner.
bx::result
optimize(const bx::context& ctx, const bx::problem& pb)
{
    constexpr double default_budget = 10.0;
    constexpr double minimal_run = 0.1;

    double budget = ctx.parameters.time_limit;
    if (!std::isfinite(budget) || budget <= 0.0) {
        bx::warning(ctx,
                    "  - auto-tune (manual) needs a time limit: {}s used.\n",
                    default_budget);
        budget = default_budget;
    }

    const double tuning = 0.5 * budget;

    // The largest power of two (4 to 32) whose first round still gives each
    // candidate a meaningful run.
    int candidates = 4;
    while (candidates < 32 &&
           tuning / (2 * candidates * std::log2(2.0 * candidates)) >=
             minimal_run)
        candidates *= 2;

    const int rounds = static_cast<int>(std::log2(candidates));
    double run_time = tuning / (candidates * rounds);

    bx::xoshiro256pp rng(
      ctx.parameters.seed > 0
        ? static_cast<bx::xoshiro256pp::result_type>(ctx.parameters.seed)
        : 1u);

    std::vector<candidate> pool;
    pool.reserve(static_cast<std::size_t>(candidates));
    pool.push_back({ ctx.parameters, {} });
    while (bx::length(pool) != candidates)
        pool.push_back({ random_parameters(ctx.parameters, rng), {} });

    bx::context internal(ctx);
    internal.start = nullptr;
    internal.update = nullptr;
    internal.finish = nullptr;
    internal.log_priority = bx::context::message_type::warning;

    if (ctx.start)
        ctx.start(ctx.parameters);

    const auto start = std::chrono::steady_clock::now();
    bx::result best;

    for (int round = 0; round != rounds; ++round) {
        for (auto& c : pool) {
            internal.parameters = c.parameters;
            internal.parameters.time_limit = run_time;
            c.result = bx::itm::optimize(internal, pb);

            if (is_better(c.result, best, pb.type))
                best = c.result;
        }

        std::stable_sort(
          pool.begin(), pool.end(), [&pb](const auto& lhs, const auto& rhs) {
              return is_better(lhs.result, rhs.result, pb.type);
          });

        bx::notice(ctx,
                   "  - auto-tune (manual) round {}: {} candidates, {:.3g}s "
                   "each, best {} with {}\n",
                   round + 1,
                   pool.size(),
                   run_time,
                   describe(pool.front().result),
                   describe(pool.front().parameters));

        pool.resize(std::max<std::size_t>(1, pool.size() / 2));
        run_time *= 2;
    }

    const double elapsed = std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - start)
                             .count();

    internal.parameters = pool.front().parameters;
    internal.parameters.time_limit = std::max(budget - elapsed, minimal_run);
    auto final_result = bx::itm::optimize(internal, pb);

    if (is_better(final_result, best, pb.type))
        best = std::move(final_result);

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
