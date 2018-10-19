/* Copyright (C) 2016-2019 INRA
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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_SOLVER_COMMON_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_SOLVER_COMMON_HPP

#include <baryonyx/core-compare>

#include <algorithm>
#include <chrono>
#include <functional>
#include <random>
#include <utility>
#include <vector>

#include "observer.hpp"
#include "private.hpp"
#include "sparse-matrix.hpp"
#include "utils.hpp"

namespace baryonyx {
namespace itm {

template<typename Solver,
         typename Float,
         typename Mode,
         typename Order,
         typename Random,
         typename Observer>
struct solver_functor
{
    std::chrono::time_point<std::chrono::steady_clock> m_begin;
    std::chrono::time_point<std::chrono::steady_clock> m_end;

    const context_ptr& m_ctx;
    Random& m_rng;

    result m_best;

    solver_functor(const context_ptr& ctx,
                   Random& rng,
                   const std::vector<std::string>& variable_names,
                   const affected_variables& affected_vars)
      : m_ctx(ctx)
      , m_rng(rng)
    {
        m_best.affected_vars = affected_vars;
        m_best.variable_name = variable_names;
    }

    result operator()(const std::vector<merged_constraint>& constraints,
                      int variables,
                      const std::unique_ptr<Float[]>& original_costs,
                      const std::unique_ptr<Float[]>& norm_costs,
                      double cost_constant)
    {
        x_type x(variables);

        int best_remaining = INT_MAX;

        auto& p = m_ctx->parameters;

        const auto kappa_step = static_cast<Float>(p.kappa_step);
        const auto kappa_max = static_cast<Float>(p.kappa_max);
        const auto alpha = static_cast<Float>(p.alpha);
        const auto theta = static_cast<Float>(p.theta);
        const auto delta =
          p.delta < 0
            ? compute_delta<Float>(m_ctx, norm_costs, theta, variables)
            : static_cast<Float>(p.delta);

        const auto pushing_k_factor = static_cast<Float>(p.pushing_k_factor);
        const auto pushing_objective_amplifier =
          static_cast<Float>(p.pushing_objective_amplifier);

        if (p.limit <= 0)
            p.limit = std::numeric_limits<int>::max();

        if (p.time_limit <= 0)
            p.time_limit = std::numeric_limits<double>::infinity();

        if (p.pushes_limit <= 0)
            p.pushes_limit = 0;

        if (p.pushing_iteration_limit <= 0)
            p.pushes_limit = 0;

        Solver slv(
          m_rng, length(constraints), variables, norm_costs, constraints);

        auto init_policy = init_solver(slv, x, p.init_policy, p.init_random);

        Order compute(slv, x, m_rng);

        m_best.variables = slv.m;
        m_best.constraints = slv.n;

        bool start_push = false;
        auto kappa = static_cast<Float>(p.kappa_min);

        init_solver(slv, x, init_policy, m_ctx->parameters.init_random);

        auto max_cost = max_cost_init(original_costs, variables, Mode());
        bounds_printer<Float, Mode> bound_print(max_cost);
        Observer obs(slv, "img", p.limit);

        m_begin = std::chrono::steady_clock::now();
        m_end = std::chrono::steady_clock::now();

        info(m_ctx, "* solver starts:\n");

        for (int i = 0; i != p.limit; ++i) {
            auto remaining = compute.run(slv, x, kappa, delta, theta);
            obs.make_observation();

            if (remaining == 0) {
                store_if_better(
                  x, results(x, original_costs, cost_constant, variables), i);
                best_remaining = remaining;
                start_push = true;
                break;
            }

            if (remaining < best_remaining) {
                store_if_better(x, remaining, i);
                best_remaining = remaining;
            }

            if (i > p.w)
                kappa += kappa_step * std::pow(static_cast<Float>(remaining) /
                                                 static_cast<Float>(slv.m),
                                               alpha);

            if (kappa > kappa_max) {
                info(m_ctx, "  - Kappa max reached: {:+.6f}\n", kappa);
                m_best.status = result_status::kappa_max_reached;
                break;
            }

            if (is_timelimit_reached()) {
                info(m_ctx, "  - Time limit reached: {} {:+.6f}\n", i, kappa);
                m_best.status = result_status::time_limit_reached;
                break;
            }
        }

        if (!start_push) {
            info(m_ctx, "  - No solution found.\n");
            m_best.status = result_status::limit_reached;
        } else {
            info(m_ctx,
                 "  - Starts pushes (limit: {} iterations: {})\n",
                 p.pushes_limit,
                 p.pushing_iteration_limit);

            for (int push = 0; push < p.pushes_limit; ++push) {
                auto remaining =
                  compute.push_and_run(slv,
                                       x,
                                       pushing_k_factor * kappa,
                                       delta,
                                       theta,
                                       pushing_objective_amplifier);

                if (remaining == 0)
                    store_if_better(
                      x,
                      results(x, original_costs, cost_constant, variables),
                      -push * p.pushing_iteration_limit - 1);

                if (is_timelimit_reached())
                    break;

                for (int iter = 0; iter < p.pushing_iteration_limit; ++iter) {
                    remaining = compute.run(slv, x, kappa, delta, theta);

                    if (remaining == 0) {
                        store_if_better(
                          x,
                          results(x, original_costs, cost_constant, variables),
                          -push * p.pushing_iteration_limit - iter - 1);
                        break;
                    }

                    if (iter > p.w)
                        kappa +=
                          kappa_step * std::pow(static_cast<Float>(remaining) /
                                                  static_cast<Float>(slv.m),
                                                alpha);

                    if (kappa > kappa_max)
                        break;

                    if (is_timelimit_reached())
                        break;
                }
            }
        }

        return m_best;
    }

private:
    bool is_timelimit_reached()
    {
        m_end = std::chrono::steady_clock::now();

        return is_time_limit(m_ctx->parameters.time_limit, m_begin, m_end);
    }

    double duration()
    {
        m_end = std::chrono::steady_clock::now();

        return compute_duration(m_begin, m_end);
    }

    template<typename Xtype>
    void store_if_better(const Xtype& x, int remaining, int i)
    {
        if (store_advance(m_best, remaining)) {
            m_best.loop = i;
            m_best.remaining_constraints = remaining;
            m_best.annoying_variable = x.upper();
            m_best.duration = duration();

            info(m_ctx,
                 "  - Constraints remaining: {} (i={} t={}s)\n",
                 remaining,
                 i,
                 m_best.duration);
        }
    }

    template<typename Xtype>
    void store_if_better(const Xtype& x, double current, int i)
    {
        if (store_solution<Mode>(m_ctx, m_best, x.data(), current)) {
            m_best.status = result_status::success;
            m_best.loop = i;
            m_best.remaining_constraints = 0;
            m_best.duration = duration();

            if (i >= 0)
                info(m_ctx,
                     "  - Solution found: {:f} (i={} t={}s)\n",
                     current,
                     i,
                     m_best.duration);
            else
                info(m_ctx,
                     "  - Solution found via push: {:f} (i={} t={}s)\n",
                     current,
                     i,
                     m_best.duration);
        }
    }
};

template<typename Solver,
         typename Float,
         typename Mode,
         typename Order,
         typename Random>
inline result
solve_problem(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "- Solver initializing\n");
    print(ctx);

    result ret;
    auto affected_vars = std::move(pb.affected_vars);

    auto constraints{ make_merged_constraints(ctx, pb) };
    if (!constraints.empty() && !pb.vars.values.empty()) {
        Random rng(init_random_generator_seed<Random>(ctx));

        auto variables = numeric_cast<int>(pb.vars.values.size());
        auto cost = make_objective_function<Float>(pb.objective, variables);
        auto norm_costs =
          normalize_costs<Float, Random>(ctx, cost, rng, variables);
        auto cost_constant = pb.objective.value;
        auto names = std::move(pb.vars.names);

        switch (ctx->parameters.observer) {
        case solver_parameters::observer_type::pnm: {
            using obs = pnm_observer<Solver, Float>;

            solver_functor<Solver, Float, Mode, Order, Random, obs> slv(
              ctx, rng, names, affected_vars);

            ret = slv(constraints, variables, cost, norm_costs, cost_constant);
        } break;
        case solver_parameters::observer_type::file: {
            using obs = file_observer<Solver, Float>;

            solver_functor<Solver, Float, Mode, Order, Random, obs> slv(
              ctx, rng, names, affected_vars);

            ret = slv(constraints, variables, cost, norm_costs, cost_constant);
        } break;
        default: {
            using obs = none_observer<Solver, Float>;

            solver_functor<Solver, Float, Mode, Order, Random, obs> slv(
              ctx, rng, names, affected_vars);

            ret = slv(constraints, variables, cost, norm_costs, cost_constant);
            break;
        }
        }

        ret.variable_name = std::move(names);
    } else {
        ret.status = result_status::success;
    }
    ret.affected_vars = std::move(affected_vars);

    return ret;
}

} // namespace itm
} // namespace baryonyx

#endif
