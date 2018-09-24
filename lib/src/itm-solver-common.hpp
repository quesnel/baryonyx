/* Copyright (C) 2016-2018 INRA
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

        m_begin = std::chrono::steady_clock::now();
        m_end = m_begin;

        int i = 0;
        int pushed = -1;
        int best_remaining = -1;
        int pushing_iteration = 0;
        const auto& p = m_ctx->parameters;

        const auto kappa_min = static_cast<Float>(p.kappa_min);
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

        auto kappa = kappa_min;

        Solver slv(
          m_rng, length(constraints), variables, norm_costs, constraints);

        init_solver(slv, x, p.init_policy, p.init_random);

        Order compute(slv, x, m_rng);

        auto max_cost = max_cost_init(original_costs, variables, Mode());
        bounds_printer<Float, Mode> bound_print(max_cost);
        Observer obs(slv, "img", p.limit);

        info(m_ctx, "* solver starts:\n");

        m_best.variables = slv.m;
        m_best.constraints = slv.n;

        bool start_pushing = false;

        for (;;) {
            int remaining = compute.run(slv, x, kappa, delta, theta);
            obs.make_observation();

            if (best_remaining == -1 || remaining < best_remaining) {
                best_remaining = remaining;
                m_best.duration = compute_duration(m_begin, m_end);
                m_best.loop = i;
                m_best.remaining_constraints = remaining;

                if (remaining == 0) {
                    m_best.status = baryonyx::result_status::success;
                    store_if_better(
                      slv.results(x, original_costs, cost_constant), x, i);
                    start_pushing = true;
                } else {

                    // bound_print(slv, m_ctx, m_best);

                    info(m_ctx,
                         "  - violated constraints: {}/{} at {}s (loop: {})\n",
                         remaining,
                         slv.m,
                         compute_duration(m_begin, m_end),
                         i);
                }
            }

#ifndef BARYONYX_FULL_OPTIMIZATION
            print_solver(slv, x, m_ctx, m_best.variable_name, p.print_level);
#endif

            if (start_pushing) {
                if (pushed == -1)
                    info(m_ctx,
                         "  - start push system (push: {} iteration: {}\n",
                         p.pushes_limit,
                         p.pushing_iteration_limit);

                if (pushing_iteration == 0 ||
                    pushing_iteration >= p.pushing_iteration_limit) {
                    pushed++;
                    pushing_iteration = 0;

                    remaining =
                      compute.push_and_run(slv,
                                           x,
                                           pushing_k_factor * kappa,
                                           delta,
                                           theta,
                                           pushing_objective_amplifier);

                    if (remaining == 0)
                        store_if_better(
                          slv.results(x, original_costs, cost_constant), x, i);
                } else {
                    ++pushing_iteration;
                }

                if (pushed > p.pushes_limit) {
                    info(
                      m_ctx,
                      "    - Push system limit reached. Solution found: {}\n",
                      m_best.solutions.back().value);
                    return m_best;
                }
            }

            if (i > p.w)
                kappa += kappa_step * std::pow(static_cast<Float>(remaining) /
                                                 static_cast<Float>(slv.m),
                                               alpha);

            ++i;
            if (p.limit > 0 && i > p.limit) {
                info(m_ctx, "  - Loop limit reached: {}\n", i);
                if (pushed == -1)
                    m_best.status = result_status::limit_reached;

                // if (context_get_integer_parameter(m_ctx, "print-level", 0) >
                // 0)
                //     print_missing_constraint(m_ctx,
                //                              slv.ap,
                //                              slv.A,
                //                              m_best.variable_value,
                //                              slv.b,
                //                              m_variable_names);

                return m_best;
            }

            if (kappa > kappa_max) {
                info(m_ctx, "  - Kappa max reached: {:+.6f}\n", kappa);
                if (pushed == -1)
                    m_best.status = result_status::kappa_max_reached;

                // if (context_get_integer_parameter(m_ctx, "print-level", 0) >
                // 0)
                //     print_missing_constraint(m_ctx,
                //                              slv.ap,
                //                              slv.A,
                //                              m_best.variable_value,
                //                              slv.b,
                //                              m_variable_names);

                return m_best;
            }

            m_end = std::chrono::steady_clock::now();
            if (is_time_limit(p.time_limit, m_begin, m_end)) {
                info(m_ctx, "  - Time limit reached: {} {:+.6f}\n", i, kappa);
                if (pushed == -1)
                    m_best.status = result_status::time_limit_reached;

                // if (context_get_integer_parameter(m_ctx, "print-level", 0) >
                // 0)
                //     print_missing_constraint(m_ctx,
                //                              slv.ap,
                //                              slv.A,
                //                              m_best.variable_value,
                //                              slv.b,
                //                              m_variable_names);

                return m_best;
            }
        }
    }

private:
    void store_if_better(double current, const x_type& x, int i)
    {
        if (store_solution<Mode>(m_ctx, m_best, x.data(), current)) {
            m_best.solutions.emplace_back(x.data(), current);
            m_best.duration = compute_duration(m_begin, m_end);
            m_best.loop = i;

            info(m_ctx,
                 "  - Best solution found: {:+.6f} (i={} t={}s)\n",
                 current,
                 m_best.loop,
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
