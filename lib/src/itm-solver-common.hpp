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
         typename Cost,
         typename Observer>
struct solver_functor
{
    std::chrono::time_point<std::chrono::steady_clock> m_begin;
    std::chrono::time_point<std::chrono::steady_clock> m_end;

    const context_ptr& m_ctx;
    random_engine& m_rng;

    result m_best;

    solver_functor(const context_ptr& ctx,
                   random_engine& rng,
                   const problem& problem)
      : m_ctx(ctx)
      , m_rng(rng)
    {
        m_best.strings = problem.strings;
        m_best.affected_vars = problem.affected_vars;
        m_best.variable_name = problem.vars.names;
    }

    result operator()(const std::vector<merged_constraint>& constraints,
                      int variables,
                      const Cost& original_costs,
                      double cost_constant)
    {
        x_type x(variables);

        int best_remaining = INT_MAX;

        auto& p = m_ctx->parameters;

        auto norm_costs = normalize_costs<Float, Cost>(
          m_ctx, original_costs, m_rng, variables);

        const auto kappa_step = static_cast<Float>(p.kappa_step);
        const auto kappa_max = static_cast<Float>(p.kappa_max);
        const auto alpha = static_cast<Float>(p.alpha);
        const auto theta = static_cast<Float>(p.theta);
        const auto delta =
          p.delta < 0
            ? compute_delta<Float, Cost>(m_ctx, norm_costs, theta, variables)
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

        solver_initializer<Solver, Float, Mode, x_type> initializer(
          slv, x, p.init_policy, p.init_random);
        Order compute(slv, x, m_rng);

        m_best.variables = slv.m;
        m_best.constraints = slv.n;

        bool start_push = false;
        auto kappa = static_cast<Float>(p.kappa_min);

        Observer obs(slv, "img", p.limit);

        m_begin = std::chrono::steady_clock::now();
        m_end = std::chrono::steady_clock::now();

        for (int i = 0; i != p.limit; ++i) {
            auto remaining = compute.run(slv, x, kappa, delta, theta);
            obs.make_observation();

            if (remaining == 0) {
                store_if_better(
                  x, original_costs.results(x, cost_constant), i);
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
                m_best.status = result_status::kappa_max_reached;
                break;
            }

            if (is_timelimit_reached()) {
                m_best.status = result_status::time_limit_reached;
                break;
            }
        }

        if (!start_push) {
            m_best.status = result_status::limit_reached;
        } else {
            for (int push = 0; push < p.pushes_limit; ++push) {
                auto remaining =
                  compute.push_and_run(slv,
                                       x,
                                       pushing_k_factor * kappa,
                                       delta,
                                       theta,
                                       pushing_objective_amplifier);

                if (remaining == 0)
                    store_if_better(x,
                                    original_costs.results(x, cost_constant),
                                    -push * p.pushing_iteration_limit - 1);

                if (is_timelimit_reached())
                    break;

                for (int iter = 0; iter < p.pushing_iteration_limit; ++iter) {
                    remaining = compute.run(slv, x, kappa, delta, theta);

                    if (remaining == 0) {
                        store_if_better(
                          x,
                          original_costs.results(x, cost_constant),
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

            if (m_ctx->update)
                m_ctx->update(m_best);
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

            if (m_ctx->update)
                m_ctx->update(m_best);
        }
    }
};

template<typename Solver,
         typename Float,
         typename Mode,
         typename Order,
         typename Cost>
inline result
solve_problem(const context_ptr& ctx, const problem& pb)
{
    if (ctx->start)
        ctx->start(ctx->parameters);

    result ret;

    auto constraints{ make_merged_constraints(ctx, pb) };
    if (!constraints.empty() && !pb.vars.values.empty()) {
        random_engine rng(init_random_generator_seed<random_engine>(ctx));

        auto variables = numeric_cast<int>(pb.vars.values.size());
        auto cost = Cost(pb.objective, variables);
        auto cost_constant = pb.objective.value;

        switch (ctx->parameters.observer) {
        case solver_parameters::observer_type::pnm: {
            using obs = pnm_observer<Solver, Float>;

            solver_functor<Solver, Float, Mode, Order, Cost, obs> slv(
              ctx, rng, pb);

            ret = slv(constraints, variables, cost, cost_constant);
        } break;
        case solver_parameters::observer_type::file: {
            using obs = file_observer<Solver, Float>;

            solver_functor<Solver, Float, Mode, Order, Cost, obs> slv(
              ctx, rng, pb);

            ret = slv(constraints, variables, cost, cost_constant);
        } break;
        default: {
            using obs = none_observer<Solver, Float>;

            solver_functor<Solver, Float, Mode, Order, Cost, obs> slv(
              ctx, rng, pb);

            ret = slv(constraints, variables, cost, cost_constant);
            break;
        }
        }

    } else {
        ret.status = result_status::success;
    }

    ret.strings = pb.strings;
    ret.variable_name = std::move(pb.vars.names);
    ret.affected_vars = std::move(pb.affected_vars);

    if (ctx->finish)
        ctx->finish(ret);

    return ret;
}

} // namespace itm
} // namespace baryonyx

#endif
