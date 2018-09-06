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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_OPTIMIZER_COMMON_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_OPTIMIZER_COMMON_HPP

#include <baryonyx/core-compare>

#include <algorithm>
#include <chrono>
#include <functional>
#include <future>
#include <random>
#include <set>
#include <thread>
#include <utility>
#include <vector>

#include "observer.hpp"
#include "private.hpp"
#include "sparse-matrix.hpp"
#include "utils.hpp"

namespace baryonyx {
namespace itm {

template<typename Float, typename Mode>
struct best_solution_recorder
{
    const context_ptr& m_ctx;
    std::mutex m_mutex;
    result m_best;

    best_solution_recorder(const context_ptr& ctx)
      : m_ctx(ctx)
    {}

    bool try_update(const result& current) noexcept
    {
        try {
            std::lock_guard<std::mutex> lock(m_mutex);

            if (current.status != result_status::success)
                return false;

            if (m_best.status != result_status::success ||
                is_better_solution(current, m_best, Mode())) {

                info(m_ctx,
                     "  - Solution found: {} (i={} t={}s)\n",
                     best_solution_value(current),
                     current.loop,
                     current.duration);

                m_best = current;

                return true;
            }
        } catch (const std::exception& e) {
            error(m_ctx, "sync optimization error: {}", e.what());
        }

        return false;
    }
};

template<typename Solver,
         typename Float,
         typename Mode,
         typename Order,
         typename Random>
struct optimize_functor
{
    std::chrono::time_point<std::chrono::steady_clock> m_begin;
    std::chrono::time_point<std::chrono::steady_clock> m_end;

    std::set<solution> m_all_solutions;

    const context_ptr& m_ctx;
    Random m_rng;
    int m_thread_id;

    result m_best;

    optimize_functor(const context_ptr& ctx,
                     unsigned thread_id,
                     typename Random::result_type seed,
                     const std::vector<std::string>& variable_names,
                     const affected_variables& affected_vars)
      : m_ctx(ctx)
      , m_rng(seed)
      , m_thread_id(thread_id)
    {
        m_best.affected_vars = affected_vars;
        m_best.variable_name = variable_names;
    }

    result operator()(best_solution_recorder<Float, Mode>& best_recorder,
                      const std::vector<merged_constraint>& constraints,
                      int variables,
                      const std::unique_ptr<Float[]>& original_costs,
                      const std::unique_ptr<Float[]>& norm_costs,
                      double cost_constant)
    {
        return (m_ctx->parameters.mode == solver_parameters::mode_type::branch)
                 ? run<x_type>(best_recorder,
                               constraints,
                               variables,
                               original_costs,
                               norm_costs,
                               cost_constant)
                 : run<x_counter_type>(best_recorder,
                                       constraints,
                                       variables,
                                       original_costs,
                                       norm_costs,
                                       cost_constant);
    }

    template<typename Xtype>
    result run(best_solution_recorder<Float, Mode>& best_recorder,
               const std::vector<merged_constraint>& constraints,
               int variables,
               const std::unique_ptr<Float[]>& original_costs,
               const std::unique_ptr<Float[]>& norm_costs,
               double cost_constant)
    {
        Xtype x(variables);

        m_begin = std::chrono::steady_clock::now();
        m_end = m_begin;

        int i = 0;
        int pushed = -1;
        int pushing_iteration = 0;
        const auto& p = m_ctx->parameters;

        auto init_policy = p.init_policy;
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

        init_policy = init_solver(slv, x, init_policy, p.init_random);

        Order compute(slv, x, m_rng);

        for (; !is_time_limit(p.time_limit, m_begin, m_end);
             m_end = std::chrono::steady_clock::now(), ++i) {

            int remaining = compute.run(slv, x, kappa, delta, theta);
            if (remaining == 0)
                store_if_better(slv.results(x, original_costs, cost_constant),
                                x,
                                i,
                                best_recorder);

            if (i > p.w)
                kappa += kappa_step * std::pow(static_cast<Float>(remaining) /
                                                 static_cast<Float>(slv.m),
                                               alpha);

            if ((p.limit > 0 && i >= p.limit) || kappa > kappa_max ||
                pushed > p.pushes_limit) {
                if (m_best.solutions.empty())
                    init_policy =
                      init_solver(slv, x, init_policy, p.init_random);
                else
                    init_policy =
                      init_solver(slv,
                                  x,
                                  m_best.solutions.back().variables,
                                  init_policy,
                                  p.init_random);

                i = 0;
                kappa = static_cast<Float>(kappa_min);
                pushed = -1;
                pushing_iteration = 0;

                continue;
            }

            if (pushed >= 0) {
                ++pushing_iteration;

                if (pushing_iteration >= p.pushing_iteration_limit) {
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
                          slv.results(x, original_costs, cost_constant),
                          x,
                          i,
                          best_recorder);
                }
            }
        }

        std::copy(m_all_solutions.begin(),
                  m_all_solutions.end(),
                  std::back_inserter(m_best.solutions));

        return m_best;
    }

private:
    template<typename Xtype>
    void store_if_better(double current,
                         const Xtype& x,
                         int i,
                         best_solution_recorder<Float, Mode>& best_recorder)
    {
        if (m_best.solutions.empty() ||
            is_better_solution(
              current, m_best.solutions.back().value, Mode())) {
            m_best.status = baryonyx::result_status::success;

            // Store only the best solution, other solutions are stored into
            // the @c std::set to avoid duplicated solutions.

            if (m_best.solutions.empty())
                m_best.solutions.emplace_back(x.data(), current);
            else
                m_best.solutions[0] = { x.data(), current };

            m_best.duration = compute_duration(m_begin, m_end);
            m_best.loop = i;
            m_best.annoying_variable = x.upper();

            best_recorder.try_update(m_best);
        }

        m_all_solutions.emplace(x.data(), current);
    }
};

//
// Get number of thread to use in optimizer from parameters list or from the
// standard thread API. If an error occurred, this function returns 1.
//
inline unsigned
get_thread_number(const baryonyx::context_ptr& ctx) noexcept
{
    unsigned ret;

    if (ctx->parameters.thread <= 0)
        ret = std::thread::hardware_concurrency();
    else
        ret = static_cast<unsigned>(ctx->parameters.thread);

    if (ret == 0)
        return 1;

    return ret;
}

template<typename Solver,
         typename Float,
         typename Mode,
         typename Order,
         typename Random>
inline result
optimize_problem(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "Optimizer initializing\n");

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

        const auto thread = get_thread_number(ctx);

        std::vector<std::thread> pool(thread);
        pool.clear();
        std::vector<std::future<result>> results(thread);
        results.clear();

        best_solution_recorder<Float, Mode> result(ctx);

        if (thread == 1)
            info(ctx, "optimizer starts with one thread\n");
        else
            info(ctx, "Optimizer starts with {} threads\n", thread);

        auto seeds = generate_seed(rng, thread);

        for (unsigned i{ 0 }; i != thread; ++i) {
            std::packaged_task<baryonyx::result()> task(
              std::bind(optimize_functor<Solver, Float, Mode, Order, Random>(
                          ctx, i, seeds[i], names, affected_vars),
                        std::ref(result),
                        std::ref(constraints),
                        variables,
                        std::ref(cost),
                        std::ref(norm_costs),
                        cost_constant));

            results.emplace_back(task.get_future());

            pool.emplace_back(std::thread(std::move(task)));
        }

        for (auto& t : pool)
            t.join();

        std::set<solution> all_solutions;

        for (unsigned i{ 0 }; i != thread; ++i) {
            auto current = results[i].get();
            if (current.status == result_status::success) {
                all_solutions.insert(current.solutions.begin(),
                                     current.solutions.end());

                if (ret.solutions.empty() ||
                    is_better_solution(current.solutions.back().value,
                                       ret.solutions.back().value,
                                       Mode()))
                    ret = current;
            }
        }

        ret.solutions.clear();

        std::copy(all_solutions.begin(),
                  all_solutions.end(),
                  std::back_inserter(ret.solutions));
    }

    return ret;
}

} // namespace itm
} // namespace baryonyx

#endif
