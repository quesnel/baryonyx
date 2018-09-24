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
#include <baryonyx/core-utils>

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
#include "result.hpp"
#include "sparse-matrix.hpp"
#include "utils.hpp"

namespace baryonyx {
namespace itm {

template<typename Float, typename Mode>
struct best_solution_recorder
{
    const context_ptr& m_ctx;
    std::chrono::time_point<std::chrono::steady_clock> m_start;
    std::mutex m_mutex;
    result m_best;

    best_solution_recorder(const context_ptr& ctx)
      : m_ctx(ctx)
      , m_start(std::chrono::steady_clock::now())
    {}

    void try_update(int remaining_constraints, int loop)
    {
        try {
            std::lock_guard<std::mutex> lock(m_mutex);

            if (store_advance(m_best, remaining_constraints)) {
                m_best.remaining_constraints = remaining_constraints;
                m_best.loop = loop;
                m_best.duration =
                  compute_duration(m_start, std::chrono::steady_clock::now());

                info(m_ctx,
                     "  - Constraints remaining: {} (i={} t={}s)\n",
                     remaining_constraints,
                     loop,
                     m_best.duration);
            }
        } catch (const std::exception& e) {
            error(m_ctx, "sync optimization error: {}", e.what());
        }
    }

    void try_update(const std::vector<bool>& solution, double value, int loop)
    {
        try {
            std::lock_guard<std::mutex> lock(m_mutex);

            if (store_solution<Mode>(m_ctx, m_best, solution, value)) {
                m_best.loop = loop;
                m_best.duration =
                  compute_duration(m_start, std::chrono::steady_clock::now());

                info(m_ctx,
                     "  - Solution found: {:G} (i={} t={}s)\n",
                     value,
                     loop,
                     m_best.duration);
            }
        } catch (const std::exception& e) {
            error(m_ctx, "sync optimization error: {}", e.what());
        }
    }
};

template<typename Solver,
         typename Float,
         typename Mode,
         typename Order,
         typename Random>
struct optimize_functor
{
    const context_ptr& m_ctx;
    Random m_rng;
    int m_thread_id;
    std::chrono::time_point<std::chrono::steady_clock> m_begin;
    std::chrono::time_point<std::chrono::steady_clock> m_end;

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
        return ((m_ctx->parameters.mode &
                 solver_parameters::mode_type::branch) ==
                solver_parameters::mode_type::branch)
                 ? run<x_counter_type>(best_recorder,
                                       constraints,
                                       variables,
                                       original_costs,
                                       norm_costs,
                                       cost_constant)
                 : run<x_type>(best_recorder,
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

        int remaining = 0;
        int best_remaining = INT_MAX;
        int i = 0;
        int pushed = -1;
        auto pushing_iteration = 0;

        auto& p = m_ctx->parameters;
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

        if (p.limit <= 0)
            p.limit = std::numeric_limits<int>::max();

        if (p.time_limit <= 0)
            p.time_limit = std::numeric_limits<double>::max();

        m_best.variables = slv.m;
        m_best.constraints = slv.n;

        m_begin = std::chrono::steady_clock::now();

        for (;;) {
            remaining = compute.run(slv, x, kappa, delta, theta);
            if (remaining == 0) {
                store_if_better(x,
                                slv.results(x, original_costs, cost_constant),
                                i,
                                best_recorder);
                pushed = 0;
            } else if (remaining < best_remaining) {
                store_if_better(x, remaining, i, best_recorder);
                best_remaining = remaining;
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
                          x,
                          slv.results(x, original_costs, cost_constant),
                          i,
                          best_recorder);
                }
            }

            if (i > p.w)
                kappa += kappa_step * std::pow(static_cast<Float>(remaining) /
                                                 static_cast<Float>(slv.m),
                                               alpha);

            ++i;
            if ((p.limit > 0 && i >= p.limit) || kappa > kappa_max ||
                pushed > p.pushes_limit) {

                if (remaining < best_remaining)
                    store_if_better(x, remaining, i, best_recorder);

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

            m_end = std::chrono::steady_clock::now();
            if (is_time_limit(p.time_limit, m_begin, m_end)) {
                if (m_best.status == result_status::uninitialized) {
                    m_best.status = result_status::time_limit_reached;
                    store_if_better(x, remaining, i, best_recorder);
                }
                break;
            }
        }

        return m_best;
    }

private:
    template<typename Xtype>
    void store_if_better(const Xtype& x,
                         int remaining,
                         int i,
                         best_solution_recorder<Float, Mode>& best_recorder)
    {
        if (store_advance(m_best, remaining)) {
            m_best.loop = i;
            m_best.remaining_constraints = remaining;
            m_best.annoying_variable = x.upper();

            best_recorder.try_update(remaining, i);
        }
    }

    template<typename Xtype>
    void store_if_better(const Xtype& x,
                         double current,
                         int i,
                         best_solution_recorder<Float, Mode>& best_recorder)
    {
        if (store_solution<Mode>(m_ctx, m_best, x.data(), current)) {
            m_best.status = result_status::success;
            m_best.loop = i;
            m_best.remaining_constraints = 0;
            m_best.annoying_variable = x.upper();

            best_recorder.try_update(x.data(), current, i);
        }
    }
};

//
// Get number of thread to use in optimizer from parameters list or
// from the standard thread API. If an error occurred, this function
// returns 1.
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
    info(ctx, "- Optimizer initializing\n");
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

        const auto thread = get_thread_number(ctx);

        std::vector<std::thread> pool(thread);
        pool.clear();
        std::vector<std::future<result>> results(thread);
        results.clear();

        best_solution_recorder<Float, Mode> result(ctx);

        if (thread == 1)
            info(ctx, "  - optimizer uses one thread\n");
        else
            info(ctx, "  - optimizer uses {} threads\n", thread);

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

        ret = std::move(result.m_best);
        std::set<solution> all_solutions;

        for (unsigned i{ 0 }; i != thread; ++i) {
            auto current = results[i].get();

            if (current.status == result_status::success)
                all_solutions.insert(current.solutions.begin(),
                                     current.solutions.end());
        }

        ret.solutions.clear();

        switch (ctx->parameters.storage) {
        case solver_parameters::storage_type::one:
            ret.solutions.push_back(*(all_solutions.rbegin()));
            break;
        case solver_parameters::storage_type::bound:
            ret.solutions.push_back(*(all_solutions.begin()));
            ret.solutions.push_back(*(all_solutions.rbegin()));
            break;
        case solver_parameters::storage_type::five: {
            int i = 0;
            for (auto& elem : all_solutions) {
                ret.solutions.push_back(*(all_solutions.rbegin()));
                ++i;

                if (i >= 5)
                    break;
            }
            break;
        }
        }
    }

    return ret;
}

} // namespace itm
} // namespace baryonyx

#endif
