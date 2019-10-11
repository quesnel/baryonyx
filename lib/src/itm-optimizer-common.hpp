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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_OPTIMIZER_COMMON_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_OPTIMIZER_COMMON_HPP

#include <baryonyx/core-compare>
#include <baryonyx/core-utils>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <functional>
#include <future>
#include <set>
#include <thread>
#include <utility>
#include <vector>

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

                if (m_ctx->update)
                    m_ctx->update(m_best);
            }
        } catch (const std::exception& e) {
            error(m_ctx, "sync optimization error: {}", e.what());
        }
    }

    void try_update(const bit_array& solution, double value, int loop)
    {
        try {
            std::lock_guard<std::mutex> lock(m_mutex);

            if (store_solution<Mode>(m_ctx, m_best, solution, value)) {
                m_best.loop = loop;
                m_best.duration =
                  compute_duration(m_start, std::chrono::steady_clock::now());

                if (m_ctx->update)
                    m_ctx->update(m_best);
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
         typename Cost>
struct optimize_functor
{
    const context_ptr& m_ctx;
    random_engine m_rng;
    int m_thread_id;

    const std::vector<std::string_view>& variable_names;
    const affected_variables& affected_vars;

    result m_best;

    optimize_functor(const context_ptr& ctx,
                     unsigned thread_id,
                     typename random_engine::result_type seed,
                     const std::vector<std::string_view>& variable_names_,
                     const affected_variables& affected_vars_)
      : m_ctx(ctx)
      , m_rng(seed)
      , m_thread_id(thread_id)
      , variable_names(variable_names_)
      , affected_vars(affected_vars_)
    {}

    result operator()(std::atomic_bool& stop_task,
                      best_solution_recorder<Float, Mode>& best_recorder,
                      const std::vector<merged_constraint>& constraints,
                      int variables,
                      const Cost& original_costs,
                      double cost_constant)
    {
        return ((m_ctx->parameters.mode &
                 solver_parameters::mode_type::branch) ==
                solver_parameters::mode_type::branch)
                 ? run<x_counter_type>(stop_task,
                                       best_recorder,
                                       constraints,
                                       variables,
                                       original_costs,
                                       cost_constant)
                 : run<x_type>(stop_task,
                               best_recorder,
                               constraints,
                               variables,
                               original_costs,
                               cost_constant);
    }

    template<typename Xtype>
    result run(std::atomic_bool& stop_task,
               best_solution_recorder<Float, Mode>& best_recorder,
               const std::vector<merged_constraint>& constraints,
               int variables,
               const Cost& original_costs,
               double cost_constant)
    {
        Xtype x(variables);

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

        if (p.pushes_limit <= 0)
            p.pushes_limit = 0;

        if (p.pushing_iteration_limit <= 0)
            p.pushes_limit = 0;

        Solver slv(
          m_rng, length(constraints), variables, norm_costs, constraints);

        solver_initializer<Solver, Float, Mode, Xtype> initializer(
          slv, x, p.init_policy, p.init_random);

        Order compute(slv, x, m_rng);

        m_best.variables = slv.m;
        m_best.constraints = slv.n;

        auto x_is_solution{ false };

        while (!stop_task.load()) {
            auto kappa = static_cast<Float>(p.kappa_min);

            initializer.reinit(slv, x, x_is_solution, m_best);

            for (int i = 0; !stop_task.load() && i != p.limit; ++i) {
                auto remaining = compute.run(slv, x, kappa, delta, theta);

                if (remaining == 0) {
                    x_is_solution = true;
                    store_if_better(x,
                                    original_costs.results(x, cost_constant),
                                    i,
                                    best_recorder);
                    best_remaining = remaining;
                    break;
                }

                x_is_solution = false;

                if (remaining < best_remaining) {
                    store_if_better(x, remaining, i, best_recorder);
                    best_remaining = remaining;
                }

                if (i > p.w)
                    kappa +=
                      kappa_step * std::pow(static_cast<Float>(remaining) /
                                              static_cast<Float>(slv.m),
                                            alpha);

                if (kappa > kappa_max)
                    break;
            }

            if (best_remaining > 0)
                continue;

            for (int push = 0; !stop_task.load() && push < p.pushes_limit;
                 ++push) {
                x_is_solution = false;

                auto remaining =
                  compute.push_and_run(slv,
                                       x,
                                       pushing_k_factor * kappa,
                                       delta,
                                       theta,
                                       pushing_objective_amplifier);

                if (remaining == 0) {
                    x_is_solution = true;
                    store_if_better(x,
                                    original_costs.results(x, cost_constant),
                                    -push * p.pushing_iteration_limit - 1,
                                    best_recorder);
                }

                for (int iter = 0;
                     !stop_task.load() && iter < p.pushing_iteration_limit;
                     ++iter) {
                    remaining = compute.run(slv, x, kappa, delta, theta);

                    if (remaining == 0) {
                        x_is_solution = true;
                        store_if_better(
                          x,
                          original_costs.results(x, cost_constant),
                          -push * p.pushing_iteration_limit - iter - 1,
                          best_recorder);
                        break;
                    }

                    if (iter > p.w)
                        kappa +=
                          kappa_step * std::pow(static_cast<Float>(remaining) /
                                                  static_cast<Float>(slv.m),
                                                alpha);

                    if (kappa > kappa_max)
                        break;
                }
            }
        }

        if (m_best.status == result_status::uninitialized && stop_task.load())
            m_best.status = result_status::time_limit_reached;

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
         typename Cost>
inline result
optimize_problem(const context_ptr& ctx, const problem& pb)
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

        const auto thread = get_thread_number(ctx);

        std::vector<std::thread> pool(thread);
        pool.clear();
        std::vector<std::future<result>> results(thread);
        results.clear();

        best_solution_recorder<Float, Mode> result(ctx);

        auto seeds = generate_seed(rng, thread);

        std::atomic_bool stop_task;

        stop_task.store(false);

        for (unsigned i{ 0 }; i != thread; ++i) {
            std::packaged_task<baryonyx::result()> task(
              std::bind(optimize_functor<Solver, Float, Mode, Order, Cost>(
                          ctx, i, seeds[i], pb.vars.names, pb.affected_vars),
                        std::ref(stop_task),
                        std::ref(result),
                        std::ref(constraints),
                        variables,
                        std::ref(cost),
                        cost_constant));

            results.emplace_back(task.get_future());

            pool.emplace_back(std::thread(std::move(task)));
        }

        // User limits the optimization process. The current thread turns into
        // sleeping mode for the specified duration time in second.  The
        // following code takes into account the first three decimals only.

        if (ctx->parameters.time_limit > 0) {
            const auto duration =
              static_cast<long>(ctx->parameters.time_limit * 1000.0);
            std::this_thread::sleep_for(std::chrono::milliseconds{ duration });
            stop_task.store(true);
        }

        for (auto& t : pool)
            t.join();

        ret = std::move(result.m_best);
        ret.affected_vars = pb.affected_vars;
        ret.variable_name = pb.vars.names;

        std::set<solution> all_solutions;

        for (unsigned i{ 0 }; i != thread; ++i) {
            auto current = results[i].get();

            if (current.status == result_status::success)
                all_solutions.insert(current.solutions.begin(),
                                     current.solutions.end());
        }

        ret.solutions.clear();

        if (!all_solutions.empty()) {
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
                    ret.solutions.push_back(elem);
                    ++i;

                    if (i >= 5)
                        break;
                }
                break;
            }
            }
        }
    }

    if (ctx->finish)
        ctx->finish(ret);

    return ret;
}

} // namespace itm
} // namespace baryonyx

#endif
