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
    mutable std::mutex m_mutex;

    std::multiset<raw_result<Mode>, raw_result_compare<Mode>> m_solutions;
    std::atomic_bool m_have_solution{ false };
    std::atomic_int m_remaining_constraints{ std::numeric_limits<int>::max() };
    double m_value{ bad_value<Mode, double>() };

    best_solution_recorder(const context_ptr& ctx)
      : m_ctx(ctx)
    {
        m_start = std::chrono::steady_clock::now();
    }

    bool copy_best_solution(bit_array& x) const noexcept
    {
        try {
            std::lock_guard<std::mutex> lock(m_mutex);

            if (m_solutions.empty() || m_solutions.begin()->x.empty())
                return false;

            x = m_solutions.begin()->x;
        } catch (const std::exception& e) {
            error(m_ctx, "copy_best_solution: {}", e.what());
            return false;
        }

        return true;
    }

    bool copy_best_solution(bit_array& x, int from, int to) const noexcept
    {
        try {
            std::lock_guard<std::mutex> lock(m_mutex);

            if (m_solutions.empty() || m_solutions.begin()->x.empty())
                return false;

            x.assign(m_solutions.begin()->x, from, to);
        } catch (const std::exception& e) {
            error(m_ctx, "copy_best_solution: {}", e.what());
            return false;
        }

        return true;
    }

    template<typename Random>
    bool copy_any_solution(bit_array& x, Random& rng) const noexcept
    {
        try {
            std::lock_guard<std::mutex> lock(m_mutex);

            if (m_solutions.empty() || m_solutions.begin()->x.empty())
                return false;

            const std::size_t zero{ 0 };
            const std::size_t max{ std::size(m_solutions) - 1 };

            std::uniform_int_distribution<std::size_t> dist(zero, max);

            auto it = m_solutions.begin();
            std::advance(it, dist(rng));
            x = it->x;
        } catch (const std::exception& e) {
            error(m_ctx, "copy_any_solution: {}", e.what());
            return false;
        }

        return true;
    }

    template<typename Random>
    bool copy_any_solution(bit_array& x, Random& rng, int from, int to) const
      noexcept
    {
        try {
            std::lock_guard<std::mutex> lock(m_mutex);

            if (m_solutions.empty() || m_solutions.begin()->x.empty())
                return false;

            const std::size_t zero{ 0 };
            const std::size_t max{ std::size(m_solutions) - 1 };

            std::uniform_int_distribution<std::size_t> dist(zero, max);

            auto it = m_solutions.begin();
            std::advance(it, dist(rng));
            x.assign(it->x, from, to);
        } catch (const std::exception& e) {
            error(m_ctx, "copy_any_solution: {}", e.what());
            return false;
        }

        return true;
    }

    void try_advance(const bit_array& solution,
                     int remaining_constraints,
                     int loop)
    {
        if (m_have_solution)
            return;

        try {
            std::lock_guard<std::mutex> lock(m_mutex);

            if (m_remaining_constraints > remaining_constraints) {
                m_remaining_constraints = remaining_constraints;
                const auto end = std::chrono::steady_clock::now();
                const auto duration = compute_duration(m_start, end);

                m_solutions.emplace(
                  solution, remaining_constraints, duration, loop);

                if (m_ctx->update)
                    m_ctx->update(remaining_constraints, 0.0, loop, duration);
            }
        } catch (const std::exception& e) {
            error(m_ctx, "sync optimization error: {}", e.what());
        }
    }

    void try_update(const bit_array& solution, double value, int loop)
    {
        try {
            std::lock_guard<std::mutex> lock(m_mutex);
            if (!m_have_solution) {
                m_solutions.clear();
                m_have_solution = true;
            }

            auto it = m_solutions.find(value);
            if (it == m_solutions.end() || it->x != solution) {
                const auto end = std::chrono::steady_clock::now();
                const auto duration = compute_duration(m_start, end);

                auto ret =
                  m_solutions.emplace(solution, value, duration, loop);

                if (m_ctx->update && ret == m_solutions.begin())
                    m_ctx->update(0, value, loop, duration);

                info(m_ctx, "  - {} solutions found\n", m_solutions.size());
            }
        } catch (const std::exception& e) {
            error(m_ctx, "sync optimization error: {}", e.what());
        }
    }
};

template<typename Solver, typename Float, typename Mode, typename Cost>
struct optimize_functor
{
    const context_ptr& m_ctx;
    random_engine m_rng;
    int m_thread_id;
    raw_result<Mode> m_best;

    optimize_functor(const context_ptr& ctx,
                     unsigned thread_id,
                     typename random_engine::result_type seed)
      : m_ctx(ctx)
      , m_rng(seed)
      , m_thread_id(thread_id)
    {}

    void operator()(const std::atomic_bool& stop_task,
                    best_solution_recorder<Float, Mode>& best_recorder,
                    const std::vector<merged_constraint>& constraints,
                    int variables,
                    const Cost& original_costs,
                    double cost_constant)
    {
        bit_array x(variables);

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

        const int w_limit = static_cast<int>(p.w);

        if (p.limit <= 0)
            p.limit = std::numeric_limits<int>::max();

        if (p.pushes_limit <= 0)
            p.pushes_limit = 0;

        if (p.pushing_iteration_limit <= 0)
            p.pushes_limit = 0;

        Solver slv(
          m_rng, length(constraints), variables, norm_costs, constraints);

        solver_initializer<Solver, Float, Mode> initializer(
          slv, p.init_policy, p.init_policy_random, p.init_random);

        compute_order compute(p.order, variables);
        compute.init(slv, x);
        bool is_a_solution = false;

        while (!stop_task.load()) {
            // auto kappa = static_cast<Float>(p.kappa_min);

            auto kappa = initializer.reinit(
              slv, x, is_a_solution, best_recorder, p.kappa_min, p.kappa_max);
            auto best_remaining = INT_MAX;

            is_a_solution = false;

            for (int i = 0; !stop_task.load() && i != p.limit; ++i) {
                auto remaining =
                  compute.run(slv, x, m_rng, kappa, delta, theta);

                if (remaining == 0) {
                    store_if_better(x,
                                    original_costs.results(x, cost_constant),
                                    i,
                                    best_recorder);
                    best_remaining = remaining;
                    is_a_solution = true;
                    break;
                }

                if (remaining < best_remaining) {
                    advance_if_better(x, remaining, i, best_recorder);
                    best_remaining = remaining;
                }

                if (i > w_limit)
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

                auto remaining =
                  compute.push_and_run(slv,
                                       x,
                                       m_rng,
                                       pushing_k_factor * kappa,
                                       delta,
                                       theta,
                                       pushing_objective_amplifier);

                if (remaining == 0) {
                    store_if_better(x,
                                    original_costs.results(x, cost_constant),
                                    -push * p.pushing_iteration_limit - 1,
                                    best_recorder);
                    is_a_solution = true;
                }

                for (int iter = 0;
                     !stop_task.load() && iter < p.pushing_iteration_limit;
                     ++iter) {
                    remaining =
                      compute.run(slv, x, m_rng, kappa, delta, theta);

                    if (remaining == 0) {
                        store_if_better(
                          x,
                          original_costs.results(x, cost_constant),
                          -push * p.pushing_iteration_limit - iter - 1,
                          best_recorder);
                        is_a_solution = true;
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
    }

private:
    void advance_if_better(const bit_array& x,
                           int remaining,
                           int i,
                           best_solution_recorder<Float, Mode>& best_recorder)
    {
        if (m_best.remaining_constraints > remaining) {
            m_best.x = x;
            m_best.loop = i;
            m_best.remaining_constraints = remaining;

            best_recorder.try_advance(x, remaining, i);
        }
    }

    void store_if_better(const bit_array& x,
                         double current,
                         int i,
                         best_solution_recorder<Float, Mode>& best_recorder)
    {
        if (is_better_solution<Mode>(current, m_best.value)) {
            m_best.x = x;
            m_best.value = current;
            m_best.loop = i;
            m_best.remaining_constraints = 0;

            best_recorder.try_update(x, current, i);
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

template<typename Solver, typename Float, typename Mode, typename Cost>
inline result
optimize_problem(const context_ptr& ctx, const problem& pb)
{
    result r;

    if (ctx->start)
        ctx->start(ctx->parameters);

    auto constraints{ make_merged_constraints(ctx, pb) };
    if (constraints.empty() || pb.vars.values.empty())
        return r;

    random_engine rng(init_random_generator_seed(ctx));

    auto variables = length(pb.vars.values);
    auto cost = Cost(pb.objective, variables);
    auto cost_constant = pb.objective.value;

    const auto thread = get_thread_number(ctx);

    std::vector<optimize_functor<Solver, Float, Mode, Cost>> functors;
    std::vector<std::thread> pool(thread);

    best_solution_recorder<Float, Mode> best_recorder(ctx);

    auto seeds = generate_seed(rng, thread);

    std::atomic_bool stop_task;
    stop_task.store(false);

    for (unsigned i = 0u; i != thread; ++i)
        functors.emplace_back(ctx, i, seeds[i]);

    for (unsigned i = 0u; i != thread; ++i)
        pool[i] = std::thread(std::ref(functors[i]),
                              std::ref(stop_task),
                              std::ref(best_recorder),
                              std::cref(constraints),
                              variables,
                              std::cref(cost),
                              cost_constant);

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

    r.strings = pb.strings;
    r.affected_vars = pb.affected_vars;
    r.variable_name = pb.vars.names;
    r.variables = variables;
    r.constraints = length(constraints);

    if (best_recorder.m_solutions.empty()) {
        r.status = result_status::time_limit_reached;
    } else if (!best_recorder.m_solutions.begin()->is_solution()) {
        r.status = result_status::time_limit_reached;

        r.duration = best_recorder.m_solutions.begin()->duration;
        r.loop = best_recorder.m_solutions.begin()->loop;
        r.remaining_constraints =
          best_recorder.m_solutions.begin()->remaining_constraints;
    } else {
        r.status = result_status::success;

        r.duration = best_recorder.m_solutions.begin()->duration;
        r.loop = best_recorder.m_solutions.begin()->loop;
        r.remaining_constraints =
          best_recorder.m_solutions.begin()->remaining_constraints;

        switch (ctx->parameters.storage) {
        case solver_parameters::storage_type::one: {
            r.solutions.resize(1);
            convert(
              *best_recorder.m_solutions.begin(), r.solutions[0], variables);
        } break;

        case solver_parameters::storage_type::bound: {
            if (best_recorder.m_solutions.size() == 1) {
                r.solutions.resize(1);
                convert(*best_recorder.m_solutions.begin(),
                        r.solutions[0],
                        variables);
            } else if (best_recorder.m_solutions.size() >= 2) {
                r.solutions.resize(2);
                convert(*best_recorder.m_solutions.begin(),
                        r.solutions[1],
                        variables);
                convert(*best_recorder.m_solutions.rbegin(),
                        r.solutions[0],
                        variables);
            }
        } break;

        case solver_parameters::storage_type::five: {
            const size_t max =
              std::min(size_t{ 5 }, std::size(best_recorder.m_solutions));
            r.solutions.resize(max);
            auto it = best_recorder.m_solutions.begin();
            for (size_t i{ 0 }; i != max; ++i)
                convert(*it++, r.solutions[max - i], variables);
        } break;
        }
    }

    if (ctx->finish)
        ctx->finish(r);

    return r;
}

} // namespace itm
} // namespace baryonyx

#endif
