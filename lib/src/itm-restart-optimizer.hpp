/* Copyright (C) 2016-2026 INRAE
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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_RESTART_OPTIMIZER_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_RESTART_OPTIMIZER_HPP

#include <algorithm>
#include <atomic>
#include <chrono>
#include <climits>
#include <cmath>
#include <memory>
#include <mutex>
#include <string_view>
#include <thread>
#include <tuple>
#include <vector>

#include "itm-common.hpp"
#include "itm-polish.hpp"
#include "itm-repair.hpp"
#include "result.hpp"
#include "utils.hpp"

namespace baryonyx::itm {

/** A change of one parameter tried by a restart thread. */
struct portfolio_variant
{
    std::string_view name;
    void (*apply)(solver_parameters& p) noexcept;
};

/**
 * The best parameters depend on the problem, as the parameter benchmark
 * showed: thread 0 keeps the user's parameters and each other thread changes
 * one of them, the most influential changes of the benchmark first.
 * Combining changes did worse than each change alone.
 */
inline constexpr portfolio_variant portfolio_variants[] = {
    { "norm=loo",
      [](solver_parameters& p) noexcept {
          p.cost_norm = solver_parameters::cost_norm_type::loo;
      } },
    { "constraint-order=random-sorting",
      [](solver_parameters& p) noexcept {
          p.order = solver_parameters::constraint_order::random_sorting;
      } },
    { "theta=0.8", [](solver_parameters& p) noexcept { p.theta = 0.8; } },
    { "delta=0.1", [](solver_parameters& p) noexcept { p.delta = 0.1; } },
    { "constraint-order=infeasibility-decr",
      [](solver_parameters& p) noexcept {
          p.order = solver_parameters::constraint_order::infeasibility_decr;
      } },
    { "norm=random",
      [](solver_parameters& p) noexcept {
          p.cost_norm = solver_parameters::cost_norm_type::random;
      } },
    { "kappa-min=0.1",
      [](solver_parameters& p) noexcept { p.kappa_min = 0.1; } },
};

/**
 * @return the index in @c portfolio_variants of the variant of
 *     @c thread_id, -1 for the user's parameters. The variants that leave
 *     the user's parameters unchanged are skipped.
 */
inline int
portfolio_variant_index(const solver_parameters& user, const unsigned thread_id)
{
    if (thread_id == 0)
        return -1;

    const auto changes = [&user](const portfolio_variant& variant) {
        auto p = user;
        variant.apply(p);
        return p.cost_norm != user.cost_norm || p.order != user.order ||
               p.theta != user.theta || p.delta != user.delta ||
               p.kappa_min != user.kappa_min;
    };

    std::vector<int> effective;
    for (int i = 0, e = length(portfolio_variants); i != e; ++i)
        if (changes(portfolio_variants[i]))
            effective.emplace_back(i);

    return effective.empty()
             ? -1
             : effective[(thread_id - 1u) % effective.size()];
}

/** Best solution shared by the restart threads. */
template<typename Mode>
class shared_best_solution
{
public:
    explicit shared_best_solution(int variables)
      : m_best(variables)
    {}

    void try_update_solution(const bit_array& x,
                             real value,
                             double duration,
                             long int loop,
                             unsigned thread_id,
                             bool by_jump)
    {
        std::scoped_lock lock(m_mutex);

        if (!m_best.is_solution()) {
            m_first_duration = duration;
            m_first_thread_id = thread_id;
            m_first_by_jump = by_jump;
        }

        if (m_best.is_other_better_solution(value)) {
            m_best.x = x;
            m_best.value = value;
            m_best.remaining_constraints = 0;
            m_best.duration = duration;
            m_best.loop = loop;
            m_thread_id = thread_id;
            m_by_jump = by_jump;
        }
    }

    bool has_solution() const
    {
        std::scoped_lock lock(m_mutex);

        return m_best.is_solution();
    }

    /** Only valid once every thread has stopped. */
    unsigned best_thread() const noexcept
    {
        return m_thread_id;
    }

    /** Only valid once every thread has stopped. */
    bool best_by_jump() const noexcept
    {
        return m_by_jump;
    }

    /** Only valid once every thread has stopped. */
    std::tuple<double, unsigned, bool> first_solution() const noexcept
    {
        return { m_first_duration, m_first_thread_id, m_first_by_jump };
    }

    void try_update_remaining(const bit_array& x,
                              int remaining_constraints,
                              double duration,
                              long int loop)
    {
        std::scoped_lock lock(m_mutex);

        if (m_best.is_other_better_solution(remaining_constraints)) {
            m_best.x = x;
            m_best.remaining_constraints = remaining_constraints;
            m_best.duration = duration;
            m_best.loop = loop;
        }
    }

    /** Copies the best feasible solution into @c x, if any. */
    bool copy_solution(bit_array& x) const
    {
        std::scoped_lock lock(m_mutex);

        if (!m_best.is_solution())
            return false;

        x = m_best.x;
        return true;
    }

    void get(int& remaining_constraints,
             real& value,
             double& duration,
             long int& loop) const
    {
        std::scoped_lock lock(m_mutex);

        remaining_constraints = m_best.remaining_constraints;
        value = m_best.value;
        duration = m_best.duration;
        loop = m_best.loop;
    }

    /** Only valid once every thread has stopped. */
    const raw_result<Mode>& best() const noexcept
    {
        return m_best;
    }

private:
    mutable std::mutex m_mutex;
    raw_result<Mode> m_best;
    unsigned m_thread_id = 0;
    bool m_by_jump = false;
    double m_first_duration = 0;
    unsigned m_first_thread_id = 0;
    bool m_first_by_jump = false;
};

/**
 * One thread of the restart optimizer: Wedelin runs from rotating initial
 * solutions, with a kappa step adapted to the outcome of each run.
 */
template<typename Solver, typename Mode, typename Cost>
struct restart_functor
{
    // A success slows the growth of kappa to improve the solutions, a
    // failure speeds it up to reach feasibility (success rule of evolution
    // strategies).
    static constexpr real kappa_step_on_success = 0.7;
    static constexpr real kappa_step_on_failure = 1.5;
    static constexpr real kappa_step_min = 1e-6;
    static constexpr real kappa_step_max = 0.1;

    // A copy: the parameters are those of the thread variant.
    context m_ctx;
    random_engine m_rng;
    unsigned m_thread_id;
    unsigned m_thread_number;
    int m_variant;

    // Searches a first solution with a weighted flip search before running
    // Wedelin: the search reaches feasibility on problems where Wedelin is
    // slow, and stops as soon as any thread has a solution.
    bool m_jump_first;

    // Read by the monitoring thread through std::atomic_ref.
    alignas(std::atomic_ref<long int>::required_alignment) long int
      m_call_number = 0;

    restart_functor(const context& ctx,
                    unsigned thread_id,
                    unsigned thread_number,
                    random_engine::result_type seed,
                    bool jump_first)
      : m_ctx(ctx)
      , m_rng(seed)
      , m_thread_id(thread_id)
      , m_thread_number(thread_number)
      , m_variant(ctx.tuning
                    ? -1
                    : portfolio_variant_index(ctx.parameters, thread_id))
      , m_jump_first(jump_first)
    {
        if (m_variant >= 0)
            portfolio_variants[m_variant].apply(m_ctx.parameters);
    }

    void operator()(const std::atomic_bool& stop_task,
                    shared_best_solution<Mode>& best,
                    const std::vector<merged_constraint>& constraints,
                    int variables,
                    const Cost& original_costs,
                    real cost_constant,
                    std::chrono::steady_clock::time_point start)
    {
        const auto& p = m_ctx.parameters;

        auto norm_costs =
          normalize_costs<Cost>(m_ctx, original_costs, m_rng, variables);

        const auto kappa_min = static_cast<real>(p.kappa_min);
        const auto kappa_max = static_cast<real>(p.kappa_max);
        const auto alpha = static_cast<real>(p.alpha);
        const auto theta = static_cast<real>(p.theta);
        const auto delta =
          p.delta < 0
            ? compute_delta<Cost>(m_ctx,
                                  norm_costs,
                                  static_cast<real>(p.kappa_step),
                                  theta,
                                  variables)
            : static_cast<real>(p.delta);
        const auto pushing_k_factor = static_cast<real>(p.pushing_k_factor);
        const auto pushing_objective_amplifier =
          static_cast<real>(p.pushing_objective_amplifier);
        const auto w_limit = warmup_loops(p.limit);

        // A run starting with a high kappa stays close to its initial
        // solution: drawing the starting kappa mixes local improvements of
        // the good solutions with deeper Lagrangian searches.
        uniform_real_distribution<double> kappa_start(
          solver_parameters::kappa_improve_start,
          solver_parameters::kappa_improve_stop);

        // Threads start from kappa steps spread around the requested one.
        const auto spread =
          (static_cast<real>(m_thread_id) -
           static_cast<real>(m_thread_number) / 2) /
          2;
        auto kappa_step = std::clamp(static_cast<real>(p.kappa_step) *
                                       std::pow(real{ 2 }, spread),
                                     kappa_step_min,
                                     kappa_step_max);

        Solver slv(
          m_rng, length(constraints), variables, norm_costs, constraints);
        compute_order compute(p.order, variables);
        bit_array x(variables);

        const auto elapsed = [start]() {
            return compute_duration(start, std::chrono::steady_clock::now());
        };

        // x stays the state of the Wedelin solver: only a copy is polished.
        polisher polish(constraints, variables);
        bit_array polished;
        repairer<Cost> repair(constraints, original_costs, variables);
        bit_array repaired;
        const auto record_solution = [&](long int loop, bool by_jump = false) {
            polished = x;
            polish(polished, original_costs);
            best.try_update_solution(polished,
                                     original_costs.results(polished,
                                                            cost_constant),
                                     elapsed(),
                                     loop,
                                     m_thread_id,
                                     by_jump);
        };

        if (m_jump_first) {
            weighted_flip_search<Cost> search(
              constraints, original_costs, variables);
            const auto flip_limit =
              repair_flip_limit(variables, length(constraints));
            const auto stop = [&]() {
                return stop_task.load() || best.has_solution();
            };

            for (unsigned attempt = 1; !stop(); ++attempt) {
                initial_solution(
                  x, best, norm_costs, constraints, attempt);

                if (search(x, m_rng, flip_limit, stop) == 0)
                    record_solution(attempt, true);
            }
        }

        const auto increase_kappa = [&](real& kappa, int remaining) {
            kappa += kappa_step * std::pow(static_cast<real>(remaining) /
                                             static_cast<real>(slv.m),
                                           alpha);
        };

        for (unsigned restart = m_thread_id; !stop_task.load(); ++restart) {
            std::atomic_ref<long int>(m_call_number)
              .fetch_add(1, std::memory_order_relaxed);

            initial_solution(x, best, norm_costs, constraints, restart);

            slv.reset();

            compute.init(slv, x);

            auto kappa = kappa_min + (kappa_max - kappa_min) *
                                       static_cast<real>(kappa_start(m_rng));
            auto remaining = INT_MAX;

            for (long int i = 0; !stop_task.load() && i != p.limit; ++i) {
                remaining =
                  compute.run(slv, x, m_rng, kappa, delta, theta);

                if (remaining == 0) {
                    record_solution(i);
                    break;
                }

                if (i > w_limit)
                    increase_kappa(kappa, remaining);

                if (kappa > kappa_max)
                    break;
            }

            if (remaining != 0 && remaining != INT_MAX &&
                is_repairable(remaining, slv.m)) {
                repaired = x;

                if (repair(repaired,
                           m_rng,
                           repair_flip_limit(variables, slv.m),
                           [&stop_task]() { return stop_task.load(); }) == 0) {
                    x = repaired;
                    remaining = 0;
                    record_solution(p.limit);
                }
            }

            if (remaining != 0) {
                if (remaining != INT_MAX)
                    best.try_update_remaining(
                      x, remaining, elapsed(), p.limit);

                kappa_step =
                  std::min(kappa_step * kappa_step_on_failure, kappa_step_max);
                continue;
            }

            kappa_step =
              std::max(kappa_step * kappa_step_on_success, kappa_step_min);

            for (int push = 0; !stop_task.load() && push < p.pushes_limit;
                 ++push) {
                remaining = compute.push_and_run(slv,
                                                 x,
                                                 m_rng,
                                                 pushing_k_factor * kappa,
                                                 delta,
                                                 theta,
                                                 pushing_objective_amplifier);

                if (remaining == 0)
                    record_solution(-push * p.pushing_iteration_limit - 1);

                for (int iter = 0;
                     !stop_task.load() && iter < p.pushing_iteration_limit;
                     ++iter) {
                    remaining =
                      compute.run(slv, x, m_rng, kappa, delta, theta);

                    if (remaining == 0) {
                        record_solution(-push * p.pushing_iteration_limit - iter - 1);
                        break;
                    }

                    increase_kappa(kappa, remaining);

                    if (kappa > kappa_max)
                        break;
                }
            }
        }
    }

private:
    // Rotates between the best solution found and the Bastert, pessimistic
    // and optimistic initializations, all slightly perturbed so that two
    // restarts never start from the same point.
    void initial_solution(bit_array& x,
                          const shared_best_solution<Mode>& best,
                          const Cost& costs,
                          const std::vector<merged_constraint>& constraints,
                          unsigned restart)
    {
        const int variables = x.size();
        bit_array empty;

        switch (restart % 4u) {
        case 0:
            if (best.copy_solution(x))
                break;
            [[fallthrough]];
        case 1:
            init_with_bastert<Cost, Mode>(x, costs, variables, 0);
            break;
        case 2:
            init_with_pre_solve<Cost, Mode>(
              x, empty, m_rng, costs, constraints, 1.0);
            break;
        default:
            init_with_pre_solve<Cost, Mode>(
              empty, x, m_rng, costs, constraints, 1.0);
            break;
        }

        uniform_real_distribution<double> rate(0.01, 0.05);
        bernoulli_distribution flip(rate(m_rng));

        for (int i = 0; i != variables; ++i)
            if (flip(m_rng))
                x.invert(i);
    }
};

/**
 * Optimizes with independent Wedelin restarts sharing only their best
 * solution: an alternative to the genetic optimizer of @c optimize_problem.
 */
template<typename Solver, typename Mode, typename Cost>
inline result
restart_optimize_problem(const context& ctx, const problem& pb)
{
    result r;

    if (ctx.start)
        ctx.start(ctx.parameters);

    auto constraints{ make_merged_constraints(ctx, pb) };
    if (constraints.empty() || pb.vars.values.empty()) {
        r.status = result_status::success;
        r.solutions.assign(1, unconstrained_solution(pb));
        r.strings = pb.strings;
        r.affected_vars = pb.affected_vars;
        r.variable_name = pb.vars.names;
        r.variables = length(pb.vars.values);
        r.remaining_constraints = 0;

        if (ctx.finish)
            ctx.finish(r);

        return r;
    }

    random_engine rng(init_random_generator_seed(ctx));

    const auto variables = length(pb.vars.values);
    const auto objective = minimization_objective(pb);
    const auto cost = Cost(objective, variables);
    const auto cost_constant = static_cast<real>(objective.value);

    const auto thread = get_thread_number(ctx);
    const auto seeds = generate_seed(rng, thread);

    // The sub-problems of lns run on 2 threads and are easy for Wedelin;
    // tuned evaluates exact configurations.
    const bool use_jump = thread >= 4u && !ctx.tuning;

    shared_best_solution<Mode> best(variables);
    std::atomic_bool stop_task{ false };

    std::vector<std::unique_ptr<restart_functor<Solver, Mode, Cost>>>
      functors;
    functors.reserve(thread);
    for (unsigned i = 0u; i != thread; ++i)
        functors.emplace_back(
          std::make_unique<restart_functor<Solver, Mode, Cost>>(
            ctx, i, thread, seeds[i], use_jump && i + 1u == thread));

    const auto start = std::chrono::steady_clock::now();

    std::vector<std::thread> pool;
    pool.reserve(thread);
    for (unsigned i = 0u; i != thread; ++i)
        pool.emplace_back(std::ref(*functors[i]),
                          std::cref(stop_task),
                          std::ref(best),
                          std::cref(constraints),
                          variables,
                          std::cref(cost),
                          cost_constant,
                          start);

    const auto limit =
      std::chrono::duration<double>(ctx.parameters.time_limit);
    auto end = start;

    do {
        const auto left = limit - (std::chrono::steady_clock::now() - start);
        std::this_thread::sleep_for(
          std::min<std::chrono::duration<double>>(std::chrono::seconds{ 1L },
                                                  left));

        long int call_number = 0;
        for (const auto& functor : functors)
            call_number += std::atomic_ref<long int>(functor->m_call_number)
                             .load(std::memory_order_relaxed);

        int remaining_constraints;
        long int loop;
        real value;
        double duration;
        best.get(remaining_constraints, value, duration, loop);

        if (ctx.update)
            ctx.update(remaining_constraints,
                       original_objective_value(pb, value),
                       loop,
                       duration,
                       call_number);

        end = std::chrono::steady_clock::now();
    } while (!is_time_limit(ctx.parameters.time_limit, start, end));

    stop_task.store(true);
    for (auto& t : pool)
        t.join();

    const auto& first = best.best();

    if (first.is_solution()) {
        const auto origin = [&](unsigned thread_id, bool by_jump) {
            const auto variant = functors[thread_id]->m_variant;
            return by_jump       ? std::string_view("feasibility jump")
                   : variant < 0 ? std::string_view("user parameters")
                                 : portfolio_variants[variant].name;
        };
        const auto [first_duration, first_thread, first_by_jump] =
          best.first_solution();

        notice(ctx,
               "  - restart optimizer: first solution at {:.3f}s by thread "
               "{} ({})\n",
               first_duration,
               first_thread,
               origin(first_thread, first_by_jump));
        notice(ctx,
               "  - restart optimizer: best solution found by thread {} ({})\n",
               best.best_thread(),
               origin(best.best_thread(), best.best_by_jump()));
    }

    r.status = first.is_solution() ? result_status::success
                                   : result_status::time_limit_reached;
    r.strings = pb.strings;
    r.affected_vars = pb.affected_vars;
    r.variable_name = pb.vars.names;
    r.variables = variables;
    r.constraints = length(constraints);
    r.duration = first.duration;
    r.loop = first.loop;
    r.remaining_constraints = first.remaining_constraints;

    r.solutions.resize(1);
    convert(first, r.solutions[0], variables);
    restore_objective_sign(pb, r);

    if (ctx.finish)
        ctx.finish(r);

    return r;
}

} // namespace baryonyx::itm

#endif
