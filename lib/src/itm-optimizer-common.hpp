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
#include <shared_mutex>
#include <thread>
#include <utility>
#include <vector>

#include "itm-common.hpp"
#include "private.hpp"
#include "result.hpp"
#include "sparse-matrix.hpp"
#include "utils.hpp"

namespace baryonyx {
namespace itm {

/**
 * @brief Local context for each thread.
 * @details Each thread need a copy of this data to avoid mutex.
 *
 */
struct local_context
{
    random_engine rng;
    std::normal_distribution<double> choose_sol_dist;
    std::normal_distribution<double> variable_p_dist;
    std::normal_distribution<double> value_p_dist;
    std::uniform_int_distribution<int> bad_solution_choose;
    std::uniform_int_distribution<std::size_t> crossover_dist;
    std::bernoulli_distribution crossover_bastert_insertion;

    const double init_kappa_improve_start;
    const double init_kappa_improve_increase;
    const double init_kappa_improve_stop;

    const unsigned thread_id;

    local_context(const context& ctx,
                  const unsigned thread_id_,
                  const random_engine::result_type seed)
      : rng(seed)
      , choose_sol_dist(
          ctx.parameters.init_crossover_solution_selection_mean,
          ctx.parameters.init_crossover_solution_selection_stddev)
      , variable_p_dist(ctx.parameters.init_mutation_variable_mean,
                        ctx.parameters.init_mutation_variable_stddev)
      , value_p_dist(ctx.parameters.init_mutation_value_mean,
                     ctx.parameters.init_mutation_value_stddev)
      , bad_solution_choose(ctx.parameters.init_population_size / 5,
                            ctx.parameters.init_population_size - 1)
      , crossover_dist(0)
      , crossover_bastert_insertion(
          ctx.parameters.init_crossover_bastert_insertion)
      , init_kappa_improve_start(ctx.parameters.init_kappa_improve_start)
      , init_kappa_improve_increase(ctx.parameters.init_kappa_improve_increase)
      , init_kappa_improve_stop(ctx.parameters.init_kappa_improve_stop)
      , thread_id(thread_id_)
    {}
};

template<typename Cost, typename Mode>
class storage
{
private:
    mutable std::shared_mutex m_indices_mutex;
    using m_indices_writer = std::unique_lock<std::shared_mutex>;
    using m_indices_reader = std::shared_lock<std::shared_mutex>;

    std::vector<int> m_indices;
    std::vector<raw_result<Mode>> m_data;
    bit_array m_bastert;
    bit_array m_random;

    const Cost& costs;
    double cost_constant;
    const int m_size;

    struct bound_indices
    {
        bound_indices(int first_, int last_)
          : first(first_)
          , last(last_)
        {}

        int first;
        int last;
    };

    bound_indices indices_bound() const noexcept
    {
        m_indices_reader lock{ m_indices_mutex };

        return bound_indices{ m_indices.front(), m_indices.back() };
    }

    void replace_result(const int id,
                        const bit_array& x,
                        const double value,
                        const double duration,
                        const std::size_t hash,
                        const long int loop,
                        const int remaining_constraints) noexcept
    {
        m_indices_reader lock{ m_indices_mutex };

        m_data[id].x = x;
        m_data[id].value = value;
        m_data[id].duration = duration;
        m_data[id].hash = hash;
        m_data[id].loop = loop;
        m_data[id].remaining_constraints = remaining_constraints;
    }

    int choose_a_bad_solution(local_context& ctx)
    {
        return ctx.bad_solution_choose(ctx.rng);
    }

    int choose_a_solution(local_context& ctx)
    {
        double value;
        do
            value = ctx.choose_sol_dist(ctx.rng);
        while (value < 0 || value > 1);

        return static_cast<int>(m_size * value);
    }

public:
    storage(random_engine& rng,
            const Cost& costs_,
            const double cost_constant_,
            const int population_size,
            const int variables,
            const std::vector<merged_constraint>& constraints_)
      : m_indices(population_size)
      , m_data(population_size)
      , m_bastert(variables)
      , m_random(variables)
      , costs(costs_)
      , cost_constant(cost_constant_)
      , m_size(population_size)
    {
        for (auto& elem : m_data)
            elem.x = bit_array(variables);

        init_with_bastert<Cost, Mode>(m_bastert, costs_, variables, 0);

        for (int i = 0, e = m_size / 2; i != e; ++i) {
            m_data[i].x = m_bastert;

            std::bernoulli_distribution dist(
              std::clamp(static_cast<double>(i) / (5 * e), 0.0, 1.0));

            for (int v = 0; v != variables; ++v)
                if (dist(rng))
                    m_data[i].x.invert(v);
        }

        for (int i = m_size / 2, e = m_size; i + 1 < e; i += 2) {
            init_with_random(m_data[i].x, rng, variables, 0.2);
            init_with_random(m_data[i + 1].x, rng, variables, 0.8);

            init_with_pre_solve<Cost, Mode>(
              m_data[i].x,
              m_data[i + 1].x,
              rng,
              costs_,
              constraints_,
              std::clamp(static_cast<double>(i) / (5 * e), 0.0, 1.0));
        }

        for (int i = 0, e = m_size; i != e; ++i) {
            m_data[i].make_hash();
            m_data[i].value = costs.results(m_data[i].x, cost_constant_);
            m_data[i].remaining_constraints = 0;

            for (int k = 0, end_k = length(constraints_); k != end_k; ++k) {
                int sum = 0;
                for (const auto& elem : constraints_[k].elements) {
                    if (m_data[i].x[elem.variable_index])
                        sum += elem.factor;
                }

                if (!(constraints_[k].min <= sum &&
                      sum <= constraints_[k].max))
                    ++m_data[i].remaining_constraints;
            }
        }

        std::iota(std::begin(m_indices), std::end(m_indices), 0);

        sort();
    }

    void show_population(const context& ctx) const
    {
        info(ctx, " Population {}:\n", m_indices.size());
        for (int i = 0; i != m_size; ++i)
            info(ctx,
                 "  - {}: value {} constraints {} hash {}\n",
                 m_indices[i],
                 m_data[m_indices[i]].value,
                 m_data[m_indices[i]].remaining_constraints,
                 m_data[m_indices[i]].hash);
    }

    void insert(local_context& ctx,
                const bit_array& x,
                const std::size_t hash,
                const int remaining_constraints,
                const double duration,
                const long int loop) noexcept
    {
        to_log(stdout,
               5u,
               "- insert advance {} (hash: {}) {}s in {} loops\n",
               remaining_constraints,
               hash,
               duration,
               loop);

        int id_to_delete = m_indices[choose_a_bad_solution(ctx)];

        to_log(stdout,
               5u,
               "- delete {} ({})\n",
               id_to_delete,
               m_data[id_to_delete].value);

        replace_result(id_to_delete,
                       x,
                       costs.results(x, cost_constant),
                       duration,
                       hash,
                       loop,
                       remaining_constraints);

        sort();
    }

    void insert(local_context& ctx,
                const bit_array& x,
                const std::size_t hash,
                const double value,
                const double duration,
                const long int loop) noexcept
    {
        to_log(stdout,
               5u,
               "- insert solution {} (hash: {}) {}s in {} loops\n",
               value,
               hash,
               duration,
               loop);

        int id_to_delete = m_indices[choose_a_bad_solution(ctx)];

        to_log(stdout,
               5u,
               "- delete {} ({})\n",
               id_to_delete,
               m_data[id_to_delete].value);

        replace_result(id_to_delete, x, value, duration, hash, loop, 0);

        sort();
    }

    bool can_be_inserted(const std::size_t hash, const int constraints) const
      noexcept
    {
        m_indices_reader lock(m_indices_mutex);

        for (int i = 0; i != m_size; ++i)
            if (m_data[i].remaining_constraints == constraints &&
                m_data[i].hash == hash)
                return false;

        return true;
    }

    bool can_be_inserted([[maybe_unused]] const std::size_t hash,
                         const double value) const noexcept
    {
        m_indices_reader lock(m_indices_mutex);

        for (int i = 0; i != m_size; ++i)
            if (m_data[i].remaining_constraints == 0 &&
                m_data[i].value == value && m_data[i].hash == hash)
                return false;

        return true;
    }

    const raw_result<Mode>& get_worst() const noexcept
    {
        int last = 1;
        const int end = length(m_indices);

        while (last != end && m_data[m_indices[last]].is_solution())
            ++last;

        return last == end ? m_data[m_indices[last - 1]]
                           : m_data[m_indices[last]];
    }

    void get_best(int& constraints_remaining,
                  double& value,
                  double& duration,
                  long int& loop) const noexcept
    {
        m_indices_reader lock(m_indices_mutex);
        int id = m_indices.front();

        constraints_remaining = m_data[id].remaining_constraints;
        value = m_data[id].value;
        duration = m_data[id].duration;
        loop = m_data[id].loop;
    }

    const raw_result<Mode>& get_best(int i) const noexcept
    {
        return m_data[m_indices[i]];
    }

    void crossover(local_context& ctx,
                   bit_array& x,
                   const bit_array& first,
                   const bit_array& second)
    {
        const auto block_size = x.block_size();
        for (int i = 0; i != block_size; ++i) {
            const auto x1 = first.block(i);
            const auto x2 = second.block(i);

            x.set_block(i, ((x1 ^ x2) & ctx.crossover_dist(ctx.rng)) ^ x1);
        }
    }

    void crossover(local_context& ctx, bit_array& x)
    {
        m_indices_reader lock(m_indices_mutex);

        if (ctx.crossover_bastert_insertion(ctx.rng)) {
            int first = m_indices[choose_a_solution(ctx)];
            std::bernoulli_distribution b(0.5);

            if (b(ctx.rng)) {
                if (b(ctx.rng)) {
                    crossover(ctx, x, m_data[first].x, m_bastert);
                    to_log(stdout,
                           7u,
                           "- crossover between {} ({}) and bastert\n",
                           first,
                           m_data[first].value);
                } else {
                    m_data[first].x = m_bastert;
                }
            } else {
                init_with_random(m_random, ctx.rng, x.size(), 0.5);
                if (b(ctx.rng)) {
                    crossover(ctx, x, m_data[first].x, m_random);
                    to_log(stdout,
                           7u,
                           "- crossover between {} ({}) and random\n",
                           first,
                           m_data[first].value);
                } else {
                    m_data[first].x = m_random;
                }
            }
        } else {
            int first = m_indices[choose_a_solution(ctx)];
            int second = m_indices[choose_a_solution(ctx)];
            while (first == second)
                second = m_indices[choose_a_solution(ctx)];

            crossover(ctx, x, m_data[first].x, m_data[second].x);

            to_log(stdout,
                   7u,
                   "- crossover between {} ({}) and {} ({})\n",
                   first,
                   m_data[first].value,
                   second,
                   m_data[second].value);
        }
    }

private:
    void sort() noexcept
    {
        m_indices_writer lock{ m_indices_mutex };

        std::sort(
          std::begin(m_indices), std::end(m_indices), [this](int i1, int i2) {
              const int cst_1 = this->m_data[i1].remaining_constraints;
              const int cst_2 = this->m_data[i2].remaining_constraints;
              const double value_1 = this->m_data[i1].value;
              const double value_2 = this->m_data[i2].value;

              if (cst_1 < cst_2)
                  return true;

              if (cst_1 == cst_2)
                  return is_better_solution<Mode>(value_1, value_2);

              return false;
          });

#ifdef BARYONYX_ENABLE_DEBUG
        to_log(stdout, 3u, "- Solutions init population:\n");
        for (int i = 0; i != m_size; ++i) {
            to_log(stdout,
                   5u,
                   "- {} id {} value {} constraint {} hash {}\n",
                   i,
                   m_indices[i],
                   m_data[m_indices[i]].value,
                   m_data[m_indices[i]].remaining_constraints,
                   m_data[m_indices[i]].hash);
        }
#endif
    }
};

template<typename Cost, typename Float, typename Mode>
struct best_solution_recorder
{
    std::chrono::time_point<std::chrono::steady_clock> m_start;
    storage<Cost, Mode> m_storage;
    std::vector<double> m_kappa_append;

    best_solution_recorder(random_engine& rng,
                           const unsigned thread_number,
                           const Cost& costs,
                           const double cost_constant,
                           const int variables,
                           const std::vector<merged_constraint>& constraints,
                           const int population_size)
      : m_storage{ rng,       costs,      cost_constant, population_size,
                   variables, constraints }
      , m_kappa_append(std::size_t{ thread_number }, 0.0)
    {
        m_start = std::chrono::steady_clock::now();
    }

    void get_best(int& constraints_remaining,
                  double& value,
                  double& duration,
                  long int& loop) const noexcept
    {
        m_storage.get_best(constraints_remaining, value, duration, loop);
    }

    void crossover(local_context& ctx, bit_array& x)
    {
        m_storage.crossover(ctx, x);
    }

    void mutation(local_context& ctx, bit_array& x)
    {
        if (ctx.value_p_dist.mean() == 0.0 && ctx.value_p_dist.stddev() == 0.0)
            return;

        double val_p, var_p;

        do
            var_p = ctx.variable_p_dist(ctx.rng);
        while (var_p <= 0.0 || var_p >= 1.0);

        do
            val_p = ctx.value_p_dist(ctx.rng);
        while (val_p < 0.0 || val_p > 1.0);

        to_log(stdout,
               7u,
               "- mutation variables {}% with "
               " {}% of set\n",
               var_p,
               val_p);

        std::bernoulli_distribution dist_var_p(var_p);
        std::bernoulli_distribution dist_value_p(val_p);

        for (int i = 0, e = x.size(); i != e; ++i) {
            if (dist_var_p(ctx.rng)) {
                to_log(stdout, 9u, "- mutate variable {}\n", i);

                x.set(i, dist_value_p(ctx.rng));
            }
        }
    }

    double reinit(local_context& ctx,
                  const bool /*is_solution*/,
                  const double kappa_min,
                  const double kappa_max,
                  bit_array& x)
    {
        to_log(stdout, 3u, "- reinitinialization thread {}.\n", ctx.thread_id);

        double kappa = kappa_min;

        if (m_kappa_append[ctx.thread_id] < ctx.init_kappa_improve_stop) {
            m_kappa_append[ctx.thread_id] += ctx.init_kappa_improve_increase;
            kappa = kappa_min +
                    (kappa_max - kappa_min) * m_kappa_append[ctx.thread_id];

            to_log(stdout, 5u, "- improve with kappa {}\n", kappa);
        } else {
            m_kappa_append[ctx.thread_id] = ctx.init_kappa_improve_start;
            crossover(ctx, x);

            to_log(stdout, 5u, "- crossover\n");
        }

        mutation(ctx, x);

        return kappa;
    }

    void try_advance(local_context& ctx,
                     const bit_array& solution,
                     const int remaining_constraints,
                     const long int loop)
    {
        auto hash = bit_array_hash()(solution);

        if (m_storage.can_be_inserted(hash, remaining_constraints)) {
            const auto end = std::chrono::steady_clock::now();
            const auto duration = compute_duration(m_start, end);

            m_storage.insert(
              ctx, solution, hash, remaining_constraints, duration, loop);
        }
    }

    void try_update(local_context& ctx,
                    const bit_array& solution,
                    const double value,
                    const long int loop)
    {
        auto hash = bit_array_hash()(solution);

        if (m_storage.can_be_inserted(hash, value)) {
            const auto end = std::chrono::steady_clock::now();
            const auto duration = compute_duration(m_start, end);

            m_storage.insert(ctx, solution, hash, value, duration, loop);
        }
    }

    void show_population(const context& ctx) const noexcept
    {
        m_storage.show_population(ctx);
    }

    const raw_result<Mode>& get_worst() const noexcept
    {
        return m_storage.get_worst();
    }

    const raw_result<Mode>& get_best(int i) const noexcept
    {
        return m_storage.get_best(i);
    }
};

template<typename Solver, typename Float, typename Mode, typename Cost>
struct optimize_functor
{
    const context& m_ctx;
    local_context m_local_ctx;
    long int m_call_number;
    unsigned m_thread_id;

    optimize_functor(const context& ctx,
                     const unsigned thread_id,
                     const typename random_engine::result_type seed)
      : m_ctx{ ctx }
      , m_local_ctx{ ctx, thread_id, seed }
      , m_call_number{ 0 }
      , m_thread_id{ thread_id }
    {}

    void operator()(const std::atomic_bool& stop_task,
                    best_solution_recorder<Cost, Float, Mode>& best_recorder,
                    const std::vector<merged_constraint>& constraints,
                    int variables,
                    const Cost& original_costs,
                    double cost_constant)
    {
        bit_array x(variables);

        auto& p = m_ctx.parameters;

        auto norm_costs = normalize_costs<Float, Cost>(
          m_ctx, original_costs, m_local_ctx.rng, variables);

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

        const auto w_limit = static_cast<long int>(p.w);

        Solver slv(m_local_ctx.rng,
                   length(constraints),
                   variables,
                   norm_costs,
                   constraints);

        compute_order compute(p.order, variables);
        bool is_a_solution = false;

        while (!stop_task.load()) {
            ++m_call_number;
            const auto kappa_start = static_cast<Float>(best_recorder.reinit(
              m_local_ctx, is_a_solution, p.kappa_min, p.kappa_max, x));
            auto kappa = kappa_start;
            compute.init(slv, x);

            auto best_remaining = INT_MAX;
            is_a_solution = false;

            for (long int i = 0; !stop_task.load() && i != p.limit; ++i) {
                auto remaining =
                  compute.run(slv, x, m_local_ctx.rng, kappa, delta, theta);

                if (remaining == 0) {
                    best_recorder.try_update(
                      m_local_ctx,
                      x,
                      original_costs.results(x, cost_constant),
                      i);

                    best_remaining = 0;
                    is_a_solution = true;
                    break;
                } else {
                    best_remaining = std::min(remaining, best_remaining);
                }

                if (i > w_limit)
                    kappa +=
                      kappa_step * std::pow(static_cast<Float>(remaining) /
                                              static_cast<Float>(slv.m),
                                            alpha);

                if (kappa > kappa_max)
                    break;
            }

            if (best_remaining > 0) {
                best_recorder.try_advance(
                  m_local_ctx, x, best_remaining, p.limit);
                continue;
            }

            for (int push = 0; !stop_task.load() && push < p.pushes_limit;
                 ++push) {

                auto remaining =
                  compute.push_and_run(slv,
                                       x,
                                       m_local_ctx.rng,
                                       pushing_k_factor,
                                       delta,
                                       theta,
                                       pushing_objective_amplifier);

                if (remaining == 0) {
                    best_recorder.try_update(
                      m_local_ctx,
                      x,
                      original_costs.results(x, cost_constant),
                      -push * p.pushing_iteration_limit - 1);
                    is_a_solution = true;
                }

                kappa = kappa_start;
                for (int iter = 0;
                     !stop_task.load() && iter < p.pushing_iteration_limit;
                     ++iter) {

                    remaining = compute.run(
                      slv, x, m_local_ctx.rng, kappa, delta, theta);

                    if (remaining == 0) {
                        best_recorder.try_update(
                          m_local_ctx,
                          x,
                          original_costs.results(x, cost_constant),
                          -push * p.pushing_iteration_limit - iter - 1);
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
};

//
// Get number of thread to use in optimizer from parameters list or
// from the standard thread API. If an error occurred, this function
// returns 1.
//
inline unsigned
get_thread_number(const baryonyx::context& ctx) noexcept
{
    unsigned ret;

    if (ctx.parameters.thread <= 0)
        ret = std::thread::hardware_concurrency();
    else
        ret = static_cast<unsigned>(ctx.parameters.thread);

    if (ret == 0)
        return 1;

    return ret;
}

template<typename Solver, typename Float, typename Mode, typename Cost>
inline result
optimize_problem(const context& ctx, const problem& pb)
{
    result r;

    if (ctx.start)
        ctx.start(ctx.parameters);

    auto constraints{ make_merged_constraints(ctx, pb) };
    if (constraints.empty() || pb.vars.values.empty()) {
        r.status = result_status::success;
        r.solutions.resize(1);
        r.solutions.back().value = pb.objective.value;
        r.strings = pb.strings;
        r.affected_vars = pb.affected_vars;
        r.variable_name = pb.vars.names;
        return r;
    }

    random_engine rng(init_random_generator_seed(ctx));

    auto variables = length(pb.vars.values);
    auto cost = Cost(pb.objective, variables);
    auto cost_constant = pb.objective.value;

    const auto thread = get_thread_number(ctx);

    std::vector<optimize_functor<Solver, Float, Mode, Cost>> functors;
    std::vector<std::thread> pool(thread);

    best_solution_recorder<Cost, Float, Mode> best_recorder(
      rng,
      thread,
      cost,
      cost_constant,
      variables,
      constraints,
      ctx.parameters.init_population_size);

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

    const auto start = std::chrono::steady_clock::now();
    auto end = start;

    do {
        std::this_thread::sleep_for(std::chrono::seconds{ 1L });

        if (ctx.update) {
            auto call_number = 0L;
            for (auto i = 0u; i != thread; ++i)
                call_number += functors[i].m_call_number;

            int constraints_remaining;
            long int loop;
            double value;
            double duration;

            best_recorder.get_best(
              constraints_remaining, value, duration, loop);

            ctx.update(
              constraints_remaining, value, loop, duration, call_number);
        }

        end = std::chrono::steady_clock::now();
    } while (!is_time_limit(ctx.parameters.time_limit, start, end));

    stop_task.store(true);

    for (auto& t : pool)
        t.join();

    r.strings = pb.strings;
    r.affected_vars = pb.affected_vars;
    r.variable_name = pb.vars.names;
    r.variables = variables;
    r.constraints = length(constraints);

    const auto& first = best_recorder.get_best(0);

    if (!first.is_solution())
        r.status = result_status::time_limit_reached;
    else
        r.status = result_status::success;

    r.duration = first.duration;
    r.loop = first.loop;
    r.remaining_constraints = first.remaining_constraints;

    switch (ctx.parameters.storage) {
    case solver_parameters::storage_type::one: {
        r.solutions.resize(1);
        convert(first, r.solutions[0], variables);
    } break;

    case solver_parameters::storage_type::bound: {
        r.solutions.resize(2);

        convert(first, r.solutions[0], variables);
        convert(best_recorder.get_worst(), r.solutions[1], variables);
    } break;

    case solver_parameters::storage_type::five: {
        r.solutions.resize(5);

        for (int i = 0; i != 5; ++i)
            convert(best_recorder.get_best(i), r.solutions[i], variables);
    } break;
    }

    best_recorder.show_population(ctx);

    if (ctx.finish)
        ctx.finish(r);

    return r;
}

} // namespace itm
} // namespace baryonyx

#endif
