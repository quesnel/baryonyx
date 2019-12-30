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

template<typename Cost, typename Mode>
class storage
{
private:
    mutable std::shared_mutex m_indices_mutex;
    using m_indices_writer = std::unique_lock<std::shared_mutex>;
    using m_indices_reader = std::shared_lock<std::shared_mutex>;

    mutable std::vector<std::shared_mutex> m_data_mutex;
    using m_data_writer = std::unique_lock<std::shared_mutex>;
    using m_data_reader = std::shared_lock<std::shared_mutex>;

    std::vector<int> m_indices;
    std::vector<raw_result<Mode>> m_data;
    bit_array m_bastert;

    std::normal_distribution<> m_choose_sol_dist;
    std::bernoulli_distribution m_crossover_bastert_insertion;

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
        m_data_writer lock{ m_data_mutex[id] };

        m_data[id].x = x;
        m_data[id].value = value;
        m_data[id].duration = duration;
        m_data[id].hash = hash;
        m_data[id].loop = loop;
        m_data[id].remaining_constraints = remaining_constraints;
    }

    int choose_a_solution(random_engine& rng)
    {
        double value;
        do
            value = m_choose_sol_dist(rng);
        while (value < 0 || value > 1);

        return static_cast<int>(m_size * value);
    }

public:
    storage(random_engine& rng,
            const Cost& costs_,
            const double cost_constant_,
            const int population_size,
            const int variables,
            const std::vector<merged_constraint>& constraints_,
            const double crossover_bastert_insertion,
            const double crossover_solution_selection_mean,
            const double crossover_solution_selection_stddev)
      : m_data_mutex(population_size)
      , m_indices(population_size)
      , m_data(population_size)
      , m_bastert(variables)
      , m_choose_sol_dist(crossover_solution_selection_mean,
                          crossover_solution_selection_stddev)
      , m_crossover_bastert_insertion(crossover_bastert_insertion)
      , costs(costs_)
      , cost_constant(cost_constant_)
      , m_size(population_size)
    {
        for (auto& elem : m_data)
            elem.x = bit_array(variables);

        init_with_bastert<Cost, Mode>(m_bastert, costs_, variables, 0);

        for (int i = 0, e = m_size / 2; i != e; ++i) {
            init_with_bastert<Cost, maximize_tag>(
              m_data[i].x, costs_, variables, 0);

            std::bernoulli_distribution dist(
              std::clamp(static_cast<double>(i) / (5 * e), 0.0, 1.0));

            for (int v = 0; v != variables; ++v)
                if (dist(rng))
                    m_data[i].x.invert(v);
        }

        for (int i = m_size / 2, e = m_size; i + 1 < e; i += 2) {
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

    void insert(const bit_array& x,
                const std::size_t hash,
                const int remaining_constraints,
                const double duration,
                const long int loop) noexcept
    {
        auto bound = indices_bound();

        to_log(stdout,
               5u,
               "- insert advance {} (hash: {}) {}s in {} loops\n",
               remaining_constraints,
               hash,
               duration,
               loop);

        to_log(stdout,
               5u,
               "- delete {} ({})\n",
               bound.last,
               m_data[bound.last].value);

        replace_result(bound.last,
                       x,
                       costs.results(x, cost_constant),
                       duration,
                       hash,
                       loop,
                       remaining_constraints);

        sort();
    }

    void insert(const bit_array& x,
                const std::size_t hash,
                const double value,
                const double duration,
                const long int loop) noexcept
    {
        auto bound = indices_bound();

        to_log(stdout,
               5u,
               "- insert solution {} (hash: {}) {}s in {} loops\n",
               value,
               hash,
               duration,
               loop);

        to_log(stdout,
               5u,
               "- delete {} ({})\n",
               bound.last,
               m_data[bound.last].value);

        replace_result(bound.last, x, value, duration, hash, loop, 0);

        sort();
    }

    bool can_be_inserted(const std::size_t hash, const int constraints) const
      noexcept
    {
        m_indices_reader lock(m_indices_mutex);

        for (int i = 0; i != m_size; ++i) {
            const auto id = m_indices[i];
            m_data_reader lock_data(m_data_mutex[id]);

            if (hash == m_data[id].hash)
                return false;

            if (m_data[id].remaining_constraints > constraints)
                return true;
        }

        return false;
    }

    bool can_be_inserted(const std::size_t hash, const double value) const
      noexcept
    {
        m_indices_reader lock(m_indices_mutex);

        for (int i = 0; i != m_size; ++i) {
            const auto id = m_indices[i];
            m_data_reader lock_data(m_data_mutex[id]);

            if (hash == m_data[id].hash)
                return false;

            if (m_data[id].remaining_constraints)
                return true;

            if (is_better_solution<Mode, double>(value, m_data[id].value))
                return true;
        }

        return false;
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
        int id;

        {
            m_indices_reader lock(m_indices_mutex);
            id = m_indices.front();
        }

        {
            m_data_reader lock_data(m_data_mutex[id]);
            constraints_remaining = m_data[id].remaining_constraints;
            value = m_data[id].value;
            duration = m_data[id].duration;
            loop = m_data[id].loop;
        }
    }

    const raw_result<Mode>& get_best(int i) const noexcept
    {
        return m_data[m_indices[i]];
    }

    void crossover(random_engine& rng,
                   bit_array& x,
                   const bit_array& first,
                   const bit_array& second)
    {
        std::uniform_int_distribution<bit_array::underlying_type> dist(0);
        std::bernoulli_distribution b_dist;

        const auto block_size = x.block_size();
        for (int i = 0; i != block_size; ++i) {
            const auto x1 = first.block(i);
            const auto x2 = second.block(i);
            const auto x_xor = x1 ^ x2;
            const auto x_rnd = dist(rng);
            const auto x_add = x_xor & x_rnd;

            x.set_block(i, x_add | (b_dist(rng) ? x1 : x2));
        }
    }

    void crossover(random_engine& rng, bit_array& x)
    {
        if (m_crossover_bastert_insertion(rng)) {
            int first = m_indices[choose_a_solution(rng)];

            m_data_reader lock_data_1{ m_data_mutex[first] };
            crossover(rng, x, m_data[first].x, m_bastert);

            to_log(stdout,
                   7u,
                   "- crossover between {} ({}) and bastert\n",
                   first,
                   m_data[first].value);

        } else {
            int first = m_indices[choose_a_solution(rng)];
            int second = m_indices[choose_a_solution(rng)];
            while (first == second)
                second = m_indices[choose_a_solution(rng)];

            m_data_reader lock_data_1{ m_data_mutex[first] };
            m_data_reader lock_data_2{ m_data_mutex[second] };
            crossover(rng, x, m_data[first].x, m_data[second].x);

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
              if (this->m_data[i1].remaining_constraints == 0) {
                  if (this->m_data[i2].remaining_constraints == 0)
                      return is_better_solution<Mode>(this->m_data[i1].value,
                                                      this->m_data[i2].value);
                  else
                      return true;
              } else {
                  if (this->m_data[i2].remaining_constraints == 0)
                      return false;
                  else {
                      if (this->m_data[i1].remaining_constraints ==
                          this->m_data[i2].remaining_constraints)
                          return is_better_solution<Mode>(
                            this->m_data[i1].value, this->m_data[i2].value);
                      else
                          return this->m_data[i1].remaining_constraints <
                                 this->m_data[i2].remaining_constraints;
                  }
              }
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
    enum init_status
    {
        init_with_best = 0,
        improve_best_1,
        improve_best_2,
        improve_best_3,
        improve_best_4,
        init_with_any,
        improve_any_1,
        improve_any_2,
        improve_any_3,
        improve_any_4,
        count
    };

    static init_status advance(init_status s) noexcept
    {
        auto current = static_cast<int>(s);
        const auto count = static_cast<int>(init_status::count);

        ++current;
        if (current >= count)
            current = 0;

        return static_cast<init_status>(current);
    }

    const context_ptr& m_ctx;
    std::chrono::time_point<std::chrono::steady_clock> m_start;
    std::vector<init_status> m_solver_state;

    storage<Cost, Mode> m_storage;

    std::normal_distribution<> m_variable_p_dist;
    std::normal_distribution<> m_value_p_dist;

    best_solution_recorder(const context_ptr& ctx,
                           const unsigned thread_number,
                           const Cost& costs,
                           const double cost_constant,
                           const int variables,
                           const std::vector<merged_constraint>& constraints,
                           random_engine& rng)
      : m_ctx(ctx)
      , m_solver_state(thread_number, init_status::init_with_best)
      , m_storage(rng,
                  costs,
                  cost_constant,
                  ctx->parameters.init_population_size,
                  variables,
                  constraints,
                  ctx->parameters.init_crossover_bastert_insertion,
                  ctx->parameters.init_crossover_solution_selection_mean,
                  ctx->parameters.init_crossover_solution_selection_stddev)
      , m_variable_p_dist(ctx->parameters.init_mutation_variable_mean,
                          ctx->parameters.init_mutation_variable_stddev)
      , m_value_p_dist(ctx->parameters.init_mutation_value_mean,
                       ctx->parameters.init_mutation_value_stddev)
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

    void mutation(random_engine& rng, bit_array& x)
    {
        m_storage.crossover(rng, x);

        if (m_value_p_dist.mean() == 0.0 && m_value_p_dist.stddev() == 0.0)
            return;

        double val_p, var_p;

        do
            var_p = m_variable_p_dist(rng);
        while (var_p <= 0.0 || var_p >= 1.0);

        do
            val_p = m_value_p_dist(rng);
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
            if (dist_var_p(rng)) {
                to_log(stdout, 9u, "- mutate variable {}\n", i);

                x.set(i, dist_value_p(rng));
            }
        }
    }

    double improve(const unsigned thread_id,
                   const double kappa_min,
                   const double kappa_max)
    {
        switch (m_solver_state[thread_id]) {
        case improve_best_2:
        case improve_any_2:
            return kappa_min + (kappa_max - kappa_min) * 0.1;
        case improve_best_3:
        case improve_any_3:
            return kappa_min + (kappa_max - kappa_min) * 0.15;
        case improve_best_4:
        case improve_any_4:
            return kappa_min + (kappa_max - kappa_min) * 0.20;
        default:
            return kappa_min;
        }
    }

    double reinit(random_engine& rng,
                  const unsigned thread_id,
                  const bool is_solution,
                  const double kappa_min,
                  const double kappa_max,
                  bit_array& x)
    {
        to_log(stdout,
               3u,
               "- reinitinialization thread {} - status {}. Is solution: {}\n",
               thread_id,
               m_solver_state[thread_id],
               is_solution);

        if (is_solution) {
            if (m_solver_state[thread_id] >= init_with_best &&
                m_solver_state[thread_id] < init_with_any)
                m_solver_state[thread_id] = init_with_any;
            else
                m_solver_state[thread_id] = init_with_best;
        }

        auto kappa = kappa_min;
        if (m_solver_state[thread_id] != init_with_best &&
            m_solver_state[thread_id] != init_with_any) {
            kappa = improve(thread_id, kappa_min, kappa_max);
            to_log(stdout, 5u, "- improve with kappa {}\n", kappa);
        } else {
            to_log(stdout, 5u, "- crossover and mutation\n");
            mutation(rng, x);
        }

        m_solver_state[thread_id] = advance(m_solver_state[thread_id]);

        return kappa;
    }

    void try_advance(const bit_array& solution,
                     const int remaining_constraints,
                     const long int loop)
    {
        auto hash = bit_array_hash()(solution);

        if (m_storage.can_be_inserted(hash, remaining_constraints)) {
            const auto end = std::chrono::steady_clock::now();
            const auto duration = compute_duration(m_start, end);

            m_storage.insert(
              solution, hash, remaining_constraints, duration, loop);
        }
    }

    void try_update(const bit_array& solution,
                    const double value,
                    const long int loop)
    {
        auto hash = bit_array_hash()(solution);

        if (m_storage.can_be_inserted(hash, value)) {
            const auto end = std::chrono::steady_clock::now();
            const auto duration = compute_duration(m_start, end);

            m_storage.insert(solution, hash, value, duration, loop);
        }
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
    const context_ptr& m_ctx;
    random_engine m_rng;
    long int m_call_number;
    unsigned m_thread_id;

    optimize_functor(const context_ptr& ctx,
                     unsigned thread_id,
                     typename random_engine::result_type seed)
      : m_ctx{ ctx }
      , m_rng{ seed }
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

        const auto w_limit = static_cast<long int>(p.w);

        if (p.limit <= 0)
            p.limit = std::numeric_limits<long int>::max();

        if (p.pushes_limit <= 0)
            p.pushes_limit = 0;

        if (p.pushing_iteration_limit <= 0)
            p.pushes_limit = 0;

        Solver slv(
          m_rng, length(constraints), variables, norm_costs, constraints);

        compute_order compute(p.order, variables);
        bool is_a_solution = false;

        while (!stop_task.load()) {
            ++m_call_number;
            const auto kappa_start = static_cast<Float>(best_recorder.reinit(
              m_rng, m_thread_id, is_a_solution, p.kappa_min, p.kappa_max, x));
            auto kappa = kappa_start;
            compute.init(slv, x);

            auto best_remaining = INT_MAX;
            is_a_solution = false;

            for (long int i = 0; !stop_task.load() && i != p.limit; ++i) {
                auto remaining =
                  compute.run(slv, x, m_rng, kappa, delta, theta);

                if (remaining == 0) {
                    best_recorder.try_update(
                      x, original_costs.results(x, cost_constant), i);
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
                best_recorder.try_advance(x, best_remaining, p.limit);
                continue;
            }

            for (int push = 0; !stop_task.load() && push < p.pushes_limit;
                 ++push) {

                auto remaining =
                  compute.push_and_run(slv,
                                       x,
                                       m_rng,
                                       pushing_k_factor,
                                       delta,
                                       theta,
                                       pushing_objective_amplifier);

                if (remaining == 0) {
                    best_recorder.try_update(
                      x,
                      original_costs.results(x, cost_constant),
                      -push * p.pushing_iteration_limit - 1);
                    is_a_solution = true;
                }

                kappa = kappa_start;
                for (int iter = 0;
                     !stop_task.load() && iter < p.pushing_iteration_limit;
                     ++iter) {

                    remaining =
                      compute.run(slv, x, m_rng, kappa, delta, theta);

                    if (remaining == 0) {
                        best_recorder.try_update(
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

    best_solution_recorder<Cost, Float, Mode> best_recorder(
      ctx, thread, cost, cost_constant, variables, constraints, rng);

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

        if (ctx->update) {
            auto call_number = 0L;
            for (auto i = 0u; i != thread; ++i)
                call_number += functors[i].m_call_number;

            int constraints_remaining;
            long int loop;
            double value;
            double duration;

            best_recorder.get_best(
              constraints_remaining, value, duration, loop);

            ctx->update(
              constraints_remaining, value, loop, duration, call_number);
        }

        end = std::chrono::steady_clock::now();
    } while (!is_time_limit(ctx->parameters.time_limit, start, end));

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

    switch (ctx->parameters.storage) {
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

    if (ctx->finish)
        ctx->finish(r);

    return r;
}

} // namespace itm
} // namespace baryonyx

#endif
