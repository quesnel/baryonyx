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

template<typename Cost, typename Mode>
void
init_with_bastert(bit_array& x,
                  const Cost& c,
                  const int variables,
                  const int value_if_0) noexcept
{
    for (int i = 0; i != variables; ++i)
        if (init_x<Mode>(c[i], value_if_0))
            x.set(i);
        else
            x.unset(i);
}

inline void
init_with_random(bit_array& x,
                 random_engine& rng,
                 const int variables,
                 const double init_ramdom) noexcept
{
    std::bernoulli_distribution dist(init_ramdom);

    for (int i = 0; i != variables; ++i)
        if (dist(rng))
            x.set(i);
        else
            x.unset(i);
}

template<typename Cost, typename Mode>
void
init_with_pre_solve(bit_array& x_pessimistic,
                    bit_array& x_optimistic,
                    random_engine& rng,
                    const Cost& c,
                    const std::vector<merged_constraint>& constraints) noexcept
{
    int max_length = 0;
    for (const auto& cst : constraints)
        if (length(cst.elements) > max_length)
            max_length = length(cst.elements);

    struct reduced_cost
    {
        float value;
        int factor;
        int id;
    };

    std::vector<reduced_cost> R(max_length);

    for (const auto& cst : constraints) {
        R.resize(cst.elements.size());
        const int r_size = length(cst.elements);

        for (int i = 0; i != r_size; ++i) {
            R[i].value = static_cast<float>(c[cst.elements[i].variable_index]);
            R[i].factor = cst.elements[i].factor;
            R[i].id = cst.elements[i].variable_index;
        }

        std::shuffle(std::begin(R), std::end(R), rng);
        std::sort(
          std::begin(R), std::end(R), [](const auto& lhs, const auto& rhs) {
              return is_better_solution<Mode>(lhs.value, rhs.value);
          });

        {
            int sum = 0;
            int best = -2;

            for (int i = -1; i < r_size; ++i) {
                if (cst.min <= sum && sum <= cst.max) {
                    best = i;
                    break;
                }

                if (i + 1 < r_size)
                    sum += R[i + 1].factor;
            }

            int i = 0;
            for (; i <= best; ++i)
                x_pessimistic.set(R[i].id);

            for (; i != r_size; ++i) {
                x_pessimistic.unset(R[i].id);
            }
        }

        {
            int sum = 0;
            int best = -2;

            for (int i = -1; i < r_size; ++i) {
                if (cst.min <= sum && sum <= cst.max)
                    best = i;

                if (best != -2 && i + 1 < r_size &&
                    stop_iterating<Mode>(R[i + 1].value))
                    break;

                if (i + 1 < r_size)
                    sum += R[i + 1].factor;
            }

            int i = 0;
            for (; i <= best; ++i)
                x_optimistic.set(R[i].id);

            for (; i != r_size; ++i)
                x_optimistic.unset(R[i].id);
        }
    }
}

template<typename Cost, typename Mode>
struct storage
{
    std::vector<raw_result<Mode>> m_data;
    std::vector<int> m_indices;
    const Cost& costs;
    double cost_constant;

    storage(random_engine& rng,
            const Cost& costs_,
            const double cost_constant_,
            const int population_size,
            const int variables,
            const std::vector<merged_constraint>& constraints_)
      : m_indices(population_size)
      , costs(costs_)
      , cost_constant(cost_constant_)
    {
        m_data.reserve(population_size);
        for (int i = 0; i < population_size; ++i)
            m_data.emplace_back(variables);

        init_with_bastert<Cost, Mode>(m_data[0].x, costs, variables, 0);
        init_with_bastert<Cost, Mode>(m_data[1].x, costs, variables, 1);

        init_with_pre_solve<Cost, Mode>(
          m_data[2].x, m_data[3].x, rng, costs, constraints_);

        for (int i = 4; i != population_size; ++i) {
            auto pc = static_cast<double>(i) /
                      static_cast<double>(population_size - 4);

            if (pc > 1.0)
                pc -= 1.0;

            init_with_random(m_data[i].x, rng, variables, pc);
        }

        std::iota(std::begin(m_indices), std::end(m_indices), 0);

        for (int i = 0; i != population_size; ++i) {
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

        sort();
    }

    void sort() noexcept
    {
        std::sort(
          std::begin(m_indices), std::end(m_indices), [this](int i1, int i2) {
              if (this->m_data[i1].remaining_constraints == 0) {
                  if (this->m_data[i2].remaining_constraints == 0)
                      return is_better_solution<Mode>(this->m_data[i1].value,
                                                      this->m_data[i2].value);
                  else
                      return true;
              } else {
                  return this->m_data[i1].remaining_constraints <
                         this->m_data[i2].remaining_constraints;

                  // if (this->m_data[i2].remaining_constraints == 0)
                  //     return false;
                  // else
                  //     return
                  //     is_better_solution<Mode>(this->m_data[i1].value,
                  //                                     this->m_data[i2].value);
              }
          });

#ifdef BARYONYX_ENABLE_DEBUG
        to_log(stdout, 2u, "Solutions init population:\n");
        for (int i = 0, e = length(m_indices); i != e; ++i)
            to_log(stdout,
                   2u,
                   "- {} id {} value {} constraint {} hash {}\n",
                   i,
                   m_indices[i],
                   m_data[m_indices[i]].value,
                   m_data[m_indices[i]].remaining_constraints,
                   m_data[m_indices[i]].hash);
#endif
    }

    bool insert(const bit_array& x,
                const std::size_t hash,
                const int remaining_constraints,
                const double duration,
                const index loop) noexcept
    {
        to_log(stdout,
               5u,
               "insert advance {} (hash: {}) {}s in {} loops\n",
               remaining_constraints,
               hash,
               duration,
               loop);

        m_data[m_indices.back()].x = x;
        m_data[m_indices.back()].value = costs.results(x, cost_constant);
        m_data[m_indices.back()].duration = duration;
        m_data[m_indices.back()].hash = hash;
        m_data[m_indices.back()].loop = loop;
        m_data[m_indices.back()].remaining_constraints = remaining_constraints;

        auto old_remaining = m_data[m_indices.front()].remaining_constraints;
        auto old_value = m_data[m_indices.front()].value;

        sort();

        return old_remaining !=
                 m_data[m_indices.front()].remaining_constraints ||
               old_value != m_data[m_indices.front()].value;
    }

    bool insert(const bit_array& x,
                const std::size_t hash,
                const double value,
                const double duration,
                const index loop) noexcept
    {
        to_log(stdout,
               5u,
               "insert solution {} (hash: {}) {}s in {} loops\n",
               value,
               hash,
               duration,
               loop);

        m_data[m_indices.back()].x = x;
        m_data[m_indices.back()].value = value;
        m_data[m_indices.back()].duration = duration;
        m_data[m_indices.back()].hash = hash;
        m_data[m_indices.back()].loop = loop;
        m_data[m_indices.back()].remaining_constraints = 0;

        auto old_remaining = m_data[m_indices.front()].remaining_constraints;
        auto old_value = m_data[m_indices.front()].value;

        sort();

        return old_remaining !=
                 m_data[m_indices.front()].remaining_constraints ||
               old_value != m_data[m_indices.front()].value;
    }

    bool already_exists(std::size_t hash) const noexcept
    {
        return std::end(m_data) !=
               std::find_if(std::begin(m_data),
                            std::end(m_data),
                            [hash](const auto& r) { return r.hash == hash; });
    }

    bool can_be_inserted(std::size_t hash, int constraints) const noexcept
    {
        to_log(stdout,
               5u,
               "constraint can_be_inserted: {} - {}\n",
               hash,
               !already_exists(hash));

        return m_data[m_indices.back()].remaining_constraints > constraints &&
               !already_exists(hash);
    }

    bool can_be_inserted(std::size_t hash, double value) const noexcept
    {
        to_log(stdout,
               5u,
               "solution can_be_inserted: {} - {}\n",
               hash,
               !already_exists(hash));

        if (m_data[m_indices.back()].remaining_constraints)
            return true;

        if (already_exists(hash))
            return false;

        return is_better_solution<Mode, double>(
          value, m_data[m_indices.back()].value);
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

    const context_ptr& m_ctx;
    std::chrono::time_point<std::chrono::steady_clock> m_start;
    std::vector<int> m_solver_state;
    mutable std::mutex m_mutex;
    random_engine& rng;

    storage<Cost, Mode> m_storage;
    std::uniform_int_distribution<int> split_at_dist;
    std::uniform_int_distribution<int> choose_sol_dist;
    std::bernoulli_distribution choose_best_or_any;
    std::bernoulli_distribution choose_mutation;

    std::atomic_int m_remaining_constraints;
    double m_value;

    best_solution_recorder(const context_ptr& ctx,
                           const unsigned thread_number,
                           const Cost& costs,
                           const double cost_constant,
                           const int variables,
                           const std::vector<merged_constraint>& constraints,
                           random_engine& rng_)
      : m_ctx(ctx)
      , m_solver_state(thread_number, 0)
      , rng(rng_)
      , m_storage(rng_,
                  costs,
                  cost_constant,
                  ctx->parameters.init_population_size,
                  variables,
                  constraints)
      , split_at_dist(1, variables - 1)
      , choose_sol_dist(1, ctx->parameters.init_population_size - 1)
      , choose_best_or_any(ctx->parameters.init_policy_random)
      , choose_mutation(ctx->parameters.init_random)
      , m_remaining_constraints(
          m_storage.m_data[m_storage.m_indices.front()].remaining_constraints)
      , m_value(m_storage.m_data[m_storage.m_indices.front()].value)
    {
        m_start = std::chrono::steady_clock::now();
    }

    void mutation(const unsigned thread_id, bit_array& x)
    {
        int first, second;
        const int split_at = split_at_dist(rng);

        if (m_solver_state[thread_id] < init_with_any) {
            first = choose_sol_dist(rng) / 10;
            second = choose_sol_dist(rng) / 5;
            if (second == first)
                ++second;
        } else {
            first = choose_sol_dist(rng);
            second = choose_sol_dist(rng);
        }

        try {
            std::lock_guard<std::mutex> lock(m_mutex);
            if (choose_best_or_any(rng)) {
                x.assign(m_storage.m_data[first].x, 0, split_at);
                x.assign(m_storage.m_data[second].x, split_at, x.size());
            } else {
                x.assign(m_storage.m_data[second].x, 0, split_at);
                x.assign(m_storage.m_data[first].x, split_at, x.size());
            }
        } catch (const std::exception& e) {
            error(m_ctx, "sync optimization error: {}", e.what());
        }

        for (int i = 0, e = x.size(); i != e; ++i)
            if (choose_mutation(rng))
                x.invert(i);

        ++m_solver_state[thread_id];
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

    double reinit(const unsigned thread_id,
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
            mutation(thread_id, x);
            to_log(stdout, 5u, "- mutation\n");
        }

        ++m_solver_state[thread_id];
        m_solver_state[thread_id] %= (count - 1);

        return kappa;
    }

    void try_advance(const bit_array& solution,
                     const int remaining_constraints,
                     const int loop)
    {
        try {
            auto hash = bit_array_hash()(solution);

            std::lock_guard<std::mutex> lock(m_mutex);

            if (m_storage.can_be_inserted(hash, remaining_constraints)) {
                const auto end = std::chrono::steady_clock::now();
                const auto duration = compute_duration(m_start, end);

                if (m_storage.insert(
                      solution, hash, remaining_constraints, duration, loop))
                    m_ctx->update(remaining_constraints, 0.0, loop, duration);
            }
        } catch (const std::exception& e) {
            error(m_ctx, "sync optimization error: {}", e.what());
        }
    }

    void try_update(const bit_array& solution,
                    const double value,
                    const int loop)
    {
        try {
            auto hash = bit_array_hash()(solution);

            std::lock_guard<std::mutex> lock(m_mutex);

            if (m_storage.can_be_inserted(hash, value)) {
                const auto end = std::chrono::steady_clock::now();
                const auto duration = compute_duration(m_start, end);

                if (m_storage.insert(solution, hash, value, duration, loop))
                    m_ctx->update(0, value, loop, duration);
            }
        } catch (const std::exception& e) {
            error(m_ctx, "sync optimization error: {}", e.what());
        }
    }

    bool have_solution() const noexcept
    {
        return m_storage.m_data[m_storage.m_indices.front()].is_solution();
    }
};

template<typename Solver, typename Float, typename Mode, typename Cost>
struct optimize_functor
{
    const context_ptr& m_ctx;
    random_engine m_rng;
    int m_thread_id;

    optimize_functor(const context_ptr& ctx,
                     unsigned thread_id,
                     typename random_engine::result_type seed)
      : m_ctx(ctx)
      , m_rng(seed)
      , m_thread_id(thread_id)
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

        const int w_limit = static_cast<int>(p.w);

        if (p.limit <= 0)
            p.limit = std::numeric_limits<int>::max();

        if (p.pushes_limit <= 0)
            p.pushes_limit = 0;

        if (p.pushing_iteration_limit <= 0)
            p.pushes_limit = 0;

        Solver slv(
          m_rng, length(constraints), variables, norm_costs, constraints);

        compute_order compute(p.order, variables);
        compute.init(slv, x);
        bool is_a_solution = false;

        while (!stop_task.load()) {
            auto kappa = static_cast<Float>(best_recorder.reinit(
              m_thread_id, is_a_solution, p.kappa_min, p.kappa_max, x));

            auto best_remaining = INT_MAX;
            is_a_solution = false;

            for (int i = 0; !stop_task.load() && i != p.limit; ++i) {
                auto remaining =
                  compute.run(slv, x, m_rng, kappa, delta, theta);

                if (remaining == 0) {
                    best_recorder.try_update(
                      x, original_costs.results(x, cost_constant), i);
                    best_remaining = 0;
                    is_a_solution = true;
                    break;
                }

                if (remaining < best_remaining)
                    best_remaining = remaining;

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
                                       pushing_k_factor * kappa,
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

    const auto& first = best_recorder.m_storage
                          .m_data[best_recorder.m_storage.m_indices.front()];

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

        int last = 1;
        const int end = length(best_recorder.m_storage.m_indices);
        while (last < end && best_recorder.m_storage
                               .m_data[best_recorder.m_storage.m_indices[last]]
                               .is_solution())
            ++last;

        convert(first, r.solutions[0], variables);
        convert(best_recorder.m_storage
                  .m_data[best_recorder.m_storage.m_indices[last]],
                r.solutions[1],
                variables);
    } break;

    case solver_parameters::storage_type::five: {
        r.solutions.resize(5);

        for (int i = 0; i != 5; ++i) {
            const auto& elem = best_recorder.m_storage
                                 .m_data[best_recorder.m_storage.m_indices[i]];
            convert(elem, r.solutions[i], variables);
        }
    } break;
    }

    if (ctx->finish)
        ctx->finish(r);

    return r;
}

} // namespace itm
} // namespace baryonyx

#endif
