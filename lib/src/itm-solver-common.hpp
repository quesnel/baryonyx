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
#include <baryonyx/core-out>

#include <fmt/format.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <future>
#include <iterator>
#include <random>
#include <set>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include "branch-and-bound-solver.hpp"
#include "fixed_array.hpp"
#include "itm.hpp"
#include "knapsack-dp-solver.hpp"
#include "matrix.hpp"
#include "private.hpp"
#include "scoped_array.hpp"
#include "utils.hpp"

#include <cassert>

namespace baryonyx {
namespace itm {

struct bound
{
    bound() = default;

    bound(int min_, int max_)
      : min(min_)
      , max(max_)
    {}

    int min;
    int max;
};

template<typename floatingpointT>
struct r_data
{
    r_data() = default;

    r_data(floatingpointT value_, int index_)
      : value(value_)
      , id(index_)
    {}

    floatingpointT value; // reduced cost value
    int id;               // index in AP matrix
};

struct c_data
{
    c_data() = default;

    c_data(int id_A_, int id_r_)
      : id_A(id_A_)
      , id_r(id_r_)
    {}

    int id_A; // index in AP matrix
    int id_r; // index in r matrix
};

template<typename iteratorT, typename randomT>
inline void
random_shuffle_unique(iteratorT begin, iteratorT end, randomT& rng) noexcept
{
    auto ret = begin++;
    for (; begin != end; ++begin) {
        if (ret->value != begin->value) {
            std::shuffle(ret, begin, rng);
            ret = begin;
        }
    }

    std::shuffle(ret, begin, rng);
}

template<typename iteratorT, typename randomT>
inline void
calculator_sort(iteratorT begin, iteratorT end, randomT& rng, minimize_tag)
{
    if (std::distance(begin, end) > 1) {
        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return lhs.value < rhs.value;
        });

        random_shuffle_unique(begin, end, rng);
    }
}

template<typename iteratorT, typename randomT>
inline void
calculator_sort(iteratorT begin, iteratorT end, randomT& rng, maximize_tag)
{
    if (std::distance(begin, end) > 1) {
        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return rhs.value < lhs.value;
        });

        random_shuffle_unique(begin, end, rng);
    }
}

template<typename floatingpointT, typename randomT>
inline bool
stop_iterating(floatingpointT value, randomT& rng, minimize_tag) noexcept
{
    if (value == 0) {
        std::bernoulli_distribution d(0.5);
        return d(rng);
    }

    return value > 0;
}

template<typename floatingpointT, typename randomT>
inline bool
stop_iterating(floatingpointT value, randomT& rng, maximize_tag) noexcept
{
    if (value == 0) {
        std::bernoulli_distribution d(0.5);
        return d(rng);
    }

    return value < 0;
}

template<typename floatingpointT>
inline bool
is_better_solution(floatingpointT lhs,
                   floatingpointT rhs,
                   minimize_tag) noexcept
{
    return lhs < rhs;
}

template<typename floatingpointT>
inline bool
is_better_solution(floatingpointT lhs,
                   floatingpointT rhs,
                   maximize_tag) noexcept
{
    return lhs > rhs;
}

template<typename floatingpointT>
inline bool
init_x(floatingpointT cost, int value_if_cost_0, minimize_tag) noexcept
{
    if (cost < 0)
        return true;

    if (cost == 0)
        return value_if_cost_0;

    return false;
}

template<typename floatingpointT>
inline bool
init_x(floatingpointT cost, int value_if_cost_0, maximize_tag) noexcept
{
    if (cost > 0)
        return true;

    if (cost == 0)
        return value_if_cost_0;

    return false;
}

template<typename apT, typename xT, typename bT>
inline bool
is_valid_solution(const apT& ap, const xT& x, const bT& b) noexcept
{
    typename apT::const_iterator it, et;

    const auto& va = ap.A();

    for (int k = 0, ek = length(b); k != ek; ++k) {
        std::tie(it, et) = ap.row(k);
        int v = 0;

        for (; it != et; ++it)
            v += va[it->value] * x[it->position];

        if (not(b(k).min <= v and v <= b(k).max))
            return false;
    }

    return true;
}

template<typename apT, typename xT, typename bT, typename C>
inline int
compute_missing_constraint(const apT& ap,
                           const xT& x,
                           const bT& b,
                           C& r) noexcept
{
    typename apT::const_iterator it, et;

    const auto& va = ap.A();

    r.clear();

    for (int k = 0, ek = length(b); k != ek; ++k) {
        std::tie(it, et) = ap.row(k);
        int v = 0;

        for (; it != et; ++it)
            v += va[it->value] * x[it->position];

        if (not(b(k).min <= v and v <= b(k).max))
            r.emplace_back(k);
    }

    return length(r);
}

template<typename apT, typename xT, typename bT>
inline void
print_missing_constraint(const baryonyx::context_ptr& ctx,
                         const apT& ap,
                         const xT& x,
                         const bT& b,
                         const std::vector<std::string>& names) noexcept
{
    std::vector<int> R;

    compute_missing_constraint(ap, x, b, R);
    info(ctx, "Constraints remaining {}:\n", length(R));

    typename apT::const_iterator it, et;
    const auto& va = ap.A();

    for (auto k : R) {
        std::tie(it, et) = ap.row(k);
        int v = 0;

        info(ctx, "{}: {} <= ", k, b(k).min);

        for (; it != et; ++it) {
            v += va[it->value] * x[it->position];

            info(ctx,
                 "{:+d} [{} ({})] ",
                 va[it->value],
                 names[it->position],
                 x[it->position]);
        }

        info(ctx, " <= {} | value: {}\n", b(k).max, v);
    }
}

/**
 * Compute a problem lower or upper bounds based on Lagrangian multipliers
 * (valid if there are equality constraints only?)
 */
template<typename floatingpointT, typename modeT>
struct bounds_printer
{
    floatingpointT bestlb;
    floatingpointT bestub;
    floatingpointT max_cost;

    template<typename c_type>
    static floatingpointT max_cost_init(const c_type& c, minimize_tag) noexcept
    {
        assert(not c.empty());

        return *std::max_element(c.cbegin(), c.cend());
    }

    template<typename c_type>
    static floatingpointT max_cost_init(const c_type& c, maximize_tag) noexcept
    {
        assert(not c.empty());

        return *std::min_element(c.cbegin(), c.cend());
    }

    template<typename c_type>
    bounds_printer(const c_type& c)
      : bestlb(std::numeric_limits<floatingpointT>::lowest())
      , bestub(std::numeric_limits<floatingpointT>::max())
      , max_cost(max_cost_init(c, modeT()))
    {}

    template<typename SolverT>
    floatingpointT init_bound(const SolverT& slv, minimize_tag)
    {
        floatingpointT b{ 0 };

        for (auto c = 0; c != slv.m; ++c)
            b += slv.pi[c] * static_cast<floatingpointT>(slv.b(c).min);

        return b;
    }

    template<typename SolverT>
    floatingpointT init_bound(const SolverT& slv, maximize_tag)
    {
        floatingpointT b{ 0 };

        for (auto c = 0; c != slv.m; ++c)
            b += slv.pi[c] * static_cast<floatingpointT>(slv.b(c).max);

        return b;
    }

    template<typename SolverT>
    floatingpointT add_bound(const SolverT& slv,
                             int j,
                             floatingpointT sum_a_pi,
                             minimize_tag)
    {
        if (slv.c[j] - sum_a_pi < 0)
            return slv.c[j] - sum_a_pi;

        return { 0 };
    }

    template<typename SolverT>
    floatingpointT add_bound(const SolverT& slv,
                             int j,
                             floatingpointT sum_a_pi,
                             maximize_tag)
    {
        if (slv.c[j] - sum_a_pi > 0)
            return slv.c[j] - sum_a_pi;

        return { 0 };
    }

    void print_bound(const context_ptr& ctx,
                     floatingpointT lower_bound,
                     floatingpointT upper_bound,
                     minimize_tag)
    {
        bool better_gap = (lower_bound > bestlb || upper_bound < bestub);

        if (upper_bound < bestub)
            bestub = upper_bound;

        if (lower_bound > bestlb)
            bestlb = lower_bound;

        if (better_gap) {
            if (bestub == static_cast<floatingpointT>(0.0))
                info(ctx, "  - Lower bound: {}   (gap: 0%)\n", bestlb);
            else
                info(ctx,
                     "  - Lower bound: {}   (gap: {}%)\n",
                     bestlb,
                     static_cast<floatingpointT>(100.) * (bestub - bestlb) /
                       bestub);
        }
    }

    void print_bound(const context_ptr& ctx,
                     floatingpointT lower_bound,
                     floatingpointT upper_bound,
                     maximize_tag)
    {
        bool better_gap = (lower_bound > bestlb || upper_bound < bestub);

        if (upper_bound < bestub)
            bestub = upper_bound;

        if (lower_bound > bestlb)
            bestlb = lower_bound;

        if (better_gap) {
            if (bestlb == static_cast<floatingpointT>(0.0))
                info(ctx, "  - Upper bound: {}   (gap: 0%)\n", bestub);
            else
                info(ctx,
                     "  - Upper bound: {}   (gap: {}%)\n",
                     bestub,
                     static_cast<floatingpointT>(100.) * (bestlb - bestub) /
                       bestlb);
        }
    }

    floatingpointT init_ub(minimize_tag)
    {
        return std::numeric_limits<floatingpointT>::max();
    }

    floatingpointT init_ub(maximize_tag)
    {
        return std::numeric_limits<floatingpointT>::lowest();
    }

    template<typename SolverT>
    void operator()(const SolverT& slv,
                    const context_ptr& ctx,
                    const baryonyx::result& best)
    {
        using AP_type = typename SolverT::AP_type;
        using const_iterator = typename AP_type::const_iterator;

        floatingpointT lb = init_bound(slv, modeT());
        floatingpointT ub = init_ub(modeT());

        if (best)
            ub = static_cast<floatingpointT>(best.value);

        for (auto j = 0; j != slv.n; ++j) {
            floatingpointT sum_a_pi = 0.;

            const_iterator ht, hend;
            std::tie(ht, hend) = slv.ap.column(j);

            for (; ht != hend; ++ht) {
                auto a = slv.ap.A()[ht->value];
                sum_a_pi += static_cast<floatingpointT>(std::abs(a)) *
                            slv.pi[ht->position];
            }

            lb += add_bound(slv, j, sum_a_pi, modeT());
        }

        lb *= max_cost; // restore original cost

        print_bound(ctx, lb, ub, modeT());
    }
};

template<typename floatingpointT, typename randomT>
struct compute_none
{
    using random_type = randomT;

    const context_ptr& m_ctx;
    std::vector<int> R;

    template<typename solverT>
    compute_none(const context_ptr& ctx, const solverT& s, randomT&)
      : m_ctx(ctx)
      , R(s.m)
    {
        compute_missing_constraint(s.ap, s.x, s.b, R);
    }

    template<typename solverT>
    int push_and_run(solverT& solver,
                     floatingpointT kappa,
                     floatingpointT delta,
                     floatingpointT theta,
                     floatingpointT objective_amplifier)
    {
        for (int k = 0, e = solver.m; k != e; ++k)
            solver.push_and_compute_update_row(
              k, kappa, delta, theta, objective_amplifier);

        return compute_missing_constraint(solver.ap, solver.x, solver.b, R);
    }

    template<typename solverT>
    int run(solverT& solver,
            floatingpointT kappa,
            floatingpointT delta,
            floatingpointT theta)
    {
        for (auto it = R.begin(), et = R.end(); it != et; ++it)
            solver.compute_update_row(*it, kappa, delta, theta);

        return compute_missing_constraint(solver.ap, solver.x, solver.b, R);
    }
};

template<typename floatingpointT, typename randomT>
struct compute_reversing
{
    using random_type = randomT;

    const context_ptr& m_ctx;
    std::vector<int> R;
    int nb = 0;

    template<typename solverT>
    compute_reversing(const context_ptr& ctx, solverT& s, randomT&)
      : m_ctx(ctx)
      , R(s.m)
    {
        compute_missing_constraint(s.ap, s.x, s.b, R);
    }

    template<typename solverT>
    int push_and_run(solverT& solver,
                     floatingpointT kappa,
                     floatingpointT delta,
                     floatingpointT theta,
                     floatingpointT objective_amplifier)
    {
        for (int k = 0, e = solver.m; k != e; ++k)
            solver.push_and_compute_update_row(
              k, kappa, delta, theta, objective_amplifier);

        return compute_missing_constraint(solver.ap, solver.x, solver.b, R);
    }

    template<typename solverT>
    int run(solverT& solver,
            floatingpointT kappa,
            floatingpointT delta,
            floatingpointT theta)
    {
        for (auto it = R.rbegin(), et = R.rend(); it != et; ++it)
            solver.compute_update_row(*it, kappa, delta, theta);

        return compute_missing_constraint(solver.ap, solver.x, solver.b, R);
    }
};

template<typename floatingpointT, typename randomT>
struct compute_random
{
    using random_type = randomT;

    const context_ptr& m_ctx;
    std::vector<int> R;
    random_type& rng;

    template<typename solverT>
    compute_random(const context_ptr& ctx, solverT& s, random_type& rng_)
      : m_ctx(ctx)
      , R(s.m)
      , rng(rng_)
    {
        compute_missing_constraint(s.ap, s.x, s.b, R);
    }

    template<typename solverT>
    int push_and_run(solverT& solver,
                     floatingpointT kappa,
                     floatingpointT delta,
                     floatingpointT theta,
                     floatingpointT objective_amplifier)
    {
        for (int k = 0, e = solver.m; k != e; ++k)
            solver.push_and_compute_update_row(
              k, kappa, delta, theta, objective_amplifier);

        return compute_missing_constraint(solver.ap, solver.x, solver.b, R);
    }

    template<typename solverT>
    int run(solverT& solver,
            floatingpointT kappa,
            floatingpointT delta,
            floatingpointT theta)
    {
        std::shuffle(R.begin(), R.end(), rng);

        for (auto it = R.begin(), et = R.end(); it != et; ++it)
            solver.compute_update_row(*it, kappa, delta, theta);

        return compute_missing_constraint(solver.ap, solver.x, solver.b, R);
    }
};

struct compute_infeasibility_incr
{};

struct compute_infeasibility_decr
{};

template<typename iteratorT>
void
sort(iteratorT begin, iteratorT end, compute_infeasibility_incr)
{
    std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
        return lhs.second < rhs.second;
    });
}

template<typename iteratorT>
void
sort(iteratorT begin, iteratorT end, compute_infeasibility_decr)
{
    std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
        return rhs.second < lhs.second;
    });
}

template<typename floatingpointT, typename randomT, typename directionT>
struct compute_infeasibility
{
    using random_type = randomT;
    using direction_type = directionT;

    const context_ptr& m_ctx;
    std::vector<std::pair<int, int>> R;
    random_type& rng;

    template<typename solverT>
    compute_infeasibility(const context_ptr& ctx,
                          solverT& s,
                          random_type& rng_)
      : m_ctx(ctx)
      , R(s.m)
      , rng(rng_)
    {
        local_compute_missing_constraint(s);
    }

    template<typename solverT>
    int local_compute_missing_constraint(solverT& solver)
    {
        R.clear();

        for (int k = 0, e = solver.m; k != e; ++k) {
            auto ak{ solver.ap.row(k) };
            int v = 0;

            for (; std::get<0>(ak) != std::get<1>(ak); ++std::get<0>(ak))
                v += solver.ap.A()[std::get<0>(ak)->value] *
                     solver.x[std::get<0>(ak)->position];

            if (solver.b(k).min > v)
                R.push_back(std::make_pair(k, solver.b(k).min - v));
            else if (solver.b(k).max < v)
                R.push_back(std::make_pair(k, v - solver.b(k).max));
        }

        return length(R);
    }

    template<typename solverT>
    int push_and_run(solverT& solver,
                     floatingpointT kappa,
                     floatingpointT delta,
                     floatingpointT theta,
                     floatingpointT objective_amplifier)
    {
        for (int k = 0, e = solver.m; k != e; ++k)
            solver.push_and_compute_update_row(
              k, kappa, delta, theta, objective_amplifier);

        return local_compute_missing_constraint(solver);
    }

    template<typename solverT>
    int run(solverT& solver,
            floatingpointT kappa,
            floatingpointT delta,
            floatingpointT theta)
    {
        baryonyx::itm::sort(R.begin(), R.end(), direction_type());

        for (auto it = R.begin(), et = R.end(); it != et; ++it)
            solver.compute_update_row(it->first, kappa, delta, theta);

        return local_compute_missing_constraint(solver);
    }
};

template<typename CostT, typename floatingpointT>
inline floatingpointT
compute_delta(const context_ptr& ctx, const CostT& c, floatingpointT theta)
{
    info(ctx, "  - delta not defined, compute it:\n");

    auto mini{ std::numeric_limits<floatingpointT>::max() };
    for (int i = 0, e = numeric_cast<int>(c.size()); i != e; ++i)
        if (c[i] != 0 and std::abs(c[i]) < mini)
            mini = std::abs(c[i]);

    if (mini == std::numeric_limits<floatingpointT>::max()) {
        info(ctx, "    - All costs equal 0. Default value is used.\n");
        mini = static_cast<floatingpointT>(0.01);
    }

    const auto ret = mini - theta * mini;

    info(ctx,
         "    - delta={} (min normalized cost:{} / theta: {})\n",
         ret,
         mini,
         theta);

    return ret;
}

template<typename SolverT,
         typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT>
struct solver_functor
{
    using c_type = typename SolverT::c_type;

    using floatingpoint_type = floatingpointT;
    using mode_type = modeT;
    using constraint_order_type = constraintOrderT;
    using random_type = randomT;

    std::chrono::time_point<std::chrono::steady_clock> m_begin;
    std::chrono::time_point<std::chrono::steady_clock> m_end;

    const context_ptr& m_ctx;
    randomT& m_rng;
    const std::vector<std::string>& m_variable_names;
    const affected_variables& m_affected_vars;

    x_type m_best_x;
    result m_best;

    solver_functor(const context_ptr& ctx,
                   randomT& rng,
                   const std::vector<std::string>& variable_names,
                   const affected_variables& affected_vars)
      : m_ctx(ctx)
      , m_rng(rng)
      , m_variable_names(variable_names)
      , m_affected_vars(affected_vars)
    {}

    result operator()(const std::vector<itm::merged_constraint>& constraints,
                      int variables,
                      const c_type& original_costs,
                      const c_type& norm_costs,
                      double cost_constant,
                      const itm::parameters& p)
    {
        m_begin = std::chrono::steady_clock::now();
        m_end = m_begin;

        int i = 0;
        int pushed = -1;
        int best_remaining = -1;
        int pushing_iteration = p.pushing_iteration_limit;

        const auto kappa_min = static_cast<floatingpoint_type>(p.kappa_min);
        const auto kappa_step = static_cast<floatingpoint_type>(p.kappa_step);
        const auto kappa_max = static_cast<floatingpoint_type>(p.kappa_max);
        const auto alpha = static_cast<floatingpoint_type>(p.alpha);
        const auto theta = static_cast<floatingpoint_type>(p.theta);
        const auto delta = p.delta < 0
                             ? compute_delta<c_type, floatingpoint_type>(
                                 m_ctx, norm_costs, theta)
                             : static_cast<floatingpoint_type>(p.delta);

        const auto pushing_k_factor =
          static_cast<floatingpoint_type>(p.pushing_k_factor);
        const auto pushing_objective_amplifier =
          static_cast<floatingpoint_type>(p.pushing_objective_amplifier);

        auto kappa = kappa_min;

        SolverT slv(m_rng,
                    variables,
                    norm_costs,
                    constraints,
                    p.init_policy,
                    p.init_random);

        constraint_order_type compute(m_ctx, slv, m_rng);

        bounds_printer<floatingpointT, modeT> bound_print(original_costs);

        info(m_ctx, "* solver starts:\n");

        for (;;) {
            int remaining = compute.run(slv, kappa, delta, theta);

            if (best_remaining == -1 or remaining < best_remaining) {
                best_remaining = remaining;
                m_best = slv.results(original_costs, cost_constant);
                m_best.loop = i;
                m_best.remaining_constraints = remaining;
                m_best.duration =
                  std::chrono::duration_cast<std::chrono::duration<double>>(
                    m_end - m_begin)
                    .count();

                bound_print(slv, m_ctx, m_best);

                info(m_ctx,
                     "  - constraints remaining: {}/{} at {}s (loop: {})\n",
                     remaining,
                     m_best.constraints,
                     m_best.duration,
                     i);
            }

#ifndef BARYONYX_FULL_OPTIMIZATION
            slv.print(m_ctx, m_variable_names, p.print_level);
#endif

            if (m_best.status == result_status::success) {
                ++pushing_iteration;

                if (pushed == -1)
                    info(m_ctx, "  - start push system:\n");

                if (pushing_iteration >= p.pushing_iteration_limit) {
                    pushed++;
                    pushing_iteration = 0;

                    info(
                      m_ctx,
                      "    - push {}: kappa * k: {} objective amplifier: {}\n",
                      pushed,
                      (pushing_k_factor * kappa),
                      (pushing_objective_amplifier));

                    remaining =
                      compute.push_and_run(slv,
                                           pushing_k_factor * kappa,
                                           delta,
                                           theta,
                                           pushing_objective_amplifier);

                    if (remaining == 0) {
                        auto current =
                          slv.results(original_costs, cost_constant);
                        current.loop = i;
                        current.remaining_constraints = 0;

                        if (store_if_better(current))
                            m_best_x = slv.x;
                    }
                }

                if (pushed > p.pushes_limit) {
                    info(
                      m_ctx,
                      "    - Push system limit reached. Solution found: {}\n",
                      m_best.value);
                    return m_best;
                }
            }

            if (i > p.w)
                kappa += kappa_step *
                         std::pow(static_cast<floatingpointT>(remaining) /
                                    static_cast<floatingpointT>(slv.m),
                                  alpha);

            if (++i > p.limit) {
                info(m_ctx, "  - Loop limit reached: {}\n", i);
                if (pushed == -1)
                    m_best.status = result_status::limit_reached;

                if (context_get_integer_parameter(m_ctx, "print-level", 0) > 0)
                    print_missing_constraint(m_ctx,
                                             slv.ap,
                                             m_best.variable_value,
                                             slv.b,
                                             m_variable_names);

                return m_best;
            }

            if (kappa > kappa_max) {
                info(m_ctx, "  - Kappa max reached: {:+.6f}\n", kappa);
                if (pushed == -1)
                    m_best.status = result_status::kappa_max_reached;

                if (context_get_integer_parameter(m_ctx, "print-level", 0) > 0)
                    print_missing_constraint(m_ctx,
                                             slv.ap,
                                             m_best.variable_value,
                                             slv.b,
                                             m_variable_names);

                return m_best;
            }

            m_end = std::chrono::steady_clock::now();
            if (is_time_limit(p.time_limit, m_begin, m_end)) {
                info(m_ctx, "  - Time limit reached: {} {:+.6f}\n", i, kappa);
                if (pushed == -1)
                    m_best.status = result_status::time_limit_reached;

                if (context_get_integer_parameter(m_ctx, "print-level", 0) > 0)
                    print_missing_constraint(m_ctx,
                                             slv.ap,
                                             m_best.variable_value,
                                             slv.b,
                                             m_variable_names);

                return m_best;
            }
        }
    }

private:
    bool store_if_better(const result& current) noexcept
    {
        if (current.status != result_status::success)
            return false;

        if (m_best.status != result_status::success or
            is_better_solution(current.value, m_best.value, mode_type())) {

            double t =
              std::chrono::duration_cast<std::chrono::duration<double>>(
                m_end - m_begin)
                .count();

            info(m_ctx,
                 "  - Solution found: {} (i={} t={}s)\n",
                 current.value,
                 current.loop,
                 t);

            m_best = current;
            m_best.duration = t;

            std::ofstream ofs("temp.sol");
            ofs << m_best;

            std::size_t i, e;

            for (i = 0, e = m_affected_vars.names.size(); i != e; ++i)
                ofs << m_affected_vars.names[i] << ':'
                    << m_affected_vars.values[i] << '\n';

            for (i = 0, e = m_variable_names.size(); i != e; ++i)
                ofs << m_variable_names[i] << ':' << m_best.variable_value[i]
                    << '\n';

            return true;
        }

        return false;
    }
};

template<typename floatingpointT, typename modeT>
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

            if (m_best.status != result_status::success or
                is_better_solution(current.value, m_best.value, modeT())) {

                info(m_ctx,
                     "  - Solution found: {} (i={} t={}s)\n",
                     current.value,
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

template<typename SolverT,
         typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT>
struct optimize_functor
{
    using c_type = typename SolverT::c_type;

    using floatingpoint_type = floatingpointT;
    using mode_type = modeT;
    using constraint_order_type = constraintOrderT;
    using random_type = randomT;

    std::chrono::time_point<std::chrono::steady_clock> m_begin;
    std::chrono::time_point<std::chrono::steady_clock> m_end;

    const context_ptr& m_ctx;
    randomT m_rng;
    int m_thread_id;
    const std::vector<std::string>& m_variable_names;
    const affected_variables& m_affected_vars;
    x_type m_best_x;
    result m_best;

    optimize_functor(const context_ptr& ctx,
                     int thread_id,
                     typename random_type::result_type seed,
                     const std::vector<std::string>& variable_names,
                     const affected_variables& affected_vars)
      : m_ctx(ctx)
      , m_rng(seed)
      , m_thread_id(thread_id)
      , m_variable_names(variable_names)
      , m_affected_vars(affected_vars)
    {}

    result operator()(
      best_solution_recorder<floatingpointT, modeT>& best_recorder,
      const std::vector<itm::merged_constraint>& constraints,
      int variables,
      const c_type& original_costs,
      const c_type& norm_costs,
      double cost_constant,
      const itm::parameters& p)
    {
        m_begin = std::chrono::steady_clock::now();
        m_end = m_begin;

        int i = 0;
        int pushed = -1;
        int pushing_iteration = 0;

        const auto kappa_min = static_cast<floatingpoint_type>(p.kappa_min);
        const auto kappa_step = static_cast<floatingpoint_type>(p.kappa_step);
        const auto kappa_max = static_cast<floatingpoint_type>(p.kappa_max);
        const auto alpha = static_cast<floatingpoint_type>(p.alpha);
        const auto theta = static_cast<floatingpoint_type>(p.theta);
        const auto delta = p.delta < 0
                             ? compute_delta<c_type, floatingpoint_type>(
                                 m_ctx, norm_costs, theta)
                             : static_cast<floatingpoint_type>(p.delta);

        const auto pushing_k_factor =
          static_cast<floatingpoint_type>(p.pushing_k_factor);
        const auto pushing_objective_amplifier =
          static_cast<floatingpoint_type>(p.pushing_objective_amplifier);

        auto kappa = kappa_min;

        SolverT slv(m_rng,
                    variables,
                    norm_costs,
                    constraints,
                    p.init_policy,
                    p.init_random);

        constraint_order_type compute(m_ctx, slv, m_rng);

        for (; not is_time_limit(p.time_limit, m_begin, m_end);
             m_end = std::chrono::steady_clock::now(), ++i) {

            int remaining = compute.run(slv, kappa, delta, theta);

            if (remaining == 0) {
                auto current = slv.results(original_costs, cost_constant);
                current.loop = i;
                current.remaining_constraints = remaining;
                if (store_if_better(current)) {
                    m_best_x = slv.x;
                    pushed = 0;

                    if (best_recorder.try_update(m_best) and pushed > 0)
                        info(m_ctx, "   * pushed solution.\n");
                }
            }

            if (i > p.w)
                kappa += kappa_step *
                         std::pow(static_cast<floatingpoint_type>(remaining) /
                                    static_cast<floatingpoint_type>(slv.m),
                                  alpha);

            if (i >= p.limit or kappa > kappa_max or pushed > p.pushes_limit) {
                slv.reinit(m_best_x, p.init_policy, p.init_random);

                i = 0;
                kappa = static_cast<floatingpoint_type>(kappa_min);
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
                                           pushing_k_factor * kappa,
                                           delta,
                                           theta,
                                           pushing_objective_amplifier);

                    if (remaining == 0) {
                        auto current =
                          slv.results(original_costs, cost_constant);
                        current.loop = i;
                        current.remaining_constraints = 0;

                        if (store_if_better(current)) {
                            m_best_x = slv.x;

                            if (best_recorder.try_update(m_best))
                                info(m_ctx, "   * pushed solution.\n");
                        }
                    }
                }
            }
        }

        return m_best;
    }

private:
    bool store_if_better(const result& current) noexcept
    {
        if (current.status != result_status::success)
            return false;

        if (m_best.status != result_status::success or
            is_better_solution(current.value, m_best.value, mode_type())) {

            double t =
              std::chrono::duration_cast<std::chrono::duration<double>>(
                m_end - m_begin)
                .count();

            m_best = current;
            m_best.duration = t;

            std::ofstream ofs(fmt::format("temp-{}.sol", m_thread_id));
            ofs << m_best;

            std::size_t i, e;

            for (i = 0, e = m_affected_vars.names.size(); i != e; ++i)
                ofs << m_affected_vars.names[i] << ':'
                    << m_affected_vars.values[i] << '\n';

            for (i = 0, e = m_variable_names.size(); i != e; ++i)
                ofs << m_variable_names[i] << ':' << m_best.variable_value[i]
                    << '\n';

            return true;
        }

        return false;
    }
};

template<typename floatingpointT, typename iteratorT, typename randomT>
inline void
random_epsilon_unique(iteratorT begin,
                      iteratorT end,
                      randomT& rng,
                      floatingpointT min,
                      floatingpointT max)
{
    assert(min != max && "rng_normalize_cost fail to define min and max");

    std::uniform_real_distribution<floatingpointT> distribution(min, max);

    for (; begin != end; ++begin)
        begin->first += distribution(rng);
}

template<typename CostT, typename floatingpointT, typename randomT>
inline CostT
rng_normalize_costs(const CostT& c, randomT& rng)
{
    std::vector<std::pair<floatingpointT, int>> r(c.size());
    int i, e;

    for (i = 0, e = numeric_cast<int>(c.size()); i != e; ++i) {
        r[i].first = c[i];
        r[i].second = i;
    }

    std::sort(r.begin(), r.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
    });

    auto begin = r.begin();
    auto end = r.end();
    auto next = r.begin()++;
    for (; begin != end; ++begin) {
        if (next->first != begin->first) {
            floatingpointT min = next->first;
            floatingpointT max;

            if (begin != end)
                max = begin->first;
            else
                max = next->first + 1;

            random_epsilon_unique(next, begin, rng, min, max);
        }

        next = begin;
    }

    // Reorder the vector according to the variable index, so, it restores the
    // initial order.

    std::sort(r.begin(), r.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.second < rhs.second;
    });

    CostT ret(c);
    for (i = 0, e = numeric_cast<int>(c.size()); i != e; ++i)
        ret[i] = r[i].first;

    // Finally we compute the l+oo norm.

    floatingpointT div = *std::max_element(c.cbegin(), c.cend());
    if (std::isnormal(div))
        for (auto& elem : ret)
            elem /= div;

    return ret;
}

/**
 * Normalizes the cost vector, i.e. divides it by its l{1,2, +oo}norm. If the
 * input vector is too small or with infinity value, the c is unchanged.
 */
template<typename CostT, typename floatingpointT, typename randomT>
inline CostT
normalize_costs(const context_ptr& ctx,
                const std::string& norm,
                const CostT& c,
                randomT& rng)
{
    if (norm == "none") {
        info(ctx, "  - No norm");
        return c;
    }

    if (norm == "rng") {
        info(ctx, "  - Compute random norm\n");
        return rng_normalize_costs<CostT, floatingpointT, randomT>(c, rng);
    }

    CostT ret(c);
    floatingpointT div{ 0 };

    if (norm == "l1") {
        info(ctx, "  - Compute l1 norm\n");
        for (auto elem : ret)
            div += std::abs(elem);
    } else if (norm == "l2") {
        info(ctx, "  - Compute l2 norm\n");
        for (auto elem : ret)
            div += elem * elem;
    } else {
        info(ctx, "  - Compute infinity-norm (default)\n");
        div = *std::max_element(c.cbegin(), c.cend());
    }

    if (std::isnormal(div))
        for (auto& elem : ret)
            elem /= div;

    return ret;
}

template<typename CostT, typename floatingpointT>
inline CostT
make_objective_function(const objective_function& obj, int n)
{
    CostT ret(n, 0);

    for (const auto& elem : obj.elements)
        ret(elem.variable_index) += static_cast<floatingpointT>(elem.factor);

    return ret;
}

template<typename randomT>
typename randomT::result_type
init_random_generator_seed(const context_ptr& ctx)
{
    auto epoch = std::chrono::system_clock::now().time_since_epoch().count();
    auto param = context_get_integer_parameter(ctx, "seed", -1);

    if (param == -1)
        return static_cast<typename randomT::result_type>(epoch);

    return static_cast<typename randomT::result_type>(param);
}

template<typename SolverT,
         typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT>
inline result
solve_problem(const context_ptr& ctx, problem& pb, const itm::parameters& p)
{
    info(ctx, "Solver initializing\n");

    result ret;
    auto affected_vars = std::move(pb.affected_vars);

    auto constraints{ itm::make_merged_constraints(ctx, pb, p) };
    if (not constraints.empty() and not pb.vars.values.empty()) {

        randomT rng(init_random_generator_seed<randomT>(ctx));

        using CostT = typename SolverT::c_type;

        auto variables = numeric_cast<int>(pb.vars.values.size());
        auto cost = make_objective_function<CostT, floatingpointT>(
          pb.objective, variables);
        auto norm_costs = normalize_costs<CostT, floatingpointT, randomT>(
          ctx, p.norm, cost, rng);
        auto cost_constant = pb.objective.value;
        auto names = std::move(pb.vars.names);

        clear(pb);

        solver_functor<SolverT,
                       floatingpointT,
                       modeT,
                       constraintOrderT,
                       randomT>
          slv(ctx, rng, names, affected_vars);

        ret = slv(constraints, variables, cost, norm_costs, cost_constant, p);

        ret.method = "inequalities_Zcoeff solver";
        ret.variable_name = std::move(names);
    } else {
        ret.status = result_status::success;
    }
    ret.affected_vars = std::move(affected_vars);

    return ret;
}

/**
 * Generates an array with unique seed value.
 */
template<typename randomT>
inline scoped_array<typename randomT::result_type>
generate_seed(randomT& rng, int thread)
{
    using type = typename randomT::result_type;

    scoped_array<type> ret(thread);

    std::uniform_int_distribution<type> dst(std::numeric_limits<type>::min(),
                                            std::numeric_limits<type>::max());

    ret[0] = numeric_cast<type>(dst(rng));

    for (int i = 1; i != thread; ++i) {
        ret[i] = numeric_cast<type>(dst(rng));

        int j = i - 1;
        while (j > 0) {
            if (ret[j] == ret[i]) {
                ret[j] = numeric_cast<type>(dst(rng));
                j = i;
            }

            --j;
        }
    }

    return ret;
}

template<typename SolverT,
         typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT>
inline result
optimize_problem(const context_ptr& ctx,
                 problem& pb,
                 const itm::parameters& p,
                 int thread)
{
    info(ctx, "Optimizer initializing\n");
    result ret;
    auto affected_vars = std::move(pb.affected_vars);

    auto constraints{ itm::make_merged_constraints(ctx, pb, p) };
    if (not constraints.empty() and not pb.vars.values.empty()) {
        randomT rng(init_random_generator_seed<randomT>(ctx));

        using CostT = typename SolverT::c_type;

        auto variables = numeric_cast<int>(pb.vars.values.size());
        auto cost =
          make_objective_function<typename SolverT::c_type, floatingpointT>(
            pb.objective, variables);
        auto norm_costs = normalize_costs<CostT, floatingpointT, randomT>(
          ctx, p.norm, cost, rng);
        auto cost_constant = pb.objective.value;
        auto names = std::move(pb.vars.names);

        clear(pb);

        std::vector<std::thread> pool(thread);
        pool.clear();
        std::vector<std::future<result>> results(thread);
        results.clear();

        best_solution_recorder<floatingpointT, modeT> result(ctx);

        if (thread == 1)
            info(ctx, "optimizer starts with one thread\n");
        else
            info(ctx, "Optimizer starts with {} threads\n", thread);

        auto seeds = generate_seed(rng, thread);

        for (int i{ 0 }; i != thread; ++i) {
            std::packaged_task<baryonyx::result()> task(
              std::bind(optimize_functor<SolverT,
                                         floatingpointT,
                                         modeT,
                                         constraintOrderT,
                                         randomT>(
                          ctx, i, seeds[i], names, affected_vars),
                        std::ref(result),
                        std::ref(constraints),
                        variables,
                        std::ref(cost),
                        std::ref(norm_costs),
                        cost_constant,
                        std::ref(p)));

            results.emplace_back(task.get_future());

            pool.emplace_back(std::thread(std::move(task)));
        }

        for (auto& t : pool)
            t.join();

        ret = results[0].get();
        for (int i{ 1 }; i != thread; ++i) {
            auto current = results[i].get();
            if (current.status == result_status::success) {
                if (ret.status != result_status::success or
                    is_better_solution(current.value, ret.value, modeT()))
                    ret = current;
            }
        }

        ret.method = "inequalities_Zcoeff optimizer";
        ret.variable_name = std::move(names);
    } else {
        ret.status = result_status::success;
    }
    ret.affected_vars = std::move(affected_vars);

    return ret;
}
}
}

#endif
