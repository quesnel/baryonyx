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

#include "fixed_array.hpp"
#include "itm-common.hpp"
#include "itm.hpp"
#include "private.hpp"
#include "sparse-matrix.hpp"
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

    floatingpointT value; ///< Reduced cost value.
    int id;               ///< Index in ap.row() vector.
};

struct c_data
{
    c_data() = default;

    c_data(int id_r_)
      : id_r(id_r_)
    {}

    int id_r; ///< Index in ap.row() vector.
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
stop_iterating(floatingpointT value, minimize_tag) noexcept
{
    return value > 0;
}

template<typename floatingpointT>
inline bool
stop_iterating(floatingpointT value, maximize_tag) noexcept
{
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

template<typename TimePoint>
double
compute_duration(const TimePoint& first, const TimePoint& last) noexcept
{
    namespace sc = std::chrono;

    auto diff = last - first;
    auto dc = sc::duration_cast<sc::duration<double, std::ratio<1>>>(diff);

    return dc.count();
}

inline std::size_t
compute_reduced_costs_vector_size(
  const std::vector<itm::merged_constraint>& csts) noexcept
{
    std::size_t rsizemax = csts[0].elements.size();

    for (std::size_t i = 1, e = csts.size(); i != e; ++i)
        rsizemax = std::max(rsizemax, csts[i].elements.size());

    return rsizemax;
}

template<typename Solver>
inline bool
is_valid_solution(const Solver& s) noexcept
{
    return s.is_valid_solution();
}

template<typename Solver, typename Container>
inline int
compute_missing_constraint(const Solver& s, Container& c)
{
    return s.compute_violated_constraints(c);
}

template<typename Solver, typename x_type>
inline itm::init_policy_type
init_solver(Solver& slv,
            const x_type& best_previous,
            itm::init_policy_type type,
            double init_random)
{
    using floatingpointT = typename Solver::floatingpoint_type;

    std::fill(slv.P.get(),
              slv.P.get() + slv.ap.length(),
              static_cast<floatingpointT>(0));
    std::fill(
      slv.pi.get(), slv.pi.get() + slv.m, static_cast<floatingpointT>(0));

    if (best_previous.empty() and type == itm::init_policy_type::best)
        type = itm::init_policy_type::bastert;

    init_random = clamp(init_random, 0.0, 1.0);
    std::bernoulli_distribution d(init_random);

    switch (type) {
    case itm::init_policy_type::bastert:
        if (init_random == 0.0 or init_random == 1.0) {
            bool value_if_cost_0 = init_random == 1.0;

            for (int i = 0; i != slv.n; ++i)
                slv.x[i] = init_x(
                  slv.c[i], value_if_cost_0, typename Solver::mode_type());
        } else {
            for (int i = 0; i != slv.n; ++i)
                slv.x[i] =
                  init_x(slv.c[i], d(slv.rng), typename Solver::mode_type());
        }

        type = itm::init_policy_type::random;
        break;
    case itm::init_policy_type::random:
        for (int i = 0; i != slv.n; ++i)
            slv.x[i] = d(slv.rng);

        type = itm::init_policy_type::best;
        break;
    case itm::init_policy_type::best:
        for (int i = 0; i != slv.n; ++i)
            slv.x[i] = (d(slv.rng)) ? (best_previous[i]) : (!best_previous[i]);

        type = itm::init_policy_type::bastert;
        break;
    }

    return type;
}

template<typename Solver>
inline void
print_solver(const Solver& slv,
             const context_ptr& ctx,
             const std::vector<std::string>& names,
             int print_level)
{
    if (print_level <= 0)
        return;

    debug(ctx, "  - X: {} to {}\n", 0, slv.n);
    for (int i = 0; i != slv.n; ++i)
        debug(ctx,
              "    - {} {}={}/c_i:{}\n",
              i,
              names[i],
              (slv.x[i] ? 1 : 0),
              slv.c[i]);
    debug(ctx, "\n");

    for (int k = 0; k != slv.m; ++k) {
        sparse_matrix<int>::const_row_iterator it, et;

        std::tie(it, et) = slv.ap.row(k);
        int v = 0;

        for (; it != et; ++it)
            v += slv.factor(it->value) * slv.x[it->column];

        bool valid = slv.bound_min(k) <= v and v <= slv.bound_max(k);

        debug(ctx,
              "C {}:{} (Lmult: {})\n",
              k,
              (valid ? "   valid" : "violated"),
              slv.pi[k]);
    }
}

template<typename Solver>
inline void
print_missing_constraint(const baryonyx::context_ptr& ctx,
                         const Solver& slv,
                         const std::vector<std::string>& names) noexcept
{
    std::vector<int> R;

    compute_missing_constraint(slv, R);
    info(ctx, "Constraints remaining: {}\n", length(R));

    sparse_matrix<int>::const_row_iterator it, et;

    for (auto k : R) {
        std::tie(it, et) = slv.ap.row(k);
        int v = 0;

        info(ctx, "{}: {} <= ", k, slv.bound_min(k));

        for (; it != et; ++it) {
            v += slv.factor(it->value) * slv.x[it->column];

            info(ctx,
                 "{:+d} [{} ({})] ",
                 slv.factor(it->value),
                 names[it->column],
                 slv.x[it->column]);
        }

        info(ctx, " <= {} | value: {}\n", slv.bound_max(k));
    }
}

template<typename floatingpointT>
floatingpointT
max_cost_init(const std::unique_ptr<floatingpointT[]>& c,
              int n,
              minimize_tag) noexcept
{
    auto it = std::max_element(c.get(), c.get() + n);
    if (it == c.get() + n)
        return *it;

    return std::numeric_limits<floatingpointT>::max();
}

template<typename floatingpointT>
floatingpointT
max_cost_init(const std::unique_ptr<floatingpointT[]>& c,
              int n,
              maximize_tag) noexcept
{
    auto it = std::min_element(c.get(), c.get() + n);
    if (it == c.get() + n)
        return *it;

    return std::numeric_limits<floatingpointT>::lowest();
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

    bounds_printer(floatingpointT max_cost_init)
      : bestlb(std::numeric_limits<floatingpointT>::lowest())
      , bestub(std::numeric_limits<floatingpointT>::max())
      , max_cost(max_cost_init)
    {}

    template<typename SolverT>
    floatingpointT init_bound(const SolverT& slv)
    {
        floatingpointT b{ 0 };

        for (auto c = 0; c != slv.m; ++c)
            b += slv.pi[c] * static_cast<floatingpointT>(slv.bound_init(c));

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
        floatingpointT lb = init_bound(slv);
        floatingpointT ub = init_ub(modeT());

        if (not best.solutions.empty())
            ub = static_cast<floatingpointT>(best.solutions.back().value);

        for (auto j = 0; j != slv.n; ++j)
            lb += add_bound(slv, j, slv.compute_sum_A_pi(j), modeT());

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
        compute_missing_constraint(s, R);
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

        return compute_missing_constraint(solver, R);
    }

    template<typename solverT>
    int run(solverT& solver,
            floatingpointT kappa,
            floatingpointT delta,
            floatingpointT theta)
    {
        for (auto it = R.begin(), et = R.end(); it != et; ++it)
            solver.compute_update_row(*it, kappa, delta, theta);

        return compute_missing_constraint(solver, R);
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
        compute_missing_constraint(s, R);
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

        return compute_missing_constraint(solver, R);
    }

    template<typename solverT>
    int run(solverT& solver,
            floatingpointT kappa,
            floatingpointT delta,
            floatingpointT theta)
    {
        for (auto it = R.rbegin(), et = R.rend(); it != et; ++it)
            solver.compute_update_row(*it, kappa, delta, theta);

        return compute_missing_constraint(solver, R);
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
        compute_missing_constraint(s, R);
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

        return compute_missing_constraint(solver, R);
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

        return compute_missing_constraint(solver, R);
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
            sparse_matrix<int>::const_row_iterator it, et;
            std::tie(it, et) = solver.ap.row(k);
            int v = 0;

            for (; it != et; ++it)
                v += solver.factor(it->value) * solver.x[it->column];

            if (solver.bound_min(k) > v)
                R.push_back(std::make_pair(k, solver.bound_min(k) - v));
            else if (solver.bound_max(k) < v)
                R.push_back(std::make_pair(k, v - solver.bound_max(k)));
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
        for (int k = 0; k != solver.m; ++k)
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

template<typename floatingpointT>
inline floatingpointT
compute_delta(const context_ptr& ctx,
              const std::unique_ptr<floatingpointT[]>& c,
              floatingpointT theta,
              int n)
{
    info(ctx, "  - delta not defined, compute it:\n");

    auto mini{ std::numeric_limits<floatingpointT>::max() };
    for (int i = 0; i != n; ++i)
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

inline int
element_number(const std::vector<itm::merged_constraint>& csts)
{
    std::size_t ret{ 0 };

    for (const auto& elem : csts)
        ret += elem.elements.size();

    return numeric_cast<int>(ret);
}

template<typename SolverT,
         typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT>
struct solver_functor
{
    using floatingpoint_type = floatingpointT;
    using mode_type = modeT;
    using constraint_order_type = constraintOrderT;
    using random_type = randomT;

    std::chrono::time_point<std::chrono::steady_clock> m_begin;
    std::chrono::time_point<std::chrono::steady_clock> m_end;

    const context_ptr& m_ctx;
    randomT& m_rng;

    result m_best;

    solver_functor(const context_ptr& ctx,
                   randomT& rng,

                   const std::vector<std::string>& variable_names,
                   const affected_variables& affected_vars)
      : m_ctx(ctx)
      , m_rng(rng)
    {
        m_best.affected_vars = affected_vars;
        m_best.variable_name = variable_names;
    }

    result operator()(const std::vector<itm::merged_constraint>& constraints,
                      int variables,
                      const std::unique_ptr<floatingpointT[]>& original_costs,
                      const std::unique_ptr<floatingpointT[]>& norm_costs,
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
                             ? compute_delta<floatingpoint_type>(
                                 m_ctx, norm_costs, theta, variables)
                             : static_cast<floatingpoint_type>(p.delta);

        const auto pushing_k_factor =
          static_cast<floatingpoint_type>(p.pushing_k_factor);
        const auto pushing_objective_amplifier =
          static_cast<floatingpoint_type>(p.pushing_objective_amplifier);

        auto kappa = kappa_min;

        SolverT slv(m_rng,
                    length(constraints),
                    variables,
                    norm_costs,
                    constraints,
                    p.init_policy,
                    p.init_random);

        constraint_order_type compute(m_ctx, slv, m_rng);

        auto max_cost = max_cost_init(original_costs, variables, mode_type());
        bounds_printer<floatingpointT, modeT> bound_print(max_cost);

        info(m_ctx, "* solver starts:\n");

        m_best.variables = slv.m;
        m_best.constraints = slv.n;

        bool start_pushing = false;

        for (;;) {
            int remaining = compute.run(slv, kappa, delta, theta);

            if (best_remaining == -1 or remaining < best_remaining) {
                best_remaining = remaining;
                m_best.duration = compute_duration(m_begin, m_end);
                m_best.loop = i;
                m_best.remaining_constraints = remaining;

                if (remaining == 0) {
                    m_best.status = baryonyx::result_status::success;
                    store_if_better(
                      slv.results(original_costs, cost_constant), slv.x, i);
                    start_pushing = true;
                } else {

                    // bound_print(slv, m_ctx, m_best);

                    info(m_ctx,
                         "  - violated constraints: {}/{} at {}s (loop: {})\n",
                         remaining,
                         slv.m,
                         compute_duration(m_begin, m_end),
                         i);
                }
            }

#ifndef BARYONYX_FULL_OPTIMIZATION
            print_solver(slv, m_ctx, m_best.variable_name, p.print_level);
#endif

            if (start_pushing) {
                ++pushing_iteration;

                if (pushed == -1)
                    info(m_ctx, "  - start push system:\n");

                if (pushing_iteration >= p.pushing_iteration_limit) {
                    pushed++;
                    pushing_iteration = 0;

                    remaining =
                      compute.push_and_run(slv,
                                           pushing_k_factor * kappa,
                                           delta,
                                           theta,
                                           pushing_objective_amplifier);

                    if (remaining == 0)
                        store_if_better(
                          slv.results(original_costs, cost_constant),
                          slv.x,
                          i);
                }

                if (pushed > p.pushes_limit) {
                    info(
                      m_ctx,
                      "    - Push system limit reached. Solution found: {}\n",
                      best_solution_value(m_best));
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

                // if (context_get_integer_parameter(m_ctx, "print-level", 0) >
                // 0)
                //     print_missing_constraint(m_ctx,
                //                              slv.ap,
                //                              slv.A,
                //                              m_best.variable_value,
                //                              slv.b,
                //                              m_variable_names);

                return m_best;
            }

            if (kappa > kappa_max) {
                info(m_ctx, "  - Kappa max reached: {:+.6f}\n", kappa);
                if (pushed == -1)
                    m_best.status = result_status::kappa_max_reached;

                // if (context_get_integer_parameter(m_ctx, "print-level", 0) >
                // 0)
                //     print_missing_constraint(m_ctx,
                //                              slv.ap,
                //                              slv.A,
                //                              m_best.variable_value,
                //                              slv.b,
                //                              m_variable_names);

                return m_best;
            }

            m_end = std::chrono::steady_clock::now();
            if (is_time_limit(p.time_limit, m_begin, m_end)) {
                info(m_ctx, "  - Time limit reached: {} {:+.6f}\n", i, kappa);
                if (pushed == -1)
                    m_best.status = result_status::time_limit_reached;

                // if (context_get_integer_parameter(m_ctx, "print-level", 0) >
                // 0)
                //     print_missing_constraint(m_ctx,
                //                              slv.ap,
                //                              slv.A,
                //                              m_best.variable_value,
                //                              slv.b,
                //                              m_variable_names);

                return m_best;
            }
        }
    }

private:
    void store_if_better(double current, const x_type& x, int i)
    {
        if (m_best.solutions.empty() or
            is_better_solution(
              current, m_best.solutions.back().value, modeT())) {
            m_best.solutions.emplace_back(x, current);
            m_best.duration = compute_duration(m_begin, m_end);
            m_best.loop = i;

            info(m_ctx,
                 "  - Solution found: {} (i={} t={}s)\n",
                 current,
                 m_best.loop,
                 m_best.duration);

#if 0
            std::ofstream ofs("temp.sol");
            ofs << m_best << m_best.affected_vars
                << best_solution_writer(m_best);
#endif
        } else {
            m_best.solutions.emplace_back(x, current);

            std::swap(*(m_best.solutions.rbegin()),
                      *(m_best.solutions.rbegin() + 1));
        }
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
                is_better_solution(current, m_best, modeT())) {

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

template<typename SolverT,
         typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT>
struct optimize_functor
{
    using floatingpoint_type = floatingpointT;
    using mode_type = modeT;
    using constraint_order_type = constraintOrderT;
    using random_type = randomT;

    std::chrono::time_point<std::chrono::steady_clock> m_begin;
    std::chrono::time_point<std::chrono::steady_clock> m_end;

    std::set<solution> m_all_solutions;

    const context_ptr& m_ctx;
    randomT m_rng;
    int m_thread_id;

    result m_best;

    optimize_functor(const context_ptr& ctx,
                     int thread_id,
                     typename random_type::result_type seed,
                     const std::vector<std::string>& variable_names,
                     const affected_variables& affected_vars)
      : m_ctx(ctx)
      , m_rng(seed)
      , m_thread_id(thread_id)
    {
        m_best.affected_vars = affected_vars;
        m_best.variable_name = variable_names;
    }

    result operator()(
      best_solution_recorder<floatingpointT, modeT>& best_recorder,
      const std::vector<itm::merged_constraint>& constraints,
      int variables,
      const std::unique_ptr<floatingpointT[]>& original_costs,
      const std::unique_ptr<floatingpointT[]>& norm_costs,
      double cost_constant,
      const itm::parameters& p)
    {
        m_begin = std::chrono::steady_clock::now();
        m_end = m_begin;

        int i = 0;
        int pushed = -1;
        int pushing_iteration = 0;

        itm::init_policy_type init_policy = p.init_policy;
        const auto kappa_min = static_cast<floatingpoint_type>(p.kappa_min);
        const auto kappa_step = static_cast<floatingpoint_type>(p.kappa_step);
        const auto kappa_max = static_cast<floatingpoint_type>(p.kappa_max);
        const auto alpha = static_cast<floatingpoint_type>(p.alpha);
        const auto theta = static_cast<floatingpoint_type>(p.theta);
        const auto delta = p.delta < 0
                             ? compute_delta<floatingpoint_type>(
                                 m_ctx, norm_costs, theta, variables)
                             : static_cast<floatingpoint_type>(p.delta);

        const auto pushing_k_factor =
          static_cast<floatingpoint_type>(p.pushing_k_factor);
        const auto pushing_objective_amplifier =
          static_cast<floatingpoint_type>(p.pushing_objective_amplifier);

        auto kappa = kappa_min;

        SolverT slv(m_rng,
                    length(constraints),
                    variables,
                    norm_costs,
                    constraints,
                    p.init_policy,
                    p.init_random);

        constraint_order_type compute(m_ctx, slv, m_rng);

        for (; not is_time_limit(p.time_limit, m_begin, m_end);
             m_end = std::chrono::steady_clock::now(), ++i) {

            int remaining = compute.run(slv, kappa, delta, theta);
            if (remaining == 0)
                store_if_better(slv.results(original_costs, cost_constant),
                                slv.x,
                                i,
                                best_recorder);

            if (i > p.w)
                kappa += kappa_step *
                         std::pow(static_cast<floatingpoint_type>(remaining) /
                                    static_cast<floatingpoint_type>(slv.m),
                                  alpha);

            if (i >= p.limit or kappa > kappa_max or pushed > p.pushes_limit) {
                if (m_best.solutions.empty())
                    init_policy =
                      init_solver(slv, x_type(), init_policy, p.init_random);
                else
                    init_policy =
                      init_solver(slv,
                                  m_best.solutions.back().variables,
                                  init_policy,
                                  p.init_random);

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

                    if (remaining == 0)
                        store_if_better(
                          slv.results(original_costs, cost_constant),
                          slv.x,
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
    void store_if_better(
      double current,
      const x_type& x,
      int i,
      best_solution_recorder<floatingpointT, modeT>& best_recorder)
    {
        if (m_best.solutions.empty() or
            is_better_solution(
              current, m_best.solutions.back().value, modeT())) {
            m_best.status = baryonyx::result_status::success;

            // Store only the best solution, other solutions are stored into
            // the @c std::set to avoid duplicated solutions.

            if (m_best.solutions.empty())
                m_best.solutions.emplace_back(x, current);
            else
                m_best.solutions[0] = { x, current };

            m_best.duration = compute_duration(m_begin, m_end);
            m_best.loop = i;

            best_recorder.try_update(m_best);

#if 0
                std::ofstream ofs(fmt::format("temp-{}.sol", m_thread_id));
                ofs << m_best << m_affected_vars
                    << best_solution_writer(m_best);
#endif
        }

        m_all_solutions.emplace(x, current);
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

template<typename floatingpointT, typename randomT>
inline std::unique_ptr<floatingpointT[]>
rng_normalize_costs(const std::unique_ptr<floatingpointT[]>& c,
                    randomT& rng,
                    int n)
{
    std::vector<std::pair<floatingpointT, int>> r(n);

    for (int i = 0; i != n; ++i) {
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

    // Reorder the vector according to the variable index, so, it restores
    // the initial order.

    std::sort(r.begin(), r.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.second < rhs.second;
    });

    auto ret = std::make_unique<floatingpointT[]>(n);
    for (int i = 0; i != n; ++i)
        ret[i] = r[i].first;

    // Finally we compute the l+oo norm.

    floatingpointT div = *std::max_element(c.get(), c.get() + n);
    if (std::isnormal(div))
        for (int i = 0; i != n; ++i)
            ret[i] /= div;

    return ret;
}

/**
 * Normalizes the cost vector, i.e. divides it by its l{1,2, +oo}norm. If
 * the input vector is too small or with infinity value, the c is
 * unchanged.
 */
template<typename floatingpointT, typename randomT>
inline std::unique_ptr<floatingpointT[]>
normalize_costs(const context_ptr& ctx,
                const std::string& norm,
                const std::unique_ptr<floatingpointT[]>& c,
                randomT& rng,
                int n)
{
    auto ret = std::make_unique<floatingpointT[]>(n);
    std::copy(c.get(), c.get() + n, ret.get());

    if (norm == "none") {
        info(ctx, "  - No norm");
        return ret;
    }

    if (norm == "rng") {
        info(ctx, "  - Compute random norm\n");
        return rng_normalize_costs<floatingpointT, randomT>(c, rng, n);
    }

    floatingpointT div{ 0 };

    if (norm == "l1") {
        info(ctx, "  - Compute l1 norm\n");
        for (int i = 0; i != n; ++i)
            div += std::abs(ret[i]);
    } else if (norm == "l2") {
        info(ctx, "  - Compute l2 norm\n");
        for (int i = 0; i != n; ++i)
            div += ret[i] * ret[i];
    } else {
        info(ctx, "  - Compute infinity-norm (default)\n");
        div = *std::max_element(c.get(), c.get() + n);
    }

    if (std::isnormal(div))
        for (int i = 0; i != n; ++i)
            ret[i] /= div;

    return ret;
}

template<typename floatingpointT>
inline std::unique_ptr<floatingpointT[]>
make_objective_function(const objective_function& obj, int n)
{
    auto ret = std::make_unique<floatingpointT[]>(n);

    for (const auto& elem : obj.elements) {
        assert(0 <= n and elem.variable_index < n);
        ret[elem.variable_index] += static_cast<floatingpointT>(elem.factor);
    }

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

        auto variables = numeric_cast<int>(pb.vars.values.size());
        auto cost =
          make_objective_function<floatingpointT>(pb.objective, variables);
        auto norm_costs = normalize_costs<floatingpointT, randomT>(
          ctx, p.norm, cost, rng, variables);
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
inline std::unique_ptr<typename randomT::result_type[]>
generate_seed(randomT& rng, int thread)
{
    using type = typename randomT::result_type;

    auto ret = std::make_unique<type[]>(thread);

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

        auto variables = numeric_cast<int>(pb.vars.values.size());
        auto cost =
          make_objective_function<floatingpointT>(pb.objective, variables);
        auto norm_costs = normalize_costs<floatingpointT, randomT>(
          ctx, p.norm, cost, rng, variables);
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

        std::set<solution> all_solutions;

        for (int i = 0; i != thread; ++i) {
            auto current = results[i].get();
            if (current.status == result_status::success) {
                all_solutions.insert(current.solutions.begin(),
                                     current.solutions.end());

                if (ret.solutions.empty() or
                    is_better_solution(current.solutions.back().value,
                                       ret.solutions.back().value,
                                       modeT()))
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
}
}

#endif
