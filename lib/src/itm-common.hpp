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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_COMMON_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_COMMON_HPP

#include <baryonyx/core>

#include "private.hpp"
#include "utils.hpp"
#include "sparse-matrix.hpp"

#include <algorithm>
#include <chrono>
#include <iterator>
#include <random>
#include <tuple>
#include <vector>

#include <cassert>

namespace baryonyx {
namespace itm {

/**
 * x_type is a std::vector<bool> instead of a baryonyx::fixed_array<bool> to
 * use the specialized version of vector, which is used for elements of type
 * bool and optimizes for space.
 */
using x_type = std::vector<bool>;

struct maximize_tag
{};

struct minimize_tag
{};

inline double
best_solution_value(const baryonyx::result& res) noexcept
{
    assert(res.status == baryonyx::result_status::success);
    assert(not res.solutions.empty());

    return res.solutions.back().value;
}

inline bool
is_better_solution(const baryonyx::result& lhs,
                   const baryonyx::result& rhs,
                   maximize_tag) noexcept
{
    return best_solution_value(lhs) > best_solution_value(rhs);
}

inline bool
is_better_solution(const baryonyx::result& lhs,
                   const baryonyx::result& rhs,
                   minimize_tag) noexcept
{
    return best_solution_value(lhs) < best_solution_value(rhs);
}

struct merged_constraint
{
    merged_constraint(std::vector<function_element> elements_,
                      int min_,
                      int max_,
                      int id_)
      : elements(std::move(elements_))
      , min(min_)
      , max(max_)
      , id(id_)
    {}

    std::vector<function_element> elements;
    int min;
    int max;
    int id;
};

std::vector<merged_constraint>
make_merged_constraints(const context_ptr& ctx, const problem& pb);

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
  const std::vector<merged_constraint>& csts) noexcept
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
inline solver_parameters::init_policy_type
init_solver(Solver& slv,
            const x_type& best_previous,
            solver_parameters::init_policy_type type,
            double init_random)
{
    using floatingpointT = typename Solver::floatingpoint_type;

    std::fill(slv.P.get(),
              slv.P.get() + slv.ap.length(),
              static_cast<floatingpointT>(0));
    std::fill(
      slv.pi.get(), slv.pi.get() + slv.m, static_cast<floatingpointT>(0));

    if (best_previous.empty() and
        type == solver_parameters::init_policy_type::best)
        type = solver_parameters::init_policy_type::bastert;

    init_random = clamp(init_random, 0.0, 1.0);
    std::bernoulli_distribution d(init_random);

    switch (type) {
    case solver_parameters::init_policy_type::bastert:
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

        type = solver_parameters::init_policy_type::random;
        break;
    case solver_parameters::init_policy_type::random:
        for (int i = 0; i != slv.n; ++i)
            slv.x[i] = d(slv.rng);

        type = solver_parameters::init_policy_type::best;
        break;
    case solver_parameters::init_policy_type::best:
        for (int i = 0; i != slv.n; ++i)
            slv.x[i] = (d(slv.rng)) ? (best_previous[i]) : (!best_previous[i]);

        type = solver_parameters::init_policy_type::bastert;
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

namespace details {
inline int
constraint(std::vector<int>::iterator it)
{
    return *it;
}

inline int
constraint(std::vector<int>::reverse_iterator it)
{
    return *it;
}

inline int
constraint(std::vector<std::pair<int, int>>::iterator it)
{
    return it->first;
}

} // namespace details

template<typename Iterator>
int
constraint(Iterator it)
{
    return details::constraint(it);
}

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
        solver.push_and_compute_update_row(
          R.begin(), R.end(), kappa, delta, theta, objective_amplifier);

        return compute_missing_constraint(solver, R);
    }

    template<typename solverT>
    int run(solverT& solver,
            floatingpointT kappa,
            floatingpointT delta,
            floatingpointT theta)
    {
        solver.compute_update_row(R.begin(), R.end(), kappa, delta, theta);

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
        solver.push_and_compute_update_row(
          R.rbegin(), R.rend(), kappa, delta, theta, objective_amplifier);

        return compute_missing_constraint(solver, R);
    }

    template<typename solverT>
    int run(solverT& solver,
            floatingpointT kappa,
            floatingpointT delta,
            floatingpointT theta)
    {
        solver.compute_update_row(R.rbegin(), R.rend(), kappa, delta, theta);

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
        std::shuffle(R.begin(), R.end(), rng);

        solver.push_and_compute_update_row(
          R.begin(), R.end(), kappa, delta, theta, objective_amplifier);

        return compute_missing_constraint(solver, R);
    }

    template<typename solverT>
    int run(solverT& solver,
            floatingpointT kappa,
            floatingpointT delta,
            floatingpointT theta)
    {
        std::shuffle(R.begin(), R.end(), rng);

        solver.compute_update_row(R.begin(), R.end(), kappa, delta, theta);

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
        baryonyx::itm::sort(R.begin(), R.end(), direction_type());

        solver.push_and_compute_update_row(
          R.begin(), R.end(), kappa, delta, theta, objective_amplifier);

        return local_compute_missing_constraint(solver);
    }

    template<typename solverT>
    int run(solverT& solver,
            floatingpointT kappa,
            floatingpointT delta,
            floatingpointT theta)
    {
        baryonyx::itm::sort(R.begin(), R.end(), direction_type());

        solver.compute_update_row(R.begin(), R.end(), kappa, delta, theta);

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
element_number(const std::vector<merged_constraint>& csts)
{
    std::size_t ret{ 0 };

    for (const auto& elem : csts)
        ret += elem.elements.size();

    return numeric_cast<int>(ret);
}

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
                const std::unique_ptr<floatingpointT[]>& c,
                randomT& rng,
                int n)
{
    auto ret = std::make_unique<floatingpointT[]>(n);
    std::copy(c.get(), c.get() + n, ret.get());

    if (ctx->parameters.cost_norm == solver_parameters::cost_norm_type::none) {
        info(ctx, "  - No norm");
        return ret;
    }

    if (ctx->parameters.cost_norm ==
        solver_parameters::cost_norm_type::random) {
        info(ctx, "  - Compute random norm\n");
        return rng_normalize_costs<floatingpointT, randomT>(c, rng, n);
    }

    floatingpointT div{ 0 };

    if (ctx->parameters.cost_norm == solver_parameters::cost_norm_type::l1) {
        info(ctx, "  - Compute l1 norm\n");
        for (int i = 0; i != n; ++i)
            div += std::abs(ret[i]);
    } else if (ctx->parameters.cost_norm ==
               solver_parameters::cost_norm_type::l2) {
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
    auto param = ctx->parameters.seed;

    if (param <= 0)
        return static_cast<typename randomT::result_type>(epoch);

    return static_cast<typename randomT::result_type>(param);
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

} // namespace itm
} // namespace baryonyx

#endif
