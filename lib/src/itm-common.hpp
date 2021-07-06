/* Copyright (C) 2016-2021 INRAE
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

#include "bit-array.hpp"
#include "debug.hpp"
#include "private.hpp"
#include "problem.hpp"
#include "sparse-matrix.hpp"
#include "utils.hpp"

#include <algorithm>
#include <chrono>
#include <numeric>
#include <random>
#include <thread>
#include <vector>

namespace baryonyx {

// We statically define this PRNG for all subsystems in Baryonyx.
using random_engine = std::default_random_engine;

constexpr real Zero{ 0 };
constexpr real One{ 1 };
constexpr real Two{ 2 };
constexpr real Middle{ (One + Two) / Two };

namespace itm {

struct maximize_tag
{};

struct minimize_tag
{};

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
make_merged_constraints(const context& ctx, const problem& pb);

template<typename Solver, typename Xtype>
bool
is_valid_constraint(const Solver& slv, int k, const Xtype& x)
{
    typename sparse_matrix<int>::const_row_iterator it, et;

    std::tie(it, et) = slv.ap.row(k);
    int v = 0;

    for (; it != et; ++it)
        v += slv.factor(it->value) * x[it->column];

    return slv.bound_min(k) <= v && v <= slv.bound_max(k);
}

template<typename Solver, typename Xtype>
bool
is_valid_solution(const Solver& slv, const Xtype& x)
{
    for (int k = 0; k != slv.m; ++k)
        if (!is_valid_constraint(slv, k, x))
            return false;

    return true;
}

template<typename Solver, typename Xtype>
int
compute_violated_constraints(const Solver& slv,
                             const Xtype& x,
                             std::vector<int>& out)
{
    out.clear();

    for (int k = 0; k != slv.m; ++k)
        if (!is_valid_constraint(slv, k, x))
            out.emplace_back(k);

    return length(out);
}

template<typename iteratorT>
inline void
random_shuffle_unique(iteratorT begin,
                      iteratorT end,
                      random_engine& rng) noexcept
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

template<typename Mode, typename Iterator>
inline void
calculator_sort(Iterator begin, Iterator end, random_engine& rng)
{
    if (std::distance(begin, end) > 1) {
        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            if constexpr (std::is_same_v<Mode, minimize_tag>)
                return lhs.value < rhs.value;
            else
                return rhs.value < lhs.value;
        });

        random_shuffle_unique(begin, end, rng);
    }
}

template<typename Mode>
inline bool
stop_iterating(real value, random_engine& rng) noexcept
{
    if (value == 0) {
        std::bernoulli_distribution d(0.5);
        return d(rng);
    }

    if constexpr (std::is_same_v<Mode, minimize_tag>)
        return value > 0;
    else
        return value < 0;
}

template<typename Mode>
constexpr real
bad_value() noexcept
{
    if constexpr (std::is_same_v<Mode, minimize_tag>)
        return std::numeric_limits<real>::infinity();
    else
        return -std::numeric_limits<real>::infinity();
}

template<typename Mode>
constexpr bool
is_bad_value(real value) noexcept
{
    return value == bad_value<Mode, real>();
}

template<typename Mode>
inline bool
stop_iterating(real value) noexcept
{
    if constexpr (std::is_same_v<Mode, minimize_tag>)
        return value > 0;
    else
        return value < 0;
}

template<typename Mode>
inline bool
is_better_solution(real lhs, real rhs) noexcept
{
    if constexpr (std::is_same_v<Mode, minimize_tag>)
        return lhs < rhs;
    else
        return lhs > rhs;
}

template<typename Mode>
inline bool
init_x(real cost, int value_if_cost_0) noexcept
{
    if constexpr (std::is_same_v<Mode, minimize_tag>) {
        if (cost < 0)
            return true;

        if (cost == 0)
            return value_if_cost_0;
    } else {
        if (cost > 0)
            return true;

        if (cost == 0)
            return value_if_cost_0;
    }

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

/// Test if the bit sign between two reals are differents.
inline bool
is_signbit_change(real lhs, real rhs) noexcept
{
    return std::signbit(lhs) != std::signbit(rhs);
}

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
                    const std::vector<merged_constraint>& constraints,
                    const double init_random) noexcept
{
    int max_length = 0;
    for (const auto& cst : constraints)
        if (length(cst.elements) > max_length)
            max_length = length(cst.elements);

    struct reduced_cost
    {
        real value;
        int factor;
        int id;
    };

    std::vector<reduced_cost> R(max_length);
    std::bernoulli_distribution dist(init_random);

    for (const auto& cst : constraints) {
        if (!dist(rng))
            continue;

        R.resize(cst.elements.size());
        const int r_size = length(cst.elements);

        for (int i = 0; i != r_size; ++i) {
            R[i].value = static_cast<real>(c[cst.elements[i].variable_index]);
            R[i].factor = cst.elements[i].factor;
            R[i].id = cst.elements[i].variable_index;
        }

        std::shuffle(std::begin(R), std::end(R), rng);
        std::sort(
          std::begin(R), std::end(R), [](const auto& lhs, const auto& rhs) {
              return is_better_solution<Mode>(lhs.value, rhs.value);
          });

        if (!x_pessimistic.empty()) {
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

        if (!x_optimistic.empty()) {
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

/// The bkmin and bkmax constraint bounds are not equal and can be assigned to
/// -infinity or +infinity. We have to scan the r vector and search a value j
/// such as b(0, k) <= Sum A(k, R[j]) < b(1, k).
///
/// @return True if the sign of the pi vector change during the affect
///     operation.
template<typename Solver, typename Xtype, typename Iterator>
bool
affect(Solver& slv,
       Xtype& x,
       Iterator it,
       int k,
       int selected,
       int r_size,
       const real kappa,
       const real delta)
{
    constexpr real one{ 1 };
    constexpr real two{ 2 };
    constexpr real middle{ (two + one) / two };

    const auto old_pi = slv.pi[k];

    auto d = delta;

    if (selected < 0) {
        // slv.pi[k] += slv.R[0].value / two;
        d += (kappa / (one - kappa)) * (slv.R[0].value / two);

        for (int i = 0; i != r_size; ++i) {
            auto var = it + slv.R[i].id;

            if (slv.R[i].is_negative_factor()) {
                x.set(var->column);
                slv.P[var->value] += d;
            } else {
                x.unset(var->column);
                slv.P[var->value] -= d;
            }
        }
    } else if (selected + 1 >= r_size) {
        // slv.pi[k] += slv.R[selected].value * middle;
        d += (kappa / (one - kappa)) * (slv.R[selected].value * middle);

        for (int i = 0; i != r_size; ++i) {
            auto var = it + slv.R[i].id;

            if (slv.R[i].is_negative_factor()) {
                x.unset(var->column);
                slv.P[var->value] -= d;
            } else {
                x.set(var->column);
                slv.P[var->value] += d;
            }
        }
    } else {
        slv.pi[k] +=
          ((slv.R[selected].value + slv.R[selected + 1].value) / two);
        d += (kappa / (one - kappa)) *
             (slv.R[selected + 1].value - slv.R[selected].value);

        int i = 0;
        for (; i <= selected; ++i) {
            auto var = it + slv.R[i].id;

            if (slv.R[i].is_negative_factor()) {
                x.unset(var->column);
                slv.P[var->value] -= d;
            } else {
                x.set(var->column);
                slv.P[var->value] += d;
            }
        }

        for (; i != r_size; ++i) {
            auto var = it + slv.R[i].id;

            if (slv.R[i].is_negative_factor()) {
                x.set(var->column);
                slv.P[var->value] += d;
            } else {
                x.unset(var->column);
                slv.P[var->value] -= d;
            }
        }
    }

    // TODO job: develops is_valid_constraint for all the solvers
    bx_expects(is_valid_constraint(slv, k, x));

    return is_signbit_change(old_pi, slv.pi[k]);
}

namespace detail {

inline int
constraint(std::vector<int>::const_iterator it)
{
    return *it;
}

inline int
constraint(std::vector<int>::const_reverse_iterator it)
{
    return *it;
}

inline int
constraint(std::vector<std::pair<int, int>>::const_iterator it)
{
    return it->first;
}
}

template<typename Iterator>
int
constraint(Iterator it)
{
    return detail::constraint(it);
}

/**
 * Compute a problem lower or upper bounds based on Lagrangian multipliers
 * (valid if there are equality constraints only?)
 */
template<typename real, typename modeT>
struct bounds_printer
{
    real bestlb;
    real bestub;
    real max_cost;

    bounds_printer(real max_cost_init)
      : bestlb(std::numeric_limits<real>::lowest())
      , bestub(std::numeric_limits<real>::max())
      , max_cost(max_cost_init)
    {}

    template<typename SolverT>
    real init_bound(const SolverT& slv)
    {
        real b{ 0 };

        for (auto c = 0; c != slv.m; ++c)
            b += slv.pi[c] * static_cast<real>(slv.bound_init(c));

        return b;
    }

    template<typename SolverT>
    real add_bound(const SolverT& slv, int j, real sum_a_pi, minimize_tag)
    {
        if (slv.c[j] - sum_a_pi < 0)
            return slv.c[j] - sum_a_pi;

        return { 0 };
    }

    template<typename SolverT>
    real add_bound(const SolverT& slv, int j, real sum_a_pi, maximize_tag)
    {
        if (slv.c[j] - sum_a_pi > 0)
            return slv.c[j] - sum_a_pi;

        return { 0 };
    }

    void print_bound(const context& ctx,
                     real lower_bound,
                     real upper_bound,
                     minimize_tag)
    {
        bool better_gap = (lower_bound > bestlb || upper_bound < bestub);

        if (upper_bound < bestub)
            bestub = upper_bound;

        if (lower_bound > bestlb)
            bestlb = lower_bound;

        if (better_gap) {
            if (bestub == static_cast<real>(0.0))
                info(ctx, "  - Lower bound: {}   (gap: 0%)\n", bestlb);
            else
                info(ctx,
                     "  - Lower bound: {}   (gap: {}%)\n",
                     bestlb,
                     static_cast<real>(100.) * (bestub - bestlb) / bestub);
        }
    }

    void print_bound(const context& ctx,
                     real lower_bound,
                     real upper_bound,
                     maximize_tag)
    {
        bool better_gap = (lower_bound > bestlb || upper_bound < bestub);

        if (upper_bound < bestub)
            bestub = upper_bound;

        if (lower_bound > bestlb)
            bestlb = lower_bound;

        if (better_gap) {
            if (bestlb == static_cast<real>(0.0))
                info(ctx, "  - Upper bound: {}   (gap: 0%)\n", bestub);
            else
                info(ctx,
                     "  - Upper bound: {}   (gap: {}%)\n",
                     bestub,
                     static_cast<real>(100.) * (bestlb - bestub) / bestlb);
        }
    }

    real init_ub(minimize_tag)
    {
        return std::numeric_limits<real>::max();
    }

    real init_ub(maximize_tag)
    {
        return std::numeric_limits<real>::lowest();
    }

    template<typename SolverT>
    void operator()(const SolverT& slv, const context& ctx, const result& best)
    {
        real lb = init_bound(slv);
        real ub = init_ub(modeT());

        if (not best.solutions.empty())
            ub = static_cast<real>(best.solutions.back().value);

        for (auto j = 0; j != slv.n; ++j)
            lb += add_bound(slv, j, slv.compute_sum_A_pi(j), modeT());

        lb *= max_cost; // restore original cost

        print_bound(ctx, lb, ub, modeT());
    }
};

struct compute_order
{
    std::vector<int> R;
    std::vector<std::pair<int, int>> m_order;
    solver_parameters::constraint_order order;
    bool use_cycle;

    compute_order(const solver_parameters::constraint_order order_,
                  const int variable_number)
      : R(static_cast<std::vector<int>::size_type>(variable_number))
      , m_order(static_cast<std::vector<int>::size_type>(variable_number))
      , order(order_ == solver_parameters::constraint_order::cycle
                ? solver_parameters::constraint_order::none
                : order_)
      , use_cycle(order_ == solver_parameters::constraint_order::cycle)
    {}

    static solver_parameters::constraint_order next_state(
      solver_parameters::constraint_order current) noexcept
    {
        using int_type = typename std::underlying_type<
          solver_parameters::constraint_order>::type;

        auto int_current = static_cast<int_type>(current);
        auto int_max =
          static_cast<int_type>(solver_parameters::constraint_order::cycle);

        ++int_current;

        return static_cast<solver_parameters::constraint_order>(
          int_current >= int_max ? 0 : int_current);
    }

    template<typename Solver, typename Xtype>
    void init(const Solver& s, const Xtype& x)
    {
        switch (order) {
        case solver_parameters::constraint_order::infeasibility_decr:
        case solver_parameters::constraint_order::infeasibility_incr:
            infeasibility_local_compute_violated_constraints(s, x);
            break;
        case solver_parameters::constraint_order::pi_sign_change:
            std::iota(R.begin(), R.end(), 0);
            break;
        case solver_parameters::constraint_order::none:
        case solver_parameters::constraint_order::reversing:
        case solver_parameters::constraint_order::random_sorting:
        case solver_parameters::constraint_order::lagrangian_decr:
        case solver_parameters::constraint_order::lagrangian_incr:
        default:
            compute_violated_constraints(s, x, R);
            break;
        }
    };

    template<typename Solver, typename Xtype>
    int push_and_run(Solver& solver,
                     Xtype& x,
                     random_engine& rng,
                     real kappa,
                     real delta,
                     real theta,
                     real objective_amplifier)
    {
        bool pi_changed = 0;
        int remaining = 0;

        if (use_cycle)
            order = next_state(order);

        switch (order) {
        case solver_parameters::constraint_order::reversing:
            solver.push_and_compute_update_row(x,
                                               R.crbegin(),
                                               R.crend(),
                                               kappa,
                                               delta,
                                               theta,
                                               objective_amplifier);

            return compute_violated_constraints(solver, x, R);
        case solver_parameters::constraint_order::random_sorting:
            std::shuffle(R.begin(), R.end(), rng);

            solver.push_and_compute_update_row(
              x, R.begin(), R.end(), kappa, delta, theta, objective_amplifier);

            return compute_violated_constraints(solver, x, R);
        case solver_parameters::constraint_order::infeasibility_decr:
            std::sort(m_order.begin(),
                      m_order.end(),
                      [](const auto& lhs, const auto& rhs) {
                          return rhs.second < lhs.second;
                      });

            solver.push_and_compute_update_row(x,
                                               m_order.cbegin(),
                                               m_order.cend(),
                                               kappa,
                                               delta,
                                               theta,
                                               objective_amplifier);

            return infeasibility_local_compute_violated_constraints(solver, x);
        case solver_parameters::constraint_order::infeasibility_incr:
            std::sort(m_order.begin(),
                      m_order.end(),
                      [](const auto& lhs, const auto& rhs) {
                          return lhs.second < rhs.second;
                      });
            solver.push_and_compute_update_row(x,
                                               m_order.cbegin(),
                                               m_order.cend(),
                                               kappa,
                                               delta,
                                               theta,
                                               objective_amplifier);

            return infeasibility_local_compute_violated_constraints(solver, x);
        case solver_parameters::constraint_order::lagrangian_decr:
            std::sort(
              R.begin(), R.end(), [&solver](const auto lhs, const auto rhs) {
                  return solver.pi[rhs] < solver.pi[lhs];
              });

            solver.push_and_compute_update_row(x,
                                               R.cbegin(),
                                               R.cend(),
                                               kappa,
                                               delta,
                                               theta,
                                               objective_amplifier);

            return compute_violated_constraints(solver, x, R);
        case solver_parameters::constraint_order::lagrangian_incr:
            std::sort(
              R.begin(), R.end(), [&solver](const auto lhs, const auto rhs) {
                  return solver.pi[lhs] < solver.pi[rhs];
              });

            solver.push_and_compute_update_row(x,
                                               R.cbegin(),
                                               R.cend(),
                                               kappa,
                                               delta,
                                               theta,
                                               objective_amplifier);

            return compute_violated_constraints(solver, x, R);
        case solver_parameters::constraint_order::pi_sign_change:
            std::shuffle(R.begin(), R.end(), rng);

            pi_changed = solver.push_and_compute_update_row(
              x, R.begin(), R.end(), kappa, delta, theta, objective_amplifier);
            remaining = local_compute_violated_constraints(solver, x);

            if (!pi_changed && remaining == 0)
                return 0;

            return remaining;
        case solver_parameters::constraint_order::none:
        default:
            solver.push_and_compute_update_row(x,
                                               R.cbegin(),
                                               R.cend(),
                                               kappa,
                                               delta,
                                               theta,
                                               objective_amplifier);
            return compute_violated_constraints(solver, x, R);
        }
    }

    template<typename Solver, typename Xtype>
    int run(Solver& solver,
            Xtype& x,
            random_engine& rng,
            real kappa,
            real delta,
            real theta)
    {
        bool pi_changed = false;
        int remaining = 0;

        switch (order) {
        case solver_parameters::constraint_order::reversing:
            solver.compute_update_row(
              x, R.crbegin(), R.crend(), kappa, delta, theta);

            return compute_violated_constraints(solver, x, R);
        case solver_parameters::constraint_order::random_sorting:
            std::shuffle(R.begin(), R.end(), rng);

            solver.compute_update_row(
              x, R.begin(), R.end(), kappa, delta, theta);

            return compute_violated_constraints(solver, x, R);
        case solver_parameters::constraint_order::infeasibility_decr:
            std::sort(m_order.begin(),
                      m_order.end(),
                      [](const auto& lhs, const auto& rhs) {
                          return rhs.second < lhs.second;
                      });
            solver.compute_update_row(
              x, m_order.cbegin(), m_order.cend(), kappa, delta, theta);

            return infeasibility_local_compute_violated_constraints(solver, x);
        case solver_parameters::constraint_order::infeasibility_incr:
            std::sort(m_order.begin(),
                      m_order.end(),
                      [](const auto& lhs, const auto& rhs) {
                          return lhs.second < rhs.second;
                      });
            solver.compute_update_row(
              x, m_order.cbegin(), m_order.cend(), kappa, delta, theta);

            return infeasibility_local_compute_violated_constraints(solver, x);
        case solver_parameters::constraint_order::lagrangian_decr:
            std::sort(
              R.begin(), R.end(), [&solver](const auto lhs, const auto rhs) {
                  return solver.pi[rhs] < solver.pi[lhs];
              });
            solver.compute_update_row(
              x, R.cbegin(), R.cend(), kappa, delta, theta);

            return compute_violated_constraints(solver, x, R);
        case solver_parameters::constraint_order::lagrangian_incr:
            std::sort(
              R.begin(), R.end(), [&solver](const auto lhs, const auto rhs) {
                  return solver.pi[lhs] < solver.pi[rhs];
              });

            solver.compute_update_row(
              x, R.cbegin(), R.cend(), kappa, delta, theta);

            return compute_violated_constraints(solver, x, R);
        case solver_parameters::constraint_order::pi_sign_change:
            std::shuffle(R.begin(), R.end(), rng);

            pi_changed = solver.compute_update_row(
              x, R.begin(), R.end(), kappa, delta, theta);
            remaining = local_compute_violated_constraints(solver, x);

            if (!pi_changed && remaining == 0)
                return 0;

            return remaining;
        case solver_parameters::constraint_order::none:
        default:
            solver.compute_update_row(
              x, R.begin(), R.end(), kappa, delta, theta);
            return compute_violated_constraints(solver, x, R);
        }
    }

    template<typename Solver, typename Xtype>
    int local_compute_violated_constraints(const Solver& slv, const Xtype& x)
    {
        int remaining = 0;
        for (int k = 0; k != slv.m; ++k)
            if (!is_valid_constraint(slv, k, x))
                ++remaining;

        return remaining;
    }

    template<typename Solver, typename Xtype>
    int infeasibility_local_compute_violated_constraints(Solver& solver,
                                                         const Xtype& x)
    {
        m_order.clear();

        for (int k = 0, e = solver.m; k != e; ++k) {
            sparse_matrix<int>::const_row_iterator it, et;
            std::tie(it, et) = solver.ap.row(k);
            int v = 0;

            for (; it != et; ++it)
                v += solver.factor(it->value) * x[it->column];

            if (solver.bound_min(k) > v)
                m_order.push_back(std::make_pair(k, solver.bound_min(k) - v));
            else if (solver.bound_max(k) < v)
                m_order.push_back(std::make_pair(k, v - solver.bound_max(k)));
        }

        return length(m_order);
    }
};

template<typename Cost>
inline real
compute_delta(const context& ctx,
              const Cost& c,
              real kappa_step,
              real theta,
              int n)
{
    info(ctx, "  - delta not defined, compute it:\n");

    const auto mini = c.min(n);
    const auto temp = mini - theta * mini;
    const auto ret = temp <= 0 ? kappa_step : temp;

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

template<typename iteratorT>
inline void
random_epsilon_unique(iteratorT begin,
                      iteratorT end,
                      random_engine& rng,
                      real min,
                      real max)
{
    bx_expects(min != max);

    std::uniform_real_distribution<real> distribution(min, max);

    for (; begin != end; ++begin)
        begin->first = distribution(rng);
}

/**
 * Normalizes the cost vector, i.e. divides it by its l{1,2, +oo}norm.
 * If the input vector is too small or with infinity value, the c is
 * unchanged.
 */
template<typename Cost>
inline Cost
normalize_costs(const context& ctx, const Cost& c, random_engine& rng, int n)
{
    Cost ret(c, n);

    switch (ctx.parameters.cost_norm) {
    case solver_parameters::cost_norm_type::none:
        info(ctx, "  - No norm");
        return ret;

    case solver_parameters::cost_norm_type::random:
        info(ctx, "  - Compute random norm\n");
        ret.make_random_norm(n, rng);
        return ret;

    case solver_parameters::cost_norm_type::l1:
        info(ctx, "  - Compute l1 norm\n");
        ret.make_l1_norm(n);
        return ret;

    case solver_parameters::cost_norm_type::l2:
        info(ctx, "  - Compute l2 norm\n");
        ret.make_l2_norm(n);
        return ret;

    default:
        info(ctx, "  - Compute infinity-norm (default)\n");
        ret.make_loo_norm(n);
        return ret;
    }
}

struct default_cost_type
{
    const baryonyx::objective_function& obj;
    std::unique_ptr<real[]> linear_elements;

    default_cost_type(const objective_function& obj_, int n)
      : obj(obj_)
      , linear_elements(std::make_unique<real[]>(n))
    {
        for (const auto& elem : obj.elements) {
            bx_ensures(0 <= elem.variable_index && elem.variable_index < n);

            linear_elements[elem.variable_index] +=
              static_cast<real>(elem.factor);
        }
    }

    default_cost_type(const default_cost_type& other, int n)
      : obj(other.obj)
      , linear_elements(std::make_unique<real[]>(n))
    {
        std::copy_n(other.linear_elements.get(), n, linear_elements.get());
    }

    template<typename Random>
    void make_random_norm(int n, Random& rng)
    {
        std::vector<std::pair<real, int>> r(n);

        for (int i = 0; i != n; ++i)
            r[i] = { linear_elements[i], i };

        std::sort(r.begin(), r.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.first < rhs.first;
        });

        auto begin = r.begin();
        auto end = r.end();
        auto next = r.begin()++;
        for (; next != end; ++next) {
            if (next->first != begin->first) {
                if (std::distance(begin, next) > 1)
                    random_epsilon_unique(
                      begin, next, rng, begin->first, next->first);

                begin = next;
            }
        }

        if (std::distance(begin, end) > 1) {
            if (begin == r.begin()) {
                random_epsilon_unique(
                  begin, end, rng, begin->first, begin->first + 1);
            } else {
                auto value = linearize(std::ptrdiff_t(0),
                                       r.front().first,
                                       std::distance(r.begin(), begin),
                                       begin->first,
                                       std::distance(r.begin(), r.end()));

                random_epsilon_unique(begin, end, rng, begin->first, value);
            }
        }

        // Reorder the vector according to the variable index, so, it
        // restores the initial order.

        std::sort(r.begin(), r.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.second < rhs.second;
        });

        for (int i = 0; i != n; ++i)
            linear_elements[i] = r[i].first;

        // Finally we compute the l+oo norm.

        real div =
          *std::max_element(linear_elements.get(), linear_elements.get() + n);
        if (std::isnormal(div))
            for (int i = 0; i != n; ++i)
                linear_elements[i] /= div;
    }

    real min(int n) const noexcept
    {
        real min = std::numeric_limits<real>::max();

        for (int i = 0; i != n; ++i)
            min = std::min(min, std::abs(linear_elements[i]));

        return min == std::numeric_limits<real>::max() ? Zero : min;
    }

    void make_l1_norm(int n)
    {
        real div = { 0 };

        for (int i = 0; i != n; ++i)
            div += std::abs(linear_elements[i]);

        if (std::isnormal(div))
            for (int i = 0; i != n; ++i)
                linear_elements[i] /= div;
    }

    void make_l2_norm(int n)
    {
        real div = { 0 };

        for (int i = 0; i != n; ++i)
            div += linear_elements[i] * linear_elements[i];

        if (std::isnormal(div))
            for (int i = 0; i != n; ++i)
                linear_elements[i] /= div;
    }

    void make_loo_norm(int n)
    {
        real div =
          *std::max_element(linear_elements.get(), linear_elements.get() + n);

        if (std::isnormal(div))
            for (int i = 0; i != n; ++i)
                linear_elements[i] /= div;
    }

    real operator[](int index) const noexcept
    {
        return linear_elements[index];
    }

    real operator()(int index,
                    [[maybe_unused]] const bit_array& x) const noexcept
    {
        return linear_elements[index];
    }

    real results(const bit_array& x, real cost_constant) const noexcept
    {
        for (int i = 0, e = length(obj.elements); i != e; ++i)
            if (x[obj.elements[i].variable_index])
                cost_constant += static_cast<real>(obj.elements[i].factor);

        return cost_constant;
    }
};

struct quadratic_cost_type
{
    struct quad
    {
        real factor;
        int id;
    };

    const baryonyx::objective_function& obj;
    std::unique_ptr<real[]> linear_elements;
    std::unique_ptr<quad[]> quadratic_elements;
    std::unique_ptr<int[]> indices;

    quadratic_cost_type(const objective_function& obj_, int n)
      : obj(obj_)
      , linear_elements(std::make_unique<real[]>(n))
      , quadratic_elements(
          std::make_unique<quad[]>(2 * obj.qelements.size() + 1))
      , indices(std::make_unique<int[]>(n + 1))
    {
        for (int i = 0, e = length(obj.elements); i != e; ++i) {
            bx_ensures(0 <= obj.elements[i].variable_index &&
                       obj.elements[i].variable_index < n);

            linear_elements[obj.elements[i].variable_index] +=
              static_cast<real>(obj.elements[i].factor);
        }

        indices[0] = 0;
        for (int var = 0; var != n; ++var) {
            indices[var + 1] = indices[var];

            for (int i = 0, e = length(obj.qelements); i != e; ++i) {
                bx_ensures(0 <= obj.qelements[i].variable_index_a &&
                           obj.qelements[i].variable_index_a < n);
                bx_ensures(0 <= obj.qelements[i].variable_index_b &&
                           obj.qelements[i].variable_index_b < n);

                if (var == obj.qelements[i].variable_index_a ||
                    var == obj.qelements[i].variable_index_b)
                    ++indices[var + 1];
            }
        }

        for (int id = 0, var = 0; var != n; ++var) {
            for (int i = 0, e = length(obj.qelements); i != e; ++i) {
                const bool is_a = obj.qelements[i].variable_index_a == var;
                const bool is_b = obj.qelements[i].variable_index_b == var;

                if (is_a || is_b) {
                    quadratic_elements[id].factor =
                      static_cast<real>(obj.qelements[i].factor);
                    quadratic_elements[id].id =
                      is_a ? obj.qelements[i].variable_index_b
                           : obj.qelements[i].variable_index_a;
                    ++id;
                }
            }
        }
    }

    quadratic_cost_type(const quadratic_cost_type& other, int n)
      : obj(other.obj)
      , linear_elements(std::make_unique<real[]>(n))
      , quadratic_elements(std::make_unique<quad[]>(other.indices[n]))
      , indices(std::make_unique<int[]>(n + 1))
    {
        std::copy_n(other.linear_elements.get(), n, linear_elements.get());
        std::copy_n(other.indices.get(), n + 1, indices.get());
        std::copy_n(other.quadratic_elements.get(),
                    other.indices[n],
                    quadratic_elements.get());
    }

    template<typename Random>
    void make_random_norm(int n, Random& rng)
    {
        std::vector<std::pair<real, int>> r(n);
        std::vector<std::pair<quad, int>> q(indices[n]);

        for (int i = 0; i != n; ++i)
            r[i] = { linear_elements[i], i };

        std::sort(r.begin(), r.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.first < rhs.first;
        });

        {
            auto begin = r.begin();
            auto end = r.end();
            auto next = r.begin()++;
            for (; next != end; ++next) {
                if (next->first != begin->first) {
                    if (std::distance(begin, next) > 1)
                        random_epsilon_unique(
                          begin, next, rng, begin->first, next->first);

                    begin = next;
                }
            }

            if (std::distance(begin, end) > 1) {
                if (begin == r.begin()) {
                    random_epsilon_unique(
                      begin, end, rng, begin->first, begin->first + 1);
                } else {
                    auto value = linearize(std::ptrdiff_t(0),
                                           r.front().first,
                                           std::distance(r.begin(), begin),
                                           begin->first,
                                           std::distance(r.begin(), r.end()));

                    random_epsilon_unique(
                      begin, end, rng, begin->first, value);
                }
            }

            // Reorder the vector according to the variable index, so,
            // it restores the initial order.

            std::sort(
              r.begin(), r.end(), [](const auto& lhs, const auto& rhs) {
                  return lhs.second < rhs.second;
              });
        }

        {
            // auto begin = q.begin();
            // auto end = q.end();
            // auto next = q.begin()++;
            // for (; next != end; ++next) {
            //     if (next->first != begin->first) {
            //         if (std::distance(begin, next) > 1)
            //             random_epsilon_unique(
            //               begin, next, rng, begin->first,
            //               next->first);

            //         begin = next;
            //     }
            // }

            // if (std::distance(begin, end) > 1) {
            //     if (begin == q.begin()) {
            //         random_epsilon_unique(
            //           begin, end, rng, begin->first, begin->first +
            //           1);
            //     } else {
            //         auto value = linearize(std::ptrdiff_t(0),
            //                                q.front().first,
            //                                std::distance(q.begin(),
            //                                begin), begin->first,
            //                                std::distance(q.begin(),
            //                                q.end()));

            //         random_epsilon_unique(
            //           begin, end, rng, begin->first, value);
            //     }
            // }

            // Reorder the vector according to the variable index, so,
            // it restores the initial order.

            // std::sort(
            //   q.begin(), q.end(), [](const auto& lhs, const auto&
            //   rhs) {
            //       return lhs.second < rhs.second;
            //   });
        }

        for (int i = 0; i != n; ++i)
            linear_elements[i] = r[i].first;
        // for (int i = 0; i != n; ++i) {
        //     quadratic_elements[i].factor = r[i].first.factor;

        // Finally we compute the l+oo norm.
        make_loo_norm(n);
    }

    real min(int n) const noexcept
    {
        real min = std::numeric_limits<real>::max();

        for (int i = 0; i != n; ++i)
            if (linear_elements[i])
                min = std::min(min, std::abs(linear_elements[i]));

        for (int i = 0, e = indices[n]; i != e; ++i)
            if (quadratic_elements[i].factor)
                min = std::min(min, std::abs(quadratic_elements[i].factor));

        return min;
    }

    void make_l1_norm(int n)
    {
        real div = { 0 };

        for (int i = 0; i != n; ++i)
            div += std::abs(linear_elements[i]);

        for (int i = 0, e = indices[n]; i != e; ++i)
            div += std::abs(quadratic_elements[i].factor);

        if (std::isnormal(div)) {
            for (int i = 0; i != n; ++i)
                linear_elements[i] /= div;

            for (int i = 0, e = indices[n]; i != e; ++i)
                quadratic_elements[i].factor /= div;
        }
    }

    void make_l2_norm(int n)
    {
        real div = { 0 };

        for (int i = 0; i != n; ++i)
            div += linear_elements[i] * linear_elements[i];
        ;

        for (int i = 0, e = indices[n]; i != e; ++i)
            div += quadratic_elements[i].factor * quadratic_elements[i].factor;

        if (std::isnormal(div)) {
            for (int i = 0; i != n; ++i)
                linear_elements[i] /= div;

            for (int i = 0, e = indices[n]; i != e; ++i)
                quadratic_elements[i].factor /= div;
        }
    }

    void make_loo_norm(int n)
    {
        real div =
          *std::max_element(linear_elements.get(), linear_elements.get() + n);

        for (int i = 0, e = indices[n]; i != e; ++i)
            div = std::max(quadratic_elements[i].factor, div);

        if (std::isnormal(div)) {
            for (int i = 0; i != n; ++i)
                linear_elements[i] /= div;

            for (int i = 0, e = indices[n]; i != e; ++i)
                quadratic_elements[i].factor /= div;
        }
    }

    real operator[](int index) const noexcept
    {
        real ret = linear_elements[index];

        auto first = indices[index];
        auto last = indices[index + 1];

        for (; first != last; ++first)
            ret += quadratic_elements[first].factor;

        return ret;
    }

    real operator()(int index, const bit_array& x) const noexcept
    {
        real ret = linear_elements[index];

        auto first = indices[index];
        auto last = indices[index + 1];

        for (; first != last; ++first)
            if (x[quadratic_elements[first].id])
                ret += quadratic_elements[first].factor;

        return ret;
    }

    real results(const bit_array& x, real cost_constant) const noexcept
    {
        for (int i = 0, e = length(obj.elements); i != e; ++i)
            if (x[obj.elements[i].variable_index])
                cost_constant += static_cast<real>(obj.elements[i].factor);

        for (int i = 0, e = length(obj.qelements); i != e; ++i)
            if (x[obj.qelements[i].variable_index_a] &&
                x[obj.qelements[i].variable_index_b])
                cost_constant += static_cast<real>(obj.qelements[i].factor);

        return cost_constant;
    }
};

inline random_engine::result_type
init_random_generator_seed(const context& ctx) noexcept
{
    auto epoch = std::chrono::system_clock::now().time_since_epoch().count();
    auto param = ctx.parameters.seed;

    if (param <= 0)
        return static_cast<random_engine::result_type>(epoch);

    return static_cast<random_engine::result_type>(param);
}

/**
 * Generates an array with unique seed value.
 */
inline std::unique_ptr<random_engine::result_type[]>
generate_seed(random_engine& rng, unsigned thread_id)
{
    using type = random_engine::result_type;

    auto ret = std::make_unique<type[]>(thread_id);

    std::uniform_int_distribution<type> dst(std::numeric_limits<type>::min(),
                                            std::numeric_limits<type>::max());

    ret[0] = numeric_cast<type>(dst(rng));

    for (unsigned i = 1; i != thread_id; ++i) {
        ret[i] = numeric_cast<type>(dst(rng));

        unsigned j = i - 1;
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

template<int f>
using real_sel = typename std::conditional<
  f == 0,
  real,
  typename std::conditional<f == 1, double, long double>::type>::type;

template<int o>
using mode_sel =
  typename std::conditional<o == 0, maximize_tag, minimize_tag>::type;

template<bool debug>
struct debug_logger
{
    std::FILE* ofs = nullptr;

    debug_logger(const std::string_view name) noexcept
    {
        if constexpr (debug) {
            char buffer[512] = { '\0' };

            auto written = fmt::format_to_n(
              buffer,
              511,
              "{}-{}.log",
              name,
              std::hash<std::thread::id>{}(std::this_thread::get_id()));

            buffer[written.size] = '\0';

            ofs = std::fopen(buffer, "w");
        } else {
            (void)name;

            ofs = nullptr;
        }
    }

    ~debug_logger() noexcept
    {
        if constexpr (debug) {
            if (ofs)
                std::fclose(ofs);
        }
    }

    template<typename... Args>
    void log([[maybe_unused]] fmt::format_string<Args...> fmt,
             [[maybe_unused]] Args&&... args) const noexcept
    {
        if constexpr (debug) {
            fmt::print(ofs, fmt, std::forward<Args>(args)...);
        }
    }

    template<typename... Args>
    void log([[maybe_unused]] unsigned indent,
             [[maybe_unused]] fmt::format_string<Args...> fmt,
             [[maybe_unused]] Args&&... args) const noexcept
    {
        if constexpr (debug) {
            fmt::print(ofs, "{:{}}", "", indent);
            fmt::print(ofs, fmt, std::forward<Args>(args)...);
        }
    }
};

} // namespace itm
} // namespace baryonyx

#endif
