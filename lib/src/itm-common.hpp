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
#include <iterator>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>

namespace baryonyx {
namespace itm {

using random_engine = std::default_random_engine;

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
make_merged_constraints(const context_ptr& ctx, const problem& pb);

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

template<typename Mode, typename Float>
inline bool
stop_iterating(Float value, random_engine& rng) noexcept
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

template<typename Mode, typename Float>
constexpr Float
bad_value() noexcept
{
    if constexpr (std::is_same_v<Mode, minimize_tag>)
        return std::numeric_limits<Float>::infinity();
    else
        return -std::numeric_limits<Float>::infinity();
}

template<typename Mode, typename Float>
constexpr bool
is_bad_value(Float value) noexcept
{
    return value == bad_value<Mode, Float>();
}

template<typename Mode, typename Float>
inline bool
stop_iterating(Float value) noexcept
{
    if constexpr (std::is_same_v<Mode, minimize_tag>)
        return value > 0;
    else
        return value < 0;
}

template<typename Mode, typename Float>
inline bool
is_better_solution(Float lhs, Float rhs) noexcept
{
    if constexpr (std::is_same_v<Mode, minimize_tag>)
        return lhs < rhs;
    else
        return lhs > rhs;
}

template<typename Mode, typename Float>
inline bool
init_x(Float cost, int value_if_cost_0) noexcept
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
template<typename Float>
inline bool
is_signbit_change(Float lhs, Float rhs) noexcept
{
    return std::signbit(lhs) != std::signbit(rhs);
}

/// The bkmin and bkmax constraint bounds are not equal and can be assigned to
/// -infinity or +infinity. We have to scan the r vector and search a value j
/// such as b(0, k) <= Sum A(k, R[j]) < b(1, k).
///
/// @return True if the sign of the pi vector change during the affect
///     operation.
template<typename Solver, typename Xtype, typename Iterator, typename Float>
bool
affect(Solver& slv,
       Xtype& x,
       Iterator it,
       int k,
       int selected,
       int r_size,
       const Float kappa,
       const Float delta)
{
    constexpr Float one{ 1 };
    constexpr Float two{ 2 };
    constexpr Float middle{ (two + one) / two };

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

template<typename Solver, typename Float, typename Mode>
class solver_initializer
{
    std::bernoulli_distribution dist;
    std::bernoulli_distribution toss_up;
    double old_best;
    bool use_cycle{ true };
    std::uint8_t iteration{ 0 };

    enum class single_automaton : std::uint8_t
    {
        pre_init,
        init,
        improve_x,
        improve_best
    };

    enum class cycle_automaton : std::uint8_t
    {
        bastert,
        pessimistic_solve,
        optimistic_solve,
    };

    static inline const std::string_view single_automaton_string[] =
      { "pre_init", "init", "improve_x", "improve_best" };

    static inline const std::string_view cycle_automaton_string[] = {
        "bastert",
        "pessimistic_solve",
        "optimistic_solve",
    };

    single_automaton single{ single_automaton::pre_init };
    cycle_automaton cycle{ cycle_automaton::pessimistic_solve };

    static constexpr short iteration_limit = 3;

    // Fully fill the solution vector X with a mix of values from the bastert
    // policy and from random value.

    void init_bastert(Solver& slv, bit_array& x) noexcept
    {
        const int value_if_cost_0 = 1;

        for (int i = 0; i != slv.n; ++i) {
            if (dist(slv.rng)) {
                if (init_x<Mode>(slv.c(i, x), value_if_cost_0))
                    x.set(i);
                else
                    x.unset(i);
            } else {
                if (toss_up(slv.rng))
                    x.set(i);
                else
                    x.unset(i);
            }
        }
    }

    void reinit_bastert(Solver& slv,
                        bit_array& x,
                        const solution& best) noexcept
    {
        const int value_if_cost_0 = 1;

        for (int i = 0; i != slv.n; ++i)
            if (dist(slv.rng)) {
                if (init_x<Mode>(slv.c(i, x), value_if_cost_0))
                    x.set(i);
                else
                    x.unset(i);
            } else {
                if (best.variables[i])
                    x.set(i);
                else
                    x.unset(i);
            }
    }

    void reinit_bastert(Solver& slv, bit_array& x) noexcept
    {
        const int value_if_cost_0 = 1;

        for (int i = 0; i != slv.n; ++i)
            if (dist(slv.rng)) {
                if (init_x<Mode>(slv.c(i, x), value_if_cost_0))
                    x.set(i);
                else
                    x.unset(i);
            }
    }

    void do_init_pessimistic_solve(Solver& slv, bit_array& x) noexcept
    {
        for (int k = 0; k != slv.m; ++k) {
            if (!dist(slv.rng))
                continue;
            auto [begin, end] = slv.ap.row(k);
            int r_size = 0;

            for (auto it = begin; it != end; ++it, ++r_size) {
                slv.R[r_size].value = slv.c(it->column, x);
                slv.R[r_size].id = r_size;
            }

            std::shuffle(slv.R.get(), slv.R.get() + r_size, slv.rng);

            std::sort(slv.R.get(),
                      slv.R.get() + r_size,
                      [](const auto& lhs, const auto& rhs) {
                          if constexpr (std::is_same_v<Mode, minimize_tag>)
                              return lhs.value < rhs.value;
                          else
                              return rhs.value < lhs.value;
                      });

            int sum = 0;
            int best = -2;

            for (int i = -1; i < r_size; ++i) {
                if (slv.bound_min(k) <= sum && sum <= slv.bound_max(k)) {
                    best = i;
                    break;
                }

                auto var = begin + slv.R[i].id;
                sum += slv.factor(var->value);
            }

            int i = 0;
            for (; i <= best; ++i) {
                auto var = begin + slv.R[i].id;
                x.set(var->column);
            }

            for (; i != r_size; ++i) {
                auto var = begin + slv.R[i].id;
                x.unset(var->column);
            }
        }
    }

    void init_pessimistic_solve(Solver& slv, bit_array& x) noexcept
    {
        for (int i = 0; i != slv.n; ++i)
            if (toss_up(slv.rng))
                x.set(i);
            else
                x.unset(i);

        do_init_pessimistic_solve(slv, x);
    }

    void reinit_pessimistic_solve(Solver& slv, bit_array& x) noexcept
    {
        do_init_pessimistic_solve(slv, x);
    }

    void reinit_pessimistic_solve(Solver& slv,
                                  bit_array& x,
                                  const solution& best) noexcept
    {
        for (int i = 0; i != slv.n; ++i)
            if (best.variables[i])
                x.set(i);
            else
                x.unset(i);

        do_init_pessimistic_solve(slv, x);
    }

    /////

    void do_init_optimistic_solve(Solver& slv, bit_array& x) noexcept
    {
        for (int k = 0; k != slv.m; ++k) {
            if (!dist(slv.rng))
                continue;

            auto [begin, end] = slv.ap.row(k);
            int r_size = 0;

            for (auto it = begin; it != end; ++it, ++r_size) {
                slv.R[r_size].value = slv.c(it->column, x);
                slv.R[r_size].id = r_size;
            }

            std::shuffle(slv.R.get(), slv.R.get() + r_size, slv.rng);

            std::sort(slv.R.get(),
                      slv.R.get() + r_size,
                      [](const auto& lhs, const auto& rhs) {
                          if constexpr (std::is_same_v<Mode, minimize_tag>)
                              return lhs.value < rhs.value;
                          else
                              return rhs.value < lhs.value;
                      });

            int sum = 0;
            int best = -2;

            for (int i = -1; i < r_size; ++i) {
                if (slv.bound_min(k) <= sum && sum <= slv.bound_max(k))
                    best = i;

                if (best != -2 && i + 1 < r_size &&
                    stop_iterating<Mode>(slv.R[i + 1].value))
                    break;

                auto var = begin + slv.R[i].id;
                sum += slv.factor(var->value);
            }

            int i = 0;
            for (; i <= best; ++i) {
                auto var = begin + slv.R[i].id;
                x.set(var->column);
            }

            for (; i != r_size; ++i) {
                auto var = begin + slv.R[i].id;
                x.unset(var->column);
            }
        }
    }

    void init_optimistic_solve(Solver& slv, bit_array& x) noexcept
    {
        for (int i = 0; i != slv.n; ++i)
            if (toss_up(slv.rng))
                x.set(i);
            else
                x.unset(i);

        do_init_optimistic_solve(slv, x);
    }

    void reinit_optimistic_solve(Solver& slv, bit_array& x) noexcept
    {
        do_init_optimistic_solve(slv, x);
    }

    void reinit_optimistic_solve(Solver& slv,
                                 bit_array& x,
                                 const solution& best) noexcept
    {
        for (int i = 0; i != slv.n; ++i)
            if (best.variables[i])
                x.set(i);
            else
                x.unset(i);

        do_init_optimistic_solve(slv, x);
    }

public:
    solver_initializer(Solver& slv,
                       bit_array& /*x*/,
                       solver_parameters::init_policy_type policy,
                       double init_policy_random,
                       double init_random)
      : dist(init_policy_random)
      , toss_up(init_random)
      , old_best(bad_value<Mode, double>())
    {
        // x.clear();
        slv.reset();

        switch (policy) {
        case solver_parameters::init_policy_type::bastert:
            cycle = cycle_automaton::bastert;
            use_cycle = false;
            break;
        case solver_parameters::init_policy_type::pessimistic_solve:
            cycle = cycle_automaton::pessimistic_solve;
            use_cycle = false;
            break;
        case solver_parameters::init_policy_type::optimistic_solve:
            cycle = cycle_automaton::optimistic_solve;
            use_cycle = false;
            break;
        case solver_parameters::init_policy_type::cycle:
            cycle = cycle_automaton::pessimistic_solve;
            use_cycle = true;
            break;
        }
    }

    void next_state(bool is_x_solution, bool is_x_better_solution)
    {
        switch (single) {
        case single_automaton::pre_init:
            single = single_automaton::init;
            break;
        case single_automaton::init:
            iteration = 0;
            single = is_x_solution && is_x_better_solution
                       ? single_automaton::improve_best
                       : single_automaton::improve_x;
            break;
        case single_automaton::improve_x:
        case single_automaton::improve_best:
            if (is_x_better_solution) {
                iteration = 0;
                single = single_automaton::improve_best;
            } else {
                ++iteration;
                if (iteration >= iteration_limit) {
                    single = single_automaton::pre_init;

                    if (use_cycle == true) {
                        if (cycle == cycle_automaton::bastert)
                            cycle = cycle_automaton::optimistic_solve;
                        else if (cycle == cycle_automaton::pessimistic_solve)
                            cycle = cycle_automaton::bastert;
                        else if (cycle == cycle_automaton::optimistic_solve)
                            cycle = cycle_automaton::pessimistic_solve;
                    }
                }
            }
            break;
        }
    }

    void reinit(Solver& slv,
                bit_array& x,
                bool is_x_solution,
                const result& best)
    {
        // x.clear();
        slv.reset();

        auto is_x_better_solution{ false };

        if (!best.solutions.empty() &&
            is_better_solution<Mode, double>(best.solutions.back().value,
                                             old_best)) {
            old_best = best.solutions.back().value;
            is_x_better_solution = true;
        }

        next_state(is_x_solution, is_x_better_solution);

        to_log(stdout,
               "  solver initialization: method {} in mode {}\n",
               cycle_automaton_string[static_cast<int>(cycle)],
               single_automaton_string[static_cast<int>(single)]);

        switch (single) {
        case single_automaton::pre_init:
        case single_automaton::init:
            switch (cycle) {
            case cycle_automaton::bastert:
                init_bastert(slv, x);
                return;
            case cycle_automaton::pessimistic_solve:
                init_pessimistic_solve(slv, x);
                break;
            case cycle_automaton::optimistic_solve:
                init_optimistic_solve(slv, x);
                break;
            }
            break;
        case single_automaton::improve_x:
            switch (cycle) {
            case cycle_automaton::bastert:
                reinit_bastert(slv, x);
                return;
            case cycle_automaton::pessimistic_solve:
                reinit_pessimistic_solve(slv, x);
                break;
            case cycle_automaton::optimistic_solve:
                reinit_optimistic_solve(slv, x);
                break;
            }
            break;
        case single_automaton::improve_best:
            switch (cycle) {
            case cycle_automaton::bastert:
                reinit_bastert(slv, x, best.solutions.back());
                return;
            case cycle_automaton::pessimistic_solve:
                reinit_pessimistic_solve(slv, x, best.solutions.back());
                break;
            case cycle_automaton::optimistic_solve:
                reinit_optimistic_solve(slv, x, best.solutions.back());
                break;
            }
            break;
        }
    }
};

template<typename Solver, typename Xtype>
inline void
print_solver(const Solver& slv,
             const Xtype& x,
             const context_ptr& ctx,
             const std::vector<std::string>& names,
             int print_level)
{
    if (print_level <= 0)
        return;

    debug(ctx, "  - X: {} to {}\n", 0, slv.n);
    for (int i = 0; i != slv.n; ++i)
        debug(ctx, "    - {} {}={}/c_i:{}\n", i, names[i], x[i], slv.c[i]);
    debug(ctx, "\n");

    for (int k = 0; k != slv.m; ++k) {
        sparse_matrix<int>::const_row_iterator it, et;

        std::tie(it, et) = slv.ap.row(k);
        int v = 0;

        for (; it != et; ++it)
            v += slv.factor(it->value) * x[it->column];

        bool valid = slv.bound_min(k) <= v && v <= slv.bound_max(k);

        debug(ctx,
              "C {}:{} (Lmult: {})\n",
              k,
              (valid ? "   valid" : "violated"),
              slv.pi[k]);
    }
}

template<typename Solver, typename Xtype>
inline void
print_missing_constraint(const context_ptr& ctx,
                         const Xtype& x,
                         const Solver& slv,
                         const std::vector<std::string>& names) noexcept
{
    std::vector<int> R;

    slv.compute_violated_constraints(x, R);
    info(ctx, "Constraints remaining: {}\n", length(R));

    sparse_matrix<int>::const_row_iterator it, et;

    for (auto k : R) {
        std::tie(it, et) = slv.ap.row(k);
        int v = 0;

        info(ctx, "{}: {} <= ", k, slv.bound_min(k));

        for (; it != et; ++it) {
            v += slv.factor(it->value) * x[it->column];

            info(ctx,
                 "{:+d} [{} ({})] ",
                 slv.factor(it->value),
                 names[it->column],
                 x[it->column]);
        }

        info(ctx, " <= {} | value: {}\n", slv.bound_max(k));
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
                    const result& best)
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

    template<typename Solver, typename Xtype, typename Float>
    int push_and_run(Solver& solver,
                     Xtype& x,
                     random_engine& rng,
                     Float kappa,
                     Float delta,
                     Float theta,
                     Float objective_amplifier)
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

    template<typename Solver, typename Xtype, typename Float>
    int run(Solver& solver,
            Xtype& x,
            random_engine& rng,
            Float kappa,
            Float delta,
            Float theta)
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

template<typename Float, typename Cost>
inline Float
compute_delta(const context_ptr& ctx, const Cost& c, Float theta, int n)
{
    info(ctx, "  - delta not defined, compute it:\n");

    const auto mini = c.min(n);
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

template<typename floatingpointT, typename iteratorT>
inline void
random_epsilon_unique(iteratorT begin,
                      iteratorT end,
                      random_engine& rng,
                      floatingpointT min,
                      floatingpointT max)
{
    bx_expects(min != max);

    std::uniform_real_distribution<floatingpointT> distribution(min, max);

    for (; begin != end; ++begin)
        begin->first = distribution(rng);
}

/**
 * Normalizes the cost vector, i.e. divides it by its l{1,2, +oo}norm.
 * If the input vector is too small or with infinity value, the c is
 * unchanged.
 */
template<typename floatingpointT, typename Cost>
inline Cost
normalize_costs(const context_ptr& ctx,
                const Cost& c,
                random_engine& rng,
                int n)
{
    Cost ret(c, n);

    switch (ctx->parameters.cost_norm) {
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

template<typename Float>
struct default_cost_type
{
    const baryonyx::objective_function& obj;
    std::unique_ptr<Float[]> linear_elements;

    default_cost_type(const objective_function& obj_, int n)
      : obj(obj_)
      , linear_elements(std::make_unique<Float[]>(n))
    {
        for (const auto& elem : obj.elements) {
            bx_ensures(0 <= elem.variable_index && elem.variable_index < n);

            linear_elements[elem.variable_index] +=
              static_cast<Float>(elem.factor);
        }
    }

    default_cost_type(const default_cost_type& other, int n)
      : obj(other.obj)
      , linear_elements(std::make_unique<Float[]>(n))
    {
        std::copy_n(other.linear_elements.get(), n, linear_elements.get());
    }

    template<typename Random>
    void make_random_norm(int n, Random& rng)
    {
        std::vector<std::pair<Float, int>> r(n);

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

        Float div =
          *std::max_element(linear_elements.get(), linear_elements.get() + n);
        if (std::isnormal(div))
            for (int i = 0; i != n; ++i)
                linear_elements[i] /= div;
    }

    Float min(int n) const noexcept
    {
        Float min = std::numeric_limits<Float>::max();

        for (int i = 0; i != n; ++i)
            if (linear_elements[i])
                min = std::min(min, std::abs(linear_elements[i]));

        return min;
    }

    void make_l1_norm(int n)
    {
        Float div = { 0 };

        for (int i = 0; i != n; ++i)
            div += std::abs(linear_elements[i]);

        if (std::isnormal(div))
            for (int i = 0; i != n; ++i)
                linear_elements[i] /= div;
    }

    void make_l2_norm(int n)
    {
        Float div = { 0 };

        for (int i = 0; i != n; ++i)
            div += linear_elements[i] * linear_elements[i];

        if (std::isnormal(div))
            for (int i = 0; i != n; ++i)
                linear_elements[i] /= div;
    }

    void make_loo_norm(int n)
    {
        Float div =
          *std::max_element(linear_elements.get(), linear_elements.get() + n);

        if (std::isnormal(div))
            for (int i = 0; i != n; ++i)
                linear_elements[i] /= div;
    }

    template<typename X>
    Float operator()(int index, [[maybe_unused]] const X& x) const noexcept
    {
        return linear_elements[index];
    }

    template<typename Xtype>
    double results(const Xtype& x, double cost_constant) const noexcept
    {
        for (int i = 0, e = length(obj.elements); i != e; ++i)
            if (x[obj.elements[i].variable_index])
                cost_constant += obj.elements[i].factor;

        return cost_constant;
    }
};

template<typename Float>
struct quadratic_cost_type
{
    struct quad
    {
        Float factor;
        int id;
    };

    const baryonyx::objective_function& obj;
    std::unique_ptr<Float[]> linear_elements;
    std::unique_ptr<quad[]> quadratic_elements;
    std::unique_ptr<int[]> indices;

    quadratic_cost_type(const objective_function& obj_, int n)
      : obj(obj_)
      , linear_elements(std::make_unique<Float[]>(n))
      , quadratic_elements(
          std::make_unique<quad[]>(2 * obj.qelements.size() + 1))
      , indices(std::make_unique<int[]>(n + 1))
    {
        for (int i = 0, e = length(obj.elements); i != e; ++i) {
            bx_ensures(0 <= obj.elements[i].variable_index &&
                       obj.elements[i].variable_index < n);

            linear_elements[obj.elements[i].variable_index] +=
              static_cast<Float>(obj.elements[i].factor);
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
                      static_cast<Float>(obj.qelements[i].factor);
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
      , linear_elements(std::make_unique<Float[]>(n))
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
        std::vector<std::pair<Float, int>> r(n);
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

    Float min(int n) const noexcept
    {
        Float min = std::numeric_limits<Float>::max();

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
        Float div = { 0 };

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
        Float div = { 0 };

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
        Float div =
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

    template<typename X>
    Float operator()(int index, const X& x) const noexcept
    {
        Float ret = linear_elements[index];

        auto first = indices[index];
        auto last = indices[index + 1];

        for (; first != last; ++first)
            if (x[quadratic_elements[first].id])
                ret += quadratic_elements[first].factor;

        return ret;
    }

    template<typename Xtype>
    double results(const Xtype& x, double cost_constant) const noexcept
    {
        for (int i = 0, e = length(obj.elements); i != e; ++i)
            if (x[obj.elements[i].variable_index])
                cost_constant += obj.elements[i].factor;

        for (int i = 0, e = length(obj.qelements); i != e; ++i)
            if (x[obj.qelements[i].variable_index_a] &&
                x[obj.qelements[i].variable_index_b])
                cost_constant += obj.qelements[i].factor;

        return cost_constant;
    }
};

inline random_engine::result_type
init_random_generator_seed(const context_ptr& ctx) noexcept
{
    auto epoch = std::chrono::system_clock::now().time_since_epoch().count();
    auto param = ctx->parameters.seed;

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
using float_sel = typename std::conditional<
  f == 0,
  float,
  typename std::conditional<f == 1, double, long double>::type>::type;

template<int o>
using mode_sel =
  typename std::conditional<o == 0, maximize_tag, minimize_tag>::type;

} // namespace itm
} // namespace baryonyx

#endif
