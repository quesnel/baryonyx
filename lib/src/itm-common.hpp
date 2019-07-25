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
#include <random>
#include <tuple>
#include <vector>

namespace baryonyx {
namespace itm {

/**
 * @brief stores vector solution in solver.
 * @details x_type uses a single @c std::vector<int8_t> instead of a
 *     @c baryonyx::fixed_array<bool> to use the specialized version of vector,
 *     which is used for elements of type bool and optimizes for space.
 *
 * @code
 */
struct x_type
{
    x_type(int size)
      : m_data(size)
    {}

    int operator[](int index) const noexcept
    {
        return m_data.get(index);
    }

    void invert(index index) noexcept
    {
        m_data.invert(index);
    }

    void set(index index)
    {
        m_data.set(index);
    }

    void unset(index index)
    {
        m_data.unset(index);
    }

    void set(index index, int v) noexcept
    {
        if (v > 0)
            m_data.set(index);
        else
            m_data.unset(index);
    }

    int upper() const noexcept
    {
        return 0;
    }

    void clear() const noexcept
    {}

    const bit_array& data() const noexcept
    {
        return m_data;
    }

private:
    bit_array m_data;
};

/**
 * @brief stores vector solution in solver and a counter.
 * @details x_type uses a single @c std::vector<int8_t> instead of a
 *     @c baryonyx::fixed_array<bool> to use the specialized version of vector,
 *     which is used for elements of type bool and optimizes for space. The
 *     counter vector is used to store number of change in vector solution
 *     during a @c update_row.
 */
struct x_counter_type
{
    x_counter_type(int size)
      : m_data(size)
      , m_counter(size, 0)
    {}

    int operator[](int index) const noexcept
    {
        return m_data.get(index);
    }

    void invert(index index) noexcept
    {
        // NOTE: only the data vector is updated. Normally, m_data and
        // m_counter have been already updated in the update_row function.

        m_data.invert(index);
    }

    void set(index index, int v) noexcept
    {
        // TODO: Maybe use integer class members to store upper and lower index
        // and make upper() and lower() function O(1).

        if (m_data.get(index) != v) {
            if (v)
                m_data.set(index);
            else
                m_data.unset(index);

            ++m_counter[index];
        }
    }

    void set(index index) noexcept
    {
        // TODO: Maybe use integer class members to store upper and lower index
        // and make upper() and lower() function O(1).

        if (!m_data.get(index)) {
            m_data.set(index);

            ++m_counter[index];
        }
    }

    void unset(index index) noexcept
    {
        // TODO: Maybe use integer class members to store upper and lower index
        // and make upper() and lower() function O(1).

        if (m_data.get(index)) {
            m_data.unset(index);

            ++m_counter[index];
        }
    }

    void clear() noexcept
    {
        std::fill(m_counter.begin(), m_counter.end(), 0);
    }

    const bit_array& data() const noexcept
    {
        return m_data;
    }

    int upper() const noexcept
    {
        int upper_index = 0;
        for (int i = 1, e = length(m_counter); i != e; ++i)
            if (m_counter[i] > m_counter[upper_index])
                upper_index = i;

        return upper_index;
    }

private:
    bit_array m_data;
    std::vector<int> m_counter;
};

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

template<typename Xtype, typename Cost>
double
results(const Xtype& x, const Cost& c, double cost_constant, int n)
{
    for (int i = 0; i != n; ++i)
        cost_constant += static_cast<double>(c[i] * x[i]);

    return cost_constant;
}

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

template<typename Mode, typename Iterator, typename Random>
inline void
calculator_sort(Iterator begin, Iterator end, Random& rng)
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

template<typename Mode, typename Float, typename Random>
inline bool
stop_iterating(Float value, Random& rng) noexcept
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
        slv.pi[k] += slv.R[0].value / two;
        d += (kappa / (one - kappa)) * (slv.R[0].value / two);

        for (int i = 0; i != r_size; ++i) {
            auto var = it + slv.R[i].id;

            if (slv.R[i].is_negative()) {
                x.set(var->column);
                slv.P[var->value] += d;
            } else {
                x.unset(var->column);
                slv.P[var->value] -= d;
            }
        }
    } else if (selected + 1 >= r_size) {
        slv.pi[k] += slv.R[selected].value * middle;
        d += (kappa / (one - kappa)) * (slv.R[selected].value * middle);

        for (int i = 0; i != r_size; ++i) {
            auto var = it + slv.R[i].id;

            if (slv.R[i].is_negative()) {
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

            if (slv.R[i].is_negative()) {
                x.unset(var->column);
                slv.P[var->value] -= d;
            } else {
                x.set(var->column);
                slv.P[var->value] += d;
            }
        }

        for (; i != r_size; ++i) {
            auto var = it + slv.R[i].id;

            if (slv.R[i].is_negative()) {
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

template<typename Solver, typename Float, typename Mode, typename Xtype>
class solver_initializer
{
    std::bernoulli_distribution dist;
    std::bernoulli_distribution toss_up{ 0.5 };
    // std::uniform_real_distribution<double> real_distribution{ 0, 1 };
    double old_best;

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
        random,
        bastert_zero_random_one,
        zero_random_one_bastert,
        random_one_bastert_zero,
        one_bastert_zero_random
    };

    single_automaton single;
    cycle_automaton cycle;

    static constexpr short iteration_limit = 3;

    void reinit(Solver& slv, Xtype& x) noexcept
    {
        for (int i = 0; i != slv.n; ++i)
            if (dist(slv.rng))
                x.set(i, 1 - x[i]);
    }

    void reinit(Solver& slv, Xtype& x, const solution& best) noexcept
    {
        for (int i = 0; i != slv.n; ++i)
            if (dist(slv.rng))
                x.set(i, 1 - best.variables[i]);
    }

    void init_zeros(Solver& slv, Xtype& x) noexcept
    {
        for (int i = 0; i != slv.n; ++i)
            x.set(i, 0);
    }

    void init_ones(Solver& slv, Xtype& x) noexcept
    {
        for (int i = 0; i != slv.n; ++i)
            x.set(i, 1);
    }

    // Fully fill the solution vector X with random values from user's
    // parameters.
    void init_random(Solver& slv, Xtype& x) noexcept
    {
        for (int i = 0; i != slv.n; ++i)
            x.set(i, dist(slv.rng));
    }

    // Fully fill the solution vector X with a mix of values from the bastert
    // policy and from random value.
    void init_bastert(Solver& slv, Xtype& x) noexcept
    {
        const int value_if_cost_0 = 1;

        for (int i = 0; i != slv.n; ++i) {
            if (dist(slv.rng))
                x.set(i, init_x<Mode>(slv.c[i], value_if_cost_0));
            else
                x.set(i, toss_up(slv.rng));
        }
    }

public:
    solver_initializer(Solver& slv,
                       Xtype& x,
                       solver_parameters::init_policy_type policy,
                       double init_random)
      : dist(init_random)
      , old_best(bad_value<Mode, double>())
      , single(single_automaton::pre_init)
      , cycle(cycle_automaton::bastert)
    {
        x.clear();
        slv.reset();

        switch (policy) {
        case solver_parameters::init_policy_type::bastert:
        case solver_parameters::init_policy_type::bastert_cycle:
            cycle = cycle_automaton::bastert;
            break;
        case solver_parameters::init_policy_type::random:
        case solver_parameters::init_policy_type::random_cycle:
            cycle = cycle_automaton::random;
            break;
        case solver_parameters::init_policy_type::best:
        case solver_parameters::init_policy_type::best_cycle:
            cycle = cycle_automaton::bastert_zero_random_one;
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
                    if (cycle == cycle_automaton::bastert_zero_random_one)
                        cycle = cycle_automaton::zero_random_one_bastert;
                    else if (cycle == cycle_automaton::zero_random_one_bastert)
                        cycle = cycle_automaton::random_one_bastert_zero;
                    else if (cycle == cycle_automaton::random_one_bastert_zero)
                        cycle = cycle_automaton::one_bastert_zero_random;
                    else if (cycle == cycle_automaton::one_bastert_zero_random)
                        cycle = cycle_automaton::bastert_zero_random_one;
                }
            }
            break;
        }
    }

    void reinit(Solver& slv, Xtype& x, bool is_x_solution, const result& best)
    {
        x.clear();
        slv.reset();

        auto is_x_better_solution{ false };

        if (!best.solutions.empty() &&
            is_better_solution<Mode, double>(old_best,
                                             best.solutions.back().value)) {
            old_best = best.solutions.back().value;
            is_x_better_solution = true;
        }

        next_state(is_x_solution, is_x_better_solution);

        switch (single) {
        case single_automaton::pre_init:
        case single_automaton::init:
            switch (cycle) {
            case cycle_automaton::bastert:
            case cycle_automaton::bastert_zero_random_one:
                init_bastert(slv, x);
                return;
            case cycle_automaton::random:
            case cycle_automaton::random_one_bastert_zero:
                init_random(slv, x);
                break;
            case cycle_automaton::zero_random_one_bastert:
                init_zeros(slv, x);
                break;
            case cycle_automaton::one_bastert_zero_random:
                init_ones(slv, x);
                break;
            }
            break;
        case single_automaton::improve_x:
            reinit(slv, x);
            break;
        case single_automaton::improve_best:
            reinit(slv, x, best.solutions.back());
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

struct compute_incremental_order
{};

struct compute_decremental_order
{};

template<typename floatingpointT, typename randomT>
struct compute_none
{
    using random_type = randomT;

    std::vector<int> R;

    template<typename solverT, typename Xtype>
    compute_none(const solverT& s, const Xtype& x, randomT&)
      : R(s.m)
    {
        compute_violated_constraints(s, x, R);
    }

    template<typename solverT, typename Xtype>
    int push_and_run(solverT& solver,
                     Xtype& x,
                     floatingpointT kappa,
                     floatingpointT delta,
                     floatingpointT theta,
                     floatingpointT objective_amplifier)
    {
        solver.push_and_compute_update_row(
          x, R.cbegin(), R.cend(), kappa, delta, theta, objective_amplifier);

        return compute_violated_constraints(solver, x, R);
    }

    template<typename solverT, typename Xtype>
    int run(solverT& solver,
            Xtype& x,
            floatingpointT kappa,
            floatingpointT delta,
            floatingpointT theta)
    {
        solver.compute_update_row(x, R.begin(), R.end(), kappa, delta, theta);

        return compute_violated_constraints(solver, x, R);
    }
};

template<typename Float, typename Random>
struct compute_pi_sign_change
{
#if 0
    std::vector<int> R;
    Random& rng;

    template<typename Solver, typename X>
    compute_pi_sign_change(const Solver& s, const X&, Random& rng_)
      : R(s.m)
      , rng(rng_)
    {
        std::iota(R.begin(), R.end(), 0);
    }

    template<typename solverT, typename Xtype>
    int push_and_run(solverT& solver,
                     Xtype& x,
                     Float kappa,
                     Float delta,
                     Float theta,
                     Float objective_amplifier)
    {
        std::shuffle(R.begin(), R.end(), rng);

        auto pi_changed = solver.push_and_compute_update_row(
          x, R.begin(), R.end(), kappa, delta, theta, objective_amplifier);

        auto remaining = compute_violated_constraints(solver, x, R);

        if (!pi_changed && remaining == 0)
            return 0;

        return compute_violated_constraints(solver, x, R);
    }

    template<typename Solver, typename X>
    int run(Solver& solver, X& x, Float kappa, Float delta, Float theta)
    {
        std::shuffle(R.begin(), R.end(), rng);

        auto pi_changed = solver.compute_update_row(
          x, R.begin(), R.end(), kappa, delta, theta);

        auto remaining = compute_violated_constraints(solver, x, R);

        if (!pi_changed && remaining == 0)
            return 0;

        return compute_violated_constraints(solver, x, R);
    }
#endif

    std::vector<int> R;
    Random& rng;

    template<typename Solver, typename X>
    compute_pi_sign_change(const Solver& s, const X&, Random& rng_)
      : R(s.m)
      , rng(rng_)
    {
        std::iota(R.begin(), R.end(), 0);
    }

    template<typename solverT, typename Xtype>
    int push_and_run(solverT& solver,
                     Xtype& x,
                     Float kappa,
                     Float delta,
                     Float theta,
                     Float objective_amplifier)
    {
        std::shuffle(R.begin(), R.end(), rng);

        auto pi_changed = solver.push_and_compute_update_row(
          x, R.begin(), R.end(), kappa, delta, theta, objective_amplifier);
        auto remaining = local_compute_violated_constraints(solver, x);

        if (!pi_changed && remaining == 0)
            return 0;

        return remaining;
    }

    template<typename Solver, typename X>
    int run(Solver& solver, X& x, Float kappa, Float delta, Float theta)
    {
        std::shuffle(R.begin(), R.end(), rng);

        auto pi_changed = solver.compute_update_row(
          x, R.begin(), R.end(), kappa, delta, theta);
        auto remaining = local_compute_violated_constraints(solver, x);

        if (!pi_changed && remaining == 0)
            return 0;

        return remaining;
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
};

template<typename floatingpointT, typename randomT>
struct compute_reversing
{
    std::vector<int> R;

    template<typename solverT, typename Xtype>
    compute_reversing(solverT& s, const Xtype& x, randomT&)
      : R(s.m)
    {
        compute_violated_constraints(s, x, R);
    }

    template<typename solverT, typename Xtype>
    int push_and_run(solverT& solver,
                     Xtype& x,
                     floatingpointT kappa,
                     floatingpointT delta,
                     floatingpointT theta,
                     floatingpointT objective_amplifier)
    {
        solver.push_and_compute_update_row(
          x, R.crbegin(), R.crend(), kappa, delta, theta, objective_amplifier);

        return compute_violated_constraints(solver, x, R);
    }

    template<typename solverT, typename Xtype>
    int run(solverT& solver,
            Xtype& x,
            floatingpointT kappa,
            floatingpointT delta,
            floatingpointT theta)
    {
        solver.compute_update_row(
          x, R.crbegin(), R.crend(), kappa, delta, theta);

        return compute_violated_constraints(solver, x, R);
    }
};

template<typename Float, typename Random, typename Order>
struct compute_lagrangian_order
{
    std::vector<int> R;

    template<typename solverT, typename Xtype>
    compute_lagrangian_order(solverT& s, const Xtype& x, Random&)
      : R(s.m)
    {
        compute_violated_constraints(s, x, R);
    }

    template<typename Iterator, typename Solver>
    static void local_sort(Iterator begin, Iterator end, const Solver& solver)
    {
        std::sort(begin, end, [&solver](const auto lhs, const auto rhs) {
            if constexpr (std::is_same_v<Order, compute_incremental_order>)
                return solver.pi[lhs] < solver.pi[rhs];
            else
                return solver.pi[rhs] < solver.pi[lhs];
        });
    }

    template<typename solverT, typename Xtype>
    int push_and_run(solverT& solver,
                     Xtype& x,
                     Float kappa,
                     Float delta,
                     Float theta,
                     Float objective_amplifier)
    {
        local_sort(R.begin(), R.end(), solver);

        solver.push_and_compute_update_row(
          x, R.cbegin(), R.cend(), kappa, delta, theta, objective_amplifier);

        return compute_violated_constraints(solver, x, R);
    }

    template<typename solverT, typename Xtype>
    int run(solverT& solver, Xtype& x, Float kappa, Float delta, Float theta)
    {
        local_sort(R.begin(), R.end(), solver);

        solver.compute_update_row(
          x, R.cbegin(), R.cend(), kappa, delta, theta);

        return compute_violated_constraints(solver, x, R);
    }
};

template<typename floatingpointT, typename randomT>
struct compute_random
{
    using random_type = randomT;

    std::vector<int> R;
    random_type& rng;

    template<typename solverT, typename Xtype>
    compute_random(solverT& s, const Xtype& x, random_type& rng_)
      : R(s.m)
      , rng(rng_)
    {
        compute_violated_constraints(s, x, R);
    }

    template<typename solverT, typename Xtype>
    int push_and_run(solverT& solver,
                     Xtype& x,
                     floatingpointT kappa,
                     floatingpointT delta,
                     floatingpointT theta,
                     floatingpointT objective_amplifier)
    {
        std::shuffle(R.begin(), R.end(), rng);

        solver.push_and_compute_update_row(
          x, R.begin(), R.end(), kappa, delta, theta, objective_amplifier);

        return compute_violated_constraints(solver, x, R);
    }

    template<typename solverT, typename Xtype>
    int run(solverT& solver,
            Xtype& x,
            floatingpointT kappa,
            floatingpointT delta,
            floatingpointT theta)
    {
        std::shuffle(R.begin(), R.end(), rng);

        solver.compute_update_row(x, R.begin(), R.end(), kappa, delta, theta);

        return compute_violated_constraints(solver, x, R);
    }
};

template<typename floatingpointT, typename randomT, typename Order>
struct compute_infeasibility
{
    using random_type = randomT;

    std::vector<std::pair<int, int>> m_order;
    random_type& rng;

    template<typename solverT, typename Xtype>
    compute_infeasibility(solverT& s, const Xtype& x, random_type& rng_)
      : m_order(s.m)
      , rng(rng_)
    {
        local_compute_violated_constraints(s, x);
    }

    template<typename solverT, typename Xtype>
    int local_compute_violated_constraints(solverT& solver, const Xtype& x)
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

    template<typename Iterator>
    static void local_sort(Iterator begin, Iterator end) noexcept
    {
        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            if constexpr (std::is_same_v<Order, compute_incremental_order>)
                return lhs.second < rhs.second;
            else
                return rhs.second < lhs.second;
        });
    }

    template<typename solverT, typename Xtype>
    int push_and_run(solverT& solver,
                     Xtype& x,
                     floatingpointT kappa,
                     floatingpointT delta,
                     floatingpointT theta,
                     floatingpointT objective_amplifier)
    {
        local_sort(m_order.begin(), m_order.end());

        solver.push_and_compute_update_row(x,
                                           m_order.cbegin(),
                                           m_order.cend(),
                                           kappa,
                                           delta,
                                           theta,
                                           objective_amplifier);

        return local_compute_violated_constraints(solver, x);
    }

    template<typename solverT, typename Xtype>
    int run(solverT& solver,
            Xtype& x,
            floatingpointT kappa,
            floatingpointT delta,
            floatingpointT theta)
    {
        local_sort(m_order.begin(), m_order.end());

        solver.compute_update_row(
          x, m_order.cbegin(), m_order.cend(), kappa, delta, theta);

        return local_compute_violated_constraints(solver, x);
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
        if (c[i] != 0 && std::abs(c[i]) < mini)
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
    bx_expects(min != max);

    std::uniform_real_distribution<floatingpointT> distribution(min, max);

    for (; begin != end; ++begin)
        begin->first = distribution(rng);
}

template<typename floatingpointT, typename randomT>
inline std::unique_ptr<floatingpointT[]>
rng_normalize_costs(const std::unique_ptr<floatingpointT[]>& c,
                    randomT& rng,
                    int n)
{
    std::vector<std::pair<floatingpointT, int>> r(n);

    for (int i = 0; i != n; ++i)
        r[i] = { c[i], i };

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
    std::copy_n(c.get(), n, ret.get());

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
        bx_ensures(0 <= n && elem.variable_index < n);
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
generate_seed(randomT& rng, unsigned thread_id)
{
    using type = typename randomT::result_type;

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

template<typename Float, typename Random, int o>
using constraint_sel = typename std::conditional<
  o == 0,
  compute_none<Float, Random>,
  typename std::conditional<
    o == 1,
    compute_reversing<Float, Random>,
    typename std::conditional<
      o == 2,
      compute_random<Float, Random>,
      typename std::conditional<
        o == 3,
        compute_infeasibility<Float, Random, compute_decremental_order>,
        typename std::conditional<
          o == 4,
          compute_infeasibility<Float, Random, compute_incremental_order>,
          typename std::conditional<
            o == 5,
            compute_lagrangian_order<Float, Random, compute_decremental_order>,
            typename std::conditional<
              o == 6,
              compute_lagrangian_order<Float,
                                       Random,
                                       compute_incremental_order>,
              compute_pi_sign_change<Float, Random>>::type>::type>::type>::
        type>::type>::type>::type;

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
