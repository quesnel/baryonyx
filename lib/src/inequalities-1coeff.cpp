/* Copyright (C) 2017 INRA
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublipnse, and/or sell copies of the Software, and to
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

#include <baryonyx/core-compare>
#include <baryonyx/core-out>

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

#include "fixed_array.hpp"
#include "itm.hpp"
#include "matrix.hpp"
#include "private.hpp"
#include "utils.hpp"

#include <cassert>
#include <utility>

namespace bx = baryonyx;

namespace {


struct maximize_tag
{
};
struct minimize_tag
{
};

struct bound
{
    bound() = default;

    bound(int min_, int max_)
      : min(min_)
      , max(max_)
    {
    }

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
    {
    }

    floatingpointT value;
    int id;
};

template<typename floatingpointT>
using AP_type = bx::SparseArray<std::int8_t, floatingpointT>;

using b_type = baryonyx::fixed_array<bound>;

template<typename floatingpointT>
using c_type = baryonyx::fixed_array<floatingpointT>;
using x_type = baryonyx::fixed_array<std::int8_t>;

template<typename floatingpointT>
using pi_type = baryonyx::fixed_array<floatingpointT>;

template<typename iteratorT, typename randomT>
static void
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
static void
calculator_sort(iteratorT begin, iteratorT end, randomT& rng, minimize_tag)
{
    if (std::distance(begin, end) > 1) {

#ifdef BARYONYX_FULL_OPTIMIZATION
        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return lhs.value < rhs.value;
        });
#else
        std::stable_sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return lhs.value < rhs.value;
        });
#endif

        random_shuffle_unique(begin, end, rng);
    }
}

template<typename iteratorT, typename randomT>
static void
calculator_sort(iteratorT begin, iteratorT end, randomT& rng, maximize_tag)
{
    if (std::distance(begin, end) > 1) {

#ifdef BARYONYX_FULL_OPTIMIZATION
        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return rhs.value < lhs.value;
        });
#else
        std::stable_sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return rhs.value < lhs.value;
        });
#endif

        random_shuffle_unique(begin, end, rng);
    }
}

template<typename floatingpointT>
static inline bool
stop_iterating(floatingpointT value, minimize_tag) noexcept
{
    return value > 0;
}

template<typename floatingpointT>
static inline bool
stop_iterating(floatingpointT value, maximize_tag) noexcept
{
    return value < 0;
}

template<typename floatingpointT>
static inline bool
is_better_solution(floatingpointT lhs,
                   floatingpointT rhs,
                   minimize_tag) noexcept
{
    return lhs < rhs;
}

static inline bool
init_x(int v, minimize_tag) noexcept
{
    return v <= 0;
}

static inline bool
init_x(int v, maximize_tag) noexcept
{
    return v >= 0;
}

template<typename floatingpointT>
static inline bool
is_better_solution(floatingpointT lhs,
                   floatingpointT rhs,
                   maximize_tag) noexcept
{
    return lhs > rhs;
}

static inline bool
is_time_limit(double limit,
              std::chrono::steady_clock::time_point begin,
              std::chrono::steady_clock::time_point end) noexcept
{
    namespace sc = std::chrono;

    if (limit <= 0)
        return false;

    return sc::duration_cast<sc::duration<double>>(end - begin).count() >
           limit;
}

template<typename floatingpointT>
static inline std::size_t
compute_r_size(const AP_type<floatingpointT>& ap, int k) noexcept
{
    auto ak{ ap.row(k) };
    const auto& va{ ap.A() };

    std::size_t r_size{ 0 };

    for (auto it = std::get<0>(ak); it != std::get<1>(ak); ++it)
        if (va[it->value] != 0)
            r_size++;

    return r_size;
}

template<typename floatingpointT>
static inline std::size_t
compute_C_size(const AP_type<floatingpointT>& ap, int k) noexcept
{
    auto ak{ ap.row(k) };
    const auto& va{ ap.A() };

    std::size_t C_size{ 0 };

    for (auto it = std::get<0>(ak); it != std::get<1>(ak); ++it)
        if (va[it->value] < 0)
            C_size++;

    return C_size;
}

template<typename floatingpointT, typename modeT, typename randomT>
struct constraint_calculator
{
    using floatingpoint_type = floatingpointT;
    using mode_type = modeT;
    using random_generator_type = randomT;

    random_generator_type& rng;
    AP_type<floatingpoint_type>& ap;
    b_type& b;
    const c_type<floatingpoint_type>& cost;
    x_type& x;
    pi_type<floatingpoint_type>& pi;
    bx::fixed_array<r_data<floatingpoint_type>> r;
    bx::fixed_array<int> C; // Stores variables with negative coefficient.
    int m;
    int n;

    constraint_calculator(random_generator_type& rng_,
                          int k,
                          int m_,
                          int n_,
                          AP_type<floatingpoint_type>& ap_,
                          b_type& b_,
                          const c_type<floatingpoint_type>& c_,
                          x_type& x_,
                          pi_type<floatingpoint_type>& pi_)
      : rng(rng_)
      , ap(ap_)
      , b(b_)
      , cost(c_)
      , x(x_)
      , pi(pi_)
      , r(compute_r_size(ap_, k))
      , C(compute_C_size(ap_, k))
      , m(m_)
      , n(n_)
    {
        auto ak{ ap.row(k) };
        const auto& va{ ap.A() };

        std::size_t r_idx{ 0 }, C_idx{ 0 };
        for (auto it = std::get<0>(ak); it != std::get<1>(ak); ++it) {

            //
            // We don't need to represent the I vector since the SparseArray
            // stores no null coefficient in each row (and column).
            // I.emplace_back(ak[i].position);
            //

            if (va[it->value] != 0)
                r[r_idx++] = { 0, it->position };

            if (va[it->value] < 0)
                C[C_idx++] = it->position;
        }
    }

    bool is_valid_solution(int k) const noexcept
    {
        int v{ 0 };
        auto ak{ ap.row(k) };
        const auto& values{ ap.A() };

        for (; std::get<0>(ak) != std::get<1>(ak); ++std::get<0>(ak))
            v += values[std::get<0>(ak)->value] * x(std::get<0>(ak)->position);

        if (not(b(k).min <= v and v <= b(k).max))
            return false;

        return true;
    }

    /**
     * @brief Reinit the @c constraint_calculator.
     * @details Reinitialize the @c constaint_calculator: reset the @c r array.
     */
    void reinit()
    {
        auto first = r.begin();
        auto last = r.end();

#define BARYONYX_ASSIGN_0(p) (p)->value = 0
        BARYONYX_UNROLL_PTR(first, last, BARYONYX_ASSIGN_0);
#undef BARYONYX_ASSIGN_0
    }

    void push(floatingpoint_type objective_amplifier)
    {
#ifndef BARYONYX_FULL_OPTIMIZATION
        for (auto& elem : r)
            elem.value += objective_amplifier * cost[elem.id];
#else
        auto first = r.begin();
        auto last = r.end();

#define BARYONYX_PUSH(p) (p)->value += cost[(p)->id] * objective_amplifier
        BARYONYX_UNROLL_PTR(first, last, BARYONYX_PUSH);
#undef BARYONYX_PUSH
#endif
    }

    void update_row(int k,
                    floatingpoint_type kappa,
                    floatingpoint_type delta,
                    floatingpoint_type theta)
    {
        auto ak{ ap.row(k) };
        const auto& va{ ap.A() };

        //
        // Decrease influence of local preferences. 0 will completely reset the
        // preference values for the current row. > 0 will keep former decision
        // in mind.
        //

        ap.mult_row_p(k, theta);

        //
        // Calculate reduced costs
        //

        auto it = std::get<0>(ak);
        for (std::size_t i{ 0 }; it != std::get<1>(ak); ++it, ++i) {
            floatingpoint_type sum_a_pi{ 0 };
            floatingpoint_type sum_a_p{ 0 };

            auto H{ ap.column(it->position) };

            for (; std::get<0>(H) != std::get<1>(H); ++std::get<0>(H)) {
                sum_a_pi += ap.A(std::get<0>(H)->position, it->position) *
                            pi(std::get<0>(H)->position);
                sum_a_p += ap.A(std::get<0>(H)->position, it->position) *
                           ap.P(std::get<0>(H)->position, it->position);
                // sum_a_pi +=
                // ap.A()[std::get<0>(H)->value] *
                // pi(std::get<0>(H)->position);
                // sum_a_p += ap.A()[std::get<0>(H)->value] *
                // ap.P()[std::get<0>(H)->value];
            }

            r[i].id = it->position;
            r[i].value = cost(it->position) - sum_a_pi - sum_a_p;
        }

        //
        // Negate reduced costs and coefficients of these variables. We need to
        // parse the row Ak[i] because we need to use r[i] not available in C.
        //

        if (not C.empty()) {
            auto it = std::get<0>(ak);
            for (std::size_t i{ 0 }; it != std::get<1>(ak); ++it, ++i) {
                if (va[it->value] < 0) {
                    r[i].value = -r[i].value;
                    ap.invert_p(k, it->position);
                }
            }

            b(k).min += C.size();
            b(k).max += C.size();
        }

        calculator_sort(r.begin(), r.end(), rng, mode_type());

        //
        // The bkmin and bkmax constraint bounds are not equal and can be
        // assigned to -infinity or +infinity. We have to scan the r vector and
        // search a value j such as b(0, k) <= Sum A(k, r[j]) < b(1, k).
        //

        int i{ 0 }, selected{ -1 }, first, second;
        const auto endi = static_cast<int>(r.size());
        int sum{ 0 };

        for (; i != endi; ++i) {
            sum += 1;

            if (b(k).min <= sum)
                break;
        }

        assert(b(k).min <= sum && "b(0, k) can not be reached, this is an "
                                  "error of the preprocessing step.");

        //
        // If all variable must be assigned to 0, we let selected assigned to 0
        // and we go to the next part otherwise, we continue to scan.
        //

        if (b(k).min <= sum and sum <= b(k).max) {
            selected = i;
            for (; i != endi; ++i) {
                sum += 1;

                if (sum <= b(k).max) {
                    if (stop_iterating(r[i].value, mode_type()))
                        break;
                    ++selected;
                } else
                    break;
            }

            assert(i != endi && "unrealisable, preprocessing error");

            first = selected;
            second = selected + 1;
        }

        if (selected < 0) {
            for (int j{ 0 }; j < endi; ++j) {
                x(r[j].id) = 0;
                ap.add_p(k, r[j].id, -delta);
            }
        } else if (second >= endi) {
            for (int j{ 0 }; j < endi; ++j) {
                x(r[j].id) = 1;
                ap.add_p(k, r[j].id, +delta);
            }
        } else {
            pi(k) += ((r[first].value + r[second].value) / 2.0);

            floatingpoint_type d =
              delta +
              ((kappa / (1.0 - kappa)) * (r[second].value - r[first].value));

            int j{ 0 };
            for (; j <= selected; ++j) {
                x(r[j].id) = 1;
                ap.add_p(k, r[j].id, +d);
            }

            for (; j != endi; ++j) {
                x(r[j].id) = 0;
                ap.add_p(k, r[j].id, -d);
            }
        }

        //
        // Clean up: correct negated costs and adjust value of negated
        // variables.
        //

        if (not C.empty()) {
            b(k).min -= C.size();
            b(k).max -= C.size();

            for (std::size_t i{ 0 }, e{ C.size() }; i != e; ++i) {
                ap.invert_p(k, C[i]);
                x[C[i]] = 1 - x[C[i]];
            }
        }
    }
};

template<typename floatingpointT>
static inline c_type<floatingpointT>
make_objective_function(const bx::objective_function& obj, int n)
{
    c_type<floatingpointT> ret(n, 0);

    for (const auto& elem : obj.elements)
        ret(elem.variable_index) += elem.factor;

    return ret;
}

template<typename floatingpointT, typename modeT, typename randomT>
struct solver
{
    using floatingpoint_type = floatingpointT;
    using mode_type = modeT;
    using random_generator_type = randomT;

    random_generator_type& rng;
    std::vector<constraint_calculator<floatingpointT, modeT, randomT>>
      row_updaters;
    int m;
    int n;
    AP_type<floatingpoint_type> ap;
    b_type b;
    c_type<floatingpoint_type> c;
    x_type x;
    pi_type<floatingpoint_type> pi;

    solver(random_generator_type& rng_,
           int n_,
           c_type<floatingpoint_type> c_,
           const std::vector<bx::itm::merged_constraint>& csts)
      : rng(rng_)
      , m(csts.size())
      , n(n_)
      , ap(m, n)
      , b(m)
      , c(std::move(c_))
      , x(n)
      , pi(m)
    {
        {
            // Compute the number of elements in the matrix A then compute for
            // each rows and columns the number of elements to correctly
            // initialize the @c `matrix` structure.

            bx::fixed_array<int> r(m, 0), c(n, 0);
            int elem{ 0 };

            for (std::size_t i{ 0 }, e{ csts.size() }; i != e; ++i) {
                for (const auto& cst : csts[i].elements) {
                    r[i]++;
                    c[cst.variable_index]++;
                    ++elem;
                }
            }

            ap.reserve(elem, r.begin(), r.end(), c.begin(), c.end());
        }

        for (std::size_t i{ 0 }, e{ csts.size() }; i != e; ++i) {
            int lower{ 0 }, upper{ 0 };

            for (const auto& cst : csts[i].elements) {
                ap.set(i, cst.variable_index, cst.factor, 0.0);

                if (cst.factor > 0)
                    ++upper;
                if (cst.factor < 0)
                    --lower;
            }

            if (csts[i].min == std::numeric_limits<int>::min())
                b(i).min = lower;
            else {
                if (csts[i].min < 0)
                    b(i).min = std::max(csts[i].min, lower);
                else
                    b(i).min = csts[i].min;
            }

            if (csts[i].max == std::numeric_limits<int>::max())
                b(i).max = upper;
            else {
                if (csts[i].max > 0)
                    b(i).max = std::min(csts[i].max, upper);
                else
                    b(i).max = csts[i].max;
            }
        }

        ap.sort();

        for (int k{ 0 }, e{ m }; k != e; ++k)
            row_updaters.emplace_back(rng_, k, m, n, ap, b, c, x, pi);

        init();
    }

    void push(floatingpoint_type objective_amplifier)
    {
        for (int k{ 0 }, e{ m }; k != e; ++k)
            row_updaters[k].push(objective_amplifier);
    }

    void reinit(random_generator_type& rng_)
    {
        std::fill(ap.P().begin(), ap.P().end(), 0);
        std::fill(pi.begin(), pi.end(), 0);

        init();

        std::bernoulli_distribution d(0.5);

        for (int i = 0; i != n; ++i)
            x(i) = d(rng_);
    }

    void reinit(random_generator_type& rng_, const x_type& best_previous)
    {
        std::fill(ap.P().begin(), ap.P().end(), 0);
        std::fill(pi.begin(), pi.end(), 0);

        init();

        x = best_previous;
        std::bernoulli_distribution d(0.5);

        for (int i = 0; i != n; ++i)
            x(i) = d(rng_);
    }

    void init()
    {
        for (int i{ 0 }, e{ n }; i != e; ++i)
            x(i) = init_x(c(i), mode_type());

        for (int k{ 0 }, e{ m }; k != e; ++k)
            row_updaters[k].reinit();
    }

    void print(const std::shared_ptr<bx::context>& ctx,
               const std::vector<std::string>& names,
               int print_level) const
    {
        if (print_level <= 0)
            return;

        ctx->debug("X: ");
        for (std::size_t i{ 0 }, e{ x.size() }; i != e; ++i)
            ctx->debug("%s=%d ", names[i].c_str(), static_cast<int>(x[i]));
        ctx->debug("\n");

        for (int k{ 0 }, ek{ m }; k != ek; ++k) {
            ctx->debug("C %d:%s ",
                       k,
                       (row_updaters[k].is_valid_solution(k) ? "   valid: "
                                                             : "violated: "));
            for (std::size_t i{ 0 }, e{ row_updaters[k].r.size() }; i != e;
                 ++i)
                ctx->debug("%s [%+.6f] ",
                           names[row_updaters[k].r[i].id].c_str(),
                           (double)row_updaters[k].r[i].value);
            ctx->debug("\n");
        }
    }

    bool is_valid_solution() const noexcept
    {
        for (int k{ 0 }, ek{ m }; k != ek; ++k) {
            int v{ 0 };
            auto ak{ ap.row(k) };
            const auto& values{ ap.A() };

            for (; std::get<0>(ak) != std::get<1>(ak); ++std::get<0>(ak))
                v += values[std::get<0>(ak)->value] *
                     x(std::get<0>(ak)->position);

            if (not(b(k).min <= v and v <= b(k).max))
                return false;
        }

        return true;
    }

    bx::result results(const c_type<floatingpoint_type>& original_costs,
                       const int cost_constant) const
    {
        bx::result ret;

        if (is_valid_solution()) {
            ret.status = bx::result_status::success;
            int value = cost_constant;

            for (int i{ 0 }, ei{ n }; i != ei; ++i)
                value += original_costs[i] * x[i];

            ret.value = static_cast<double>(value);

            ret.variable_value.resize(n, 0);

            for (int i{ 0 }, ei{ n }; i != ei; ++i)
                ret.variable_value[i] = x(i);
        }

        ret.variables = n;
        ret.constraints = m;

        return ret;
    }
};

template<typename T>
static std::size_t
cycle_avoidance_hash(const std::vector<T>& vec)
{
    std::size_t seed{ vec.size() };

    for (auto& elem : vec)
        seed ^= elem + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    return seed;
}

template<typename T>
static std::size_t
cycle_avoidance_hash(const bx::fixed_array<T>& vec)
{
    std::size_t seed{ vec.size() };

    for (auto& elem : vec)
        seed ^= elem + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    return seed;
}

template<typename T>
struct no_cycle_avoidance
{
    using value_type = T;

    no_cycle_avoidance(std::size_t limit_ = 0) noexcept
    {
        (void)limit_;
    }

    bool have_cycle(const x_type&) const noexcept
    {
        return false;
    }

    bool have_cycle(const std::vector<T>&) const noexcept
    {
        return false;
    }

    bool have_cycle(const std::vector<std::pair<int, int>>&) const noexcept
    {
        return false;
    }
};

template<typename T>
struct cycle_avoidance
{
    using value_type = T;

    std::vector<std::size_t> history;
    std::size_t limit;
    int nb{ 0 };

    cycle_avoidance(std::size_t limit_ = 48l)
      : history(limit_)
      , limit(limit_)
    {
        history.clear();
    }

    ~cycle_avoidance() = default;

    bool have_cycle()
    {
        int distance{ 2 };
        const int end{ bx::numeric_cast<int>(history.size()) };

        if (end < distance)
            return false;

        for (; distance != end; ++distance) {
            int i{ end - 1 };
            int j{ i - distance };
            int cycle_size{ 0 };

            while (j >= 0 and history[i] == history[j]) {
                cycle_size++;
                --j;
                --i;
            }

            if (cycle_size > 1) {
                history.clear();
                nb++;
                return true;
            }
        }

        if (history.size() > limit)
            history.erase(history.begin(),
                          history.begin() + (history.size() - limit));

        return false;
    }

    bool have_cycle(const x_type& x)
    {
        history.push_back(cycle_avoidance_hash(x));
        return have_cycle();
    }

    bool have_cycle(const std::vector<T>& R)
    {
        history.push_back(cycle_avoidance_hash(R));
        return have_cycle();
    }

    bool have_cycle(const std::vector<std::pair<int, int>>& R)
    {
        history.push_back(cycle_avoidance_hash(R));
        return have_cycle();
    }
};

template<typename solverT>
static std::size_t
compute_missing_constraint(solverT& solver, std::vector<int>& R)
{
    R.clear();

    for (int k{ 0 }, ek{ solver.m }; k != ek; ++k) {
        int v = 0;
        auto ak{ solver.ap.row(k) };
        const auto& values{ solver.ap.A() };

        for (; std::get<0>(ak) != std::get<1>(ak); ++std::get<0>(ak))
            v += values[std::get<0>(ak)->value] *
                 solver.x(std::get<0>(ak)->position);

        if (not(solver.b(k).min <= v and v <= solver.b(k).max))
            R.push_back(k);
    }

    return R.size();
}

#ifndef BARYONYX_FULL_OPTIMIZATION
template<typename floatingpointT>
static void
print_AP(const std::shared_ptr<bx::context>& ctx,
         const AP_type<floatingpointT>& ap,
         int k,
         int rows,
         int cols)
{
    int level = ctx->get_integer_parameter("print-level", 0l);
    if (level <= 1)
        return;

    ctx->debug("P after constraint %d computation:\n", k);
    std::vector<floatingpointT> to_show(cols);

    for (int i{ 0 }; i != rows; ++i) {
        std::fill(std::begin(to_show),
                  std::end(to_show),
                  std::numeric_limits<floatingpointT>::infinity());

        auto its = ap.row(i);

        for (; std::get<0>(its) != std::get<1>(its); ++std::get<0>(its))
            to_show[std::get<0>(its)->position] =
              ap.P()[std::get<0>(its)->value];

        for (auto elem : to_show)
            if (elem == std::numeric_limits<floatingpointT>::infinity())
                ctx->debug("          ");
            else
                ctx->debug("%+.6f ", (double)elem);

        ctx->debug("\n");
    }
}
#endif

template<typename floatingpointT, typename randomT>
struct compute_reversing
{
    using random_generator_type = randomT;

    std::shared_ptr<bx::context> m_ctx;
    std::vector<int> R;
    no_cycle_avoidance<int> detect_infeasability_cycle;
    int nb = 0;

    compute_reversing(std::shared_ptr<bx::context> ctx, randomT&)
      : m_ctx(std::move(ctx))
    {
    }

    template<typename solverT>
    std::size_t run_all(solverT& solver,
                        floatingpointT kappa,
                        floatingpointT delta,
                        floatingpointT theta)
    {
        for (int i{ solver.m - 1 }; i >= 0; --i)
            solver.row_updaters[i].update_row(i, kappa, delta, theta);

        return compute_missing_constraint(solver, R);
    }

    template<typename solverT>
    std::size_t run(solverT& solver,
                    floatingpointT kappa,
                    floatingpointT delta,
                    floatingpointT theta)
    {
        auto ret = compute_missing_constraint(solver, R);
        if (ret == 0)
            return 0;

        if (detect_infeasability_cycle.have_cycle(R))
            return run_all(solver, 0.8 * kappa, delta, theta);

        for (auto it{ R.crbegin() }, et{ R.crend() }; it != et; ++it) {
            solver.row_updaters[*it].update_row(*it, kappa, delta, theta);

#ifndef BARYONYX_FULL_OPTIMIZATION
            print_AP(m_ctx, solver.ap, *it, solver.m, solver.n);
#endif
        }

        return compute_missing_constraint(solver, R);
    }
};

template<typename floatingpointT, typename randomT>
struct compute_none
{
    using random_generator_type = randomT;

    std::shared_ptr<bx::context> m_ctx;
    std::vector<int> R;
    no_cycle_avoidance<int> detect_infeasability_cycle;

    compute_none(std::shared_ptr<bx::context> ctx, randomT&)
      : m_ctx(std::move(ctx))
    {
    }

    template<typename solverT>
    std::size_t run_all(solverT& solver,
                        floatingpointT kappa,
                        floatingpointT delta,
                        floatingpointT theta)
    {
        for (int i{ 0 }; i != solver.m; ++i)
            solver.row_updaters[i].update_row(i, kappa, delta, theta);

        return compute_missing_constraint(solver, R);
    }

    template<typename solverT>
    std::size_t run(solverT& solver,
                    floatingpointT kappa,
                    floatingpointT delta,
                    floatingpointT theta)
    {
        auto ret = compute_missing_constraint(solver, R);
        if (ret == 0)
            return 0;

        if (detect_infeasability_cycle.have_cycle(R))
            return run_all(solver, 0.8 * kappa, delta, theta);

        for (auto it{ R.cbegin() }, et{ R.cend() }; it != et; ++it) {
            solver.row_updaters[*it].update_row(*it, kappa, delta, theta);
#ifndef BARYONYX_FULL_OPTIMIZATION
            print_AP(m_ctx, solver.ap, *it, solver.m, solver.n);
#endif
        }

        return compute_missing_constraint(solver, R);
    }
};

template<typename floatingpointT, typename randomT>
struct compute_random
{
    using random_generator_type = randomT;

    std::shared_ptr<bx::context> m_ctx;
    std::vector<int> R;
    no_cycle_avoidance<int> detect_infeasability_cycle;
    random_generator_type& rng;

    compute_random(std::shared_ptr<bx::context> ctx,
                   random_generator_type& rng_)
      : m_ctx(std::move(ctx))
      , rng(rng_)
    {
    }

    template<typename solverT>
    std::size_t run_all(solverT& solver,
                        floatingpointT kappa,
                        floatingpointT delta,
                        floatingpointT theta)
    {
        R.resize(solver.m);
        std::iota(R.begin(), R.end(), 0);
        std::shuffle(R.begin(), R.end(), rng);

        for (auto it{ R.cbegin() }, et{ R.cend() }; it != et; ++it)
            solver.row_updaters[*it].update_row(*it, kappa, delta, theta);

        return compute_missing_constraint(solver, R);
    }

    template<typename solverT>
    std::size_t run(solverT& solver,
                    floatingpointT kappa,
                    floatingpointT delta,
                    floatingpointT theta)
    {
        auto ret = compute_missing_constraint(solver, R);
        if (ret == 0)
            return 0;

        if (detect_infeasability_cycle.have_cycle(R))
            return run_all(solver, 0.8 * kappa, delta, theta);

        std::shuffle(R.begin(), R.end(), rng);
        for (auto it{ R.cbegin() }, et{ R.cend() }; it != et; ++it) {
            solver.row_updaters[*it].update_row(*it, kappa, delta, theta);

#ifndef BARYONYX_FULL_OPTIMIZATION
            print_AP(m_ctx, solver.ap, *it, solver.m, solver.n);
#endif
        }

        return compute_missing_constraint(solver, R);
    }
};

struct compute_infeasibility_incr
{
};

struct compute_infeasibility_decr
{
};

template<typename floatingpointT, typename randomT, typename directionT>
struct compute_infeasibility
{
    using random_generator_type = randomT;
    using direction_type = directionT;

    std::shared_ptr<bx::context> m_ctx;
    std::vector<std::pair<int, int>> R;
    no_cycle_avoidance<int> detect_infeasability_cycle;
    random_generator_type& rng;

    template<typename iteratorT>
    void sort(iteratorT begin, iteratorT end, compute_infeasibility_incr) const
    {
        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return lhs.second < rhs.second;
        });
    }

    template<typename iteratorT>
    void sort(iteratorT begin, iteratorT end, compute_infeasibility_decr) const
    {
        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return rhs.second < lhs.second;
        });
    }

    compute_infeasibility(std::shared_ptr<bx::context> ctx,
                          random_generator_type& rng_)
      : m_ctx(std::move(ctx))
      , rng(rng_)
    {
    }

    template<typename solverT>
    std::size_t compute_missing_constraint(solverT& solver)
    {
        R.clear();

        for (int k{ 0 }, ek{ solver.m }; k != ek; ++k) {
            int v = 0;

            auto ak{ solver.ap.row(k) };
            const auto& values{ solver.ap.A() };

            for (; std::get<0>(ak) != std::get<1>(ak); ++std::get<0>(ak))
                v += values[std::get<0>(ak)->value] *
                     solver.x(std::get<0>(ak)->position);

            if (solver.b(k).min > v)
                R.emplace_back(k, solver.b(k).min - v);

            if (solver.b(k).max < v)
                R.emplace_back(k, v - solver.b(k).max);
        }

        return R.size();
    }

    template<typename solverT>
    std::size_t run_all(solverT& solver,
                        floatingpointT kappa,
                        floatingpointT delta,
                        floatingpointT theta)
    {
        R.clear();

        for (int k{ 0 }, ek{ solver.m }; k != ek; ++k) {
            int v = 0;

            auto ak{ solver.ap.row(k) };
            const auto& values{ solver.ap.A() };

            for (; std::get<0>(ak) != std::get<1>(ak); ++std::get<0>(ak))
                v += values[std::get<0>(ak)->value] *
                     solver.x(std::get<1>(ak)->position);

            if (solver.b(k).min > v)
                R.emplace_back(k, solver.b(k).min - v);
            else if (solver.b(k).max < v)
                R.emplace_back(k, v - solver.b(k).max);
            else
                R.emplace_back(k, 0);
        }

        sort(R.begin(), R.end(), direction_type());

        auto ret = R.begin();
        auto it = ret + 1;
        for (; it != R.end(); ++it) {
            if (ret->second != it->second) {
                std::shuffle(ret, it, rng);
                ret = it;
            }
        }

        std::shuffle(ret, R.end(), rng);
        for (std::size_t i{ 0 }, e = { R.size() }; i != e; ++i)
            solver.row_updaters[R[i].first].update_row(
              R[i].first, kappa, delta, theta);

        return compute_missing_constraint(solver);
    }

    template<typename solverT>
    std::size_t run(solverT& solver,
                    floatingpointT kappa,
                    floatingpointT delta,
                    floatingpointT theta)
    {
        auto ret = compute_missing_constraint(solver);
        if (ret == 0)
            return 0;

        if (detect_infeasability_cycle.have_cycle(R))
            return run_all(solver, 0.8 * kappa, delta, theta);

        sort(R.begin(), R.end(), direction_type());

        auto first = R.begin();
        auto it = first + 1;
        for (; it != R.end(); ++it) {
            if (first->second != it->second) {
                std::shuffle(first, it, rng);
                first = it;
            }
        }

        std::shuffle(first, R.end(), rng);
        for (std::size_t i{ 0 }, e = { R.size() }; i != e; ++i) {
            solver.row_updaters[R[i].first].update_row(
              R[i].first, kappa, delta, theta);

#ifndef BARYONYX_FULL_OPTIMIZATION
            print_AP(m_ctx, solver.ap, R[i].first, solver.m, solver.n);
#endif
        }

        return compute_missing_constraint(solver);
    }
};

template<typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT>
struct solver_functor
{
    using floatingpoint_type = floatingpointT;
    using mode_type = modeT;
    using constraint_order_type = constraintOrderT;
    using random_generator_type = randomT;

    std::chrono::time_point<std::chrono::steady_clock> m_begin;
    std::chrono::time_point<std::chrono::steady_clock> m_end;

    std::shared_ptr<bx::context> m_ctx;
    std::vector<std::string> m_names;
    x_type m_best_x;
    bx::result m_best;

    solver_functor(std::shared_ptr<bx::context> ctx,
                   std::vector<std::string> names)
      : m_ctx(std::move(ctx))
      , m_names(std::move(names))
    {
    }

    bx::result operator()(
      const std::vector<bx::itm::merged_constraint>& constraints,
      int variables,
      const c_type<floatingpointT>& original_costs,
      const c_type<floatingpointT>& norm_costs,
      int cost_constant,
      const bx::itm::parameters& p,
      randomT& rng)
    {
        m_begin = std::chrono::steady_clock::now();
        m_end = m_begin;

        int i{ 0 };
        int i2{ 0 };
        floatingpoint_type kappa_old{ 0 };
        floatingpoint_type kappa = p.kappa_min;
        int best_remaining{ -1 };

        solver<floatingpoint_type, mode_type, random_generator_type> slv(
          rng, variables, norm_costs, constraints);

        constraint_order_type compute(m_ctx, rng);

        m_ctx->info("* solver starts:\n");

        for (;;) {
            int remaining = compute.run(slv, kappa, p.delta, p.theta);
            auto current = slv.results(original_costs, cost_constant);
            current.loop = i;
            current.remaining_constraints = remaining;
            current.duration =
              std::chrono::duration_cast<std::chrono::duration<double>>(
                m_end - m_begin)
                .count();

            if (best_remaining == -1 or remaining < best_remaining) {
                best_remaining = remaining;
                m_best = current;

                m_ctx->info(
                  "  - constraints remaining: %d/%d at %fs (loop: %d)\n",
                  remaining,
                  current.constraints,
                  current.duration,
                  i);
            }

#ifndef BARYONYX_FULL_OPTIMIZATION
            slv.print(m_ctx, m_names, p.print_level);
#endif

            if (current.status == bx::result_status::success) {
                m_ctx->info("  - Solution found: %f\n", current.value);
                m_best = current;
                return m_best;
            }

            if (i2 <= p.w) {
                kappa = p.kappa_min;
                i2++;
            } else {
                i2 = 0;
                kappa = kappa_old +
                        p.kappa_step *
                          std::pow(
                            static_cast<floatingpointT>(remaining) /
                              static_cast<floatingpointT>(current.constraints),
                            p.alpha);
                kappa_old = kappa;
            }

            if (++i > p.limit) {
                m_ctx->info("  - Loop limit reached: %d\n", i);
                m_best.status = bx::result_status::limit_reached;
                return m_best;
            }

            if (kappa > p.kappa_max) {
                m_ctx->info("  - Kappa max reached: %+.6f\n", (double)kappa);
                m_best.status = bx::result_status::kappa_max_reached;
                return m_best;
            }

            m_end = std::chrono::steady_clock::now();
            if (is_time_limit(p.time_limit, m_begin, m_end)) {
                m_ctx->info(
                  "  - Time limit reached: %d %+.6f\n", i, (double)kappa);
                m_best.status = bx::result_status::time_limit_reached;
                return m_best;
            }
        }
    }
};

template<typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT>
struct optimize_functor
{
    using floatingpoint_type = floatingpointT;
    using mode_type = modeT;
    using constraint_order_type = constraintOrderT;
    using random_generator_type = randomT;

    std::chrono::time_point<std::chrono::steady_clock> m_begin;
    std::chrono::time_point<std::chrono::steady_clock> m_end;

    std::shared_ptr<bx::context> m_ctx;
    int m_thread_id;
    const std::vector<std::string>& m_variable_names;
    const bx::affected_variables& m_affected_vars;
    x_type m_best_x;
    bx::result m_best;

    optimize_functor(std::shared_ptr<bx::context> ctx,
                     int thread_id,
                     const std::vector<std::string>& variable_names,
                     const bx::affected_variables& affected_vars)
      : m_ctx(std::move(ctx))
      , m_thread_id(thread_id)
      , m_variable_names(variable_names)
      , m_affected_vars(affected_vars)
    {
    }

    bx::result operator()(
      const std::vector<bx::itm::merged_constraint>& constraints,
      int variables,
      const c_type<floatingpointT>& original_costs,
      const c_type<floatingpointT>& norm_costs,
      int cost_constant,
      const bx::itm::parameters& p,
      randomT& rng)
    {
        m_begin = std::chrono::steady_clock::now();
        m_end = m_begin;

        int i{ 0 };
        int i2{ 0 };
        floatingpoint_type kappa_old{ 0 };
        floatingpoint_type kappa = p.kappa_min;

        int pushed{ -1 };
        int pushing_iteration{ 0 };

        solver<floatingpoint_type, mode_type, random_generator_type> slv(
          rng, variables, norm_costs, constraints);

        constraint_order_type compute(m_ctx, rng);

        for (; not is_time_limit(p.time_limit, m_begin, m_end);
             m_end = std::chrono::steady_clock::now(), ++i) {

            int remaining = compute.run(slv, kappa, p.delta, p.theta);

            auto current = slv.results(original_costs, cost_constant);
            current.loop = i;
            current.remaining_constraints = remaining;

            if (store_if_better(current)) {
                m_best_x = slv.x;
                pushed = 0;
            }

            if (i2 <= p.w) {
                kappa = p.kappa_min;
                i2++;
            } else {
                i2 = 0;
                kappa = kappa_old +
                        p.kappa_step *
                          std::pow(static_cast<double>(remaining) /
                                     static_cast<double>(current.constraints),
                                   p.alpha);

                kappa_old = kappa;
            }

            if (i >= p.limit or kappa > p.kappa_max or
                pushed > p.pushes_limit) {
                if (m_best.status == bx::result_status::success) {
                    slv.reinit(rng, m_best_x);
                } else {
                    slv.reinit(rng);
                }

                i = 0;
                i2 = 0;
                kappa_old = 0.0;
                kappa = p.kappa_min;
                pushed = -1;
                pushing_iteration = 0;

                continue;
            }

            if (pushed >= 0) {
                ++pushing_iteration;

                if (pushing_iteration >= p.pushing_iteration_limit) {
                    pushed++;
                    pushing_iteration = 0;

                    slv.push(p.pushing_objective_amplifier);

                    remaining = compute.run_all(
                      slv, p.pushing_k_factor * kappa, p.delta, p.theta);

                    if (remaining == 0) {
                        current = slv.results(original_costs, cost_constant);
                        current.loop = i;
                        if (store_if_better(current)) {
                            lp_debug(m_ctx, "  `-> push found new solution\n");
                            m_best_x = slv.x;
                        }
                    }
                }
            }
        }

        return m_best;
    }

private:
    bool store_if_better(const bx::result& current) noexcept
    {
        if (current.status != bx::result_status::success)
            return false;

        if (m_best.status != bx::result_status::success or
            is_better_solution(current.value, m_best.value, mode_type())) {

            double t =
              std::chrono::duration_cast<std::chrono::duration<double>>(
                m_end - m_begin)
                .count();

            m_ctx->info("  - Solution found: %f (i=%d t=%fs thread:%d)\n",
                        current.value,
                        current.loop,
                        t,
                        m_thread_id);

            m_best = current;
            m_best.duration = t;

            std::ofstream ofs(bx::stringf("temp-%d.sol", m_thread_id));
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
static void
random_epsilon_unique(iteratorT begin,
                      iteratorT end,
                      randomT& rng,
                      floatingpointT min,
                      floatingpointT max)
{
    assert(min != max && "rng_normalize_cost fail to define min and max");

    std::uniform_real_distribution<floatingpointT> distribution(min, max);

    for (; begin != end; ++begin)
        begin->second += distribution(rng);
}

template<typename floatingpointT, typename randomT>
static c_type<floatingpointT>
rng_normalize_costs(const c_type<floatingpointT>& c, randomT& rng)
{
    std::vector<std::pair<floatingpointT, int>> r(c.size());
    int i, e;

    for (i = 0, e = bx::numeric_cast<int>(c.size()); i != e; ++i) {
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

    c_type<floatingpointT> ret(c);
    for (i = 0, e = bx::numeric_cast<int>(c.size()); i != e; ++i)
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
 *      input vector is too small or with infinity value, the c is unchanged.
 */
template<typename floatingpointT, typename randomT>
static c_type<floatingpointT>
normalize_costs(const std::shared_ptr<bx::context>& ctx,
                const std::string& norm,
                const c_type<floatingpointT>& c,
                randomT& rng)
{
    if (norm == "none") {
        ctx->info("  - No norm");
        return c;
    }

    if (norm == "rng") {
        ctx->info("  - Compute random norm\n");
        return rng_normalize_costs(c, rng);
    }

    c_type<floatingpointT> ret(c);
    double div{ 0 };

    if (norm == "l1") {
        ctx->info("  - Compute l1 norm\n");
        for (auto elem : ret)
            div += std::abs(elem);
    } else if (norm == "l2") {
        ctx->info("  - Compute l2 norm\n");
        for (auto elem : ret)
            div += elem * elem;
    } else {
        ctx->info("  - Compute infinity-norm (default)\n");
        div = *std::max_element(c.cbegin(), c.cend());
    }

    if (std::isnormal(div))
        for (auto& elem : ret)
            elem /= div;

    return ret;
}

template<typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT>
static bx::result
solve(std::shared_ptr<bx::context> ctx,
      bx::problem& pb,
      const bx::itm::parameters& p,
      randomT& rng)
{
    ctx->info("Solver initializing\n");

    bx::result ret;

    auto constraints{ bx::itm::make_merged_constraints(ctx, pb, p) };
    if (not constraints.empty() and not pb.vars.values.empty()) {
        auto variables = bx::numeric_cast<int>(pb.vars.values.size());
        auto cost =
          make_objective_function<floatingpointT>(pb.objective, variables);
        auto norm_costs = normalize_costs(ctx, p.norm, cost, rng);
        auto cost_constant = pb.objective.value;
        auto names = std::move(pb.vars.names);
        auto affected_vars = std::move(pb.affected_vars);

        bx::clear(pb);

        solver_functor<floatingpointT, modeT, constraintOrderT, randomT> slv(
          ctx, names);

        ret =
          slv(constraints, variables, cost, norm_costs, cost_constant, p, rng);

        ret.method = "inequalities_1coeff solver";
        ret.variable_name = std::move(names);
        ret.affected_vars = std::move(affected_vars);
    } else {
        ret.status = bx::result_status::success;
    }

    return ret;
}

template<typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT>
static bx::result
optimize(std::shared_ptr<bx::context> ctx,
         bx::problem& pb,
         const bx::itm::parameters& p,
         randomT& rng,
         int thread)
{
    bx::Expects(thread >= 1, "optimize: bad thread number");

    ctx->info("Optimizer initializing\n");

    bx::result ret;

    auto constraints{ bx::itm::make_merged_constraints(ctx, pb, p) };
    if (not constraints.empty() and not pb.vars.values.empty()) {
        auto variables = bx::numeric_cast<int>(pb.vars.values.size());
        auto cost =
          make_objective_function<floatingpointT>(pb.objective, variables);
        auto norm_costs = normalize_costs(ctx, p.norm, cost, rng);
        auto cost_constant = pb.objective.value;
        auto names = std::move(pb.vars.names);
        auto affected_vars = std::move(pb.affected_vars);

        bx::clear(pb);

        std::vector<std::thread> pool(thread);
        pool.clear();
        std::vector<std::future<bx::result>> results(thread);
        results.clear();

        if (thread == 1)
            ctx->info("optimizer starts with one thread\n");
        else
            ctx->info("Optimizer starts with %d threads\n", thread);

        for (int i{ 0 }; i != thread; ++i) {
            std::packaged_task<bx::result()> task(std::bind(
              optimize_functor<floatingpointT,
                               modeT,
                               constraintOrderT,
                               randomT>(ctx, i, names, affected_vars),
              std::ref(constraints),
              variables,
              std::ref(cost),
              std::ref(norm_costs),
              cost_constant,
              std::ref(p),
              std::ref(rng)));

            results.emplace_back(task.get_future());

            pool.emplace_back(std::thread(std::move(task)));
        }

        for (auto& t : pool)
            t.join();

        ret = results[0].get();
        for (int i{ 1 }; i != thread; ++i) {
            auto current = results[i].get();
            if (current.status == bx::result_status::success) {
                if (ret.status != bx::result_status::success or
                    is_better_solution(current.value, ret.value, modeT()))
                    ret = current;
            }
        }

        ret.method = "inequalities_1coeff optimizer";
        ret.variable_name = std::move(names);
        ret.affected_vars = std::move(affected_vars);
    } else {
        ret.status = bx::result_status::success;
    }

    return ret;
}

template<typename realT, typename modeT, typename randomT>
static bx::result
dispatch_solve(std::shared_ptr<bx::context> ctx,
               bx::problem& pb,
               const bx::itm::parameters& p,
               randomT& rng)
{

    switch (p.order) {
    case bx::itm::constraint_order::none:
        return ::solve<realT, modeT, ::compute_none<realT, randomT>, randomT>(
          ctx, pb, p, rng);
    case bx::itm::constraint_order::reversing:
        return ::
          solve<realT, modeT, ::compute_reversing<realT, randomT>, randomT>(
            ctx, pb, p, rng);
    case bx::itm::constraint_order::random_sorting:
        return ::solve<realT, modeT, ::compute_random<realT, randomT>>(
          ctx, pb, p, rng);
    case bx::itm::constraint_order::infeasibility_decr:
        return ::solve<
          realT,
          modeT,
          ::compute_infeasibility<realT, randomT, compute_infeasibility_decr>>(
          ctx, pb, p, rng);
    case bx::itm::constraint_order::infeasibility_incr:
        return ::solve<
          realT,
          modeT,
          ::compute_infeasibility<realT, randomT, compute_infeasibility_incr>>(
          ctx, pb, p, rng);
    }

    return {};
}

template<typename realT, typename modeT, typename randomT>
static bx::result
dispatch_optimize(std::shared_ptr<bx::context> ctx,
                  bx::problem& pb,
                  const bx::itm::parameters& p,
                  randomT& rng,
                  int thread)
{

    switch (p.order) {
    case bx::itm::constraint_order::none:
        return ::
          optimize<realT, modeT, ::compute_none<realT, randomT>, randomT>(
            ctx, pb, p, rng, thread);
    case bx::itm::constraint_order::reversing:
        return ::
          optimize<realT, modeT, ::compute_reversing<realT, randomT>, randomT>(
            ctx, pb, p, rng, thread);
    case bx::itm::constraint_order::random_sorting:
        return ::optimize<realT, modeT, ::compute_random<realT, randomT>>(
          ctx, pb, p, rng, thread);
    case bx::itm::constraint_order::infeasibility_decr:
        return ::optimize<
          realT,
          modeT,
          ::compute_infeasibility<realT, randomT, compute_infeasibility_decr>>(
          ctx, pb, p, rng, thread);
    case bx::itm::constraint_order::infeasibility_incr:
        return ::optimize<
          realT,
          modeT,
          ::compute_infeasibility<realT, randomT, compute_infeasibility_incr>>(
          ctx, pb, p, rng, thread);
    }

    return {};
}

} // anonymous namespace

namespace baryonyx {
namespace itm {

result
inequalities_1coeff_wedelin_solve(const std::shared_ptr<baryonyx::context>& ctx,
                                  problem& pb)
{
    parameters p(ctx);

    using random_generator_type = std::default_random_engine;
    random_generator_type::result_type seed = ctx->get_integer_parameter(
      "seed", std::chrono::system_clock::now().time_since_epoch().count());
    random_generator_type rng(seed);

    if (pb.type == baryonyx::objective_function_type::maximize) {
        switch (p.float_type) {
        case floating_point_type::float_type:
            return dispatch_solve<float,
                                  ::maximize_tag,
                                  random_generator_type>(ctx, pb, p, rng);
        case floating_point_type::double_type:
            return dispatch_solve<double,
                                  ::maximize_tag,
                                  random_generator_type>(ctx, pb, p, rng);
        case floating_point_type::longdouble_type:
            return dispatch_solve<long double,
                                  ::maximize_tag,
                                  random_generator_type>(ctx, pb, p, rng);
        }
    } else {
        switch (p.float_type) {
        case floating_point_type::float_type:
            return dispatch_solve<float,
                                  ::minimize_tag,
                                  random_generator_type>(ctx, pb, p, rng);
        case floating_point_type::double_type:
            return dispatch_solve<double,
                                  ::minimize_tag,
                                  random_generator_type>(ctx, pb, p, rng);
        case floating_point_type::longdouble_type:
            return dispatch_solve<long double,
                                  ::minimize_tag,
                                  random_generator_type>(ctx, pb, p, rng);
        }
    }

    return {};
}

result
inequalities_1coeff_wedelin_optimize(const std::shared_ptr<baryonyx::context>& ctx,
                                     problem& pb,
                                     int thread)
{
    parameters p(ctx);

    using random_generator_type = std::default_random_engine;
    random_generator_type::result_type seed = ctx->get_integer_parameter(
      "seed", std::chrono::system_clock::now().time_since_epoch().count());

    random_generator_type rng(seed);

    if (pb.type == baryonyx::objective_function_type::maximize) {
        switch (p.float_type) {
        case floating_point_type::float_type:
            return dispatch_optimize<float,
                                     ::maximize_tag,
                                     random_generator_type>(
              ctx, pb, p, rng, thread);
        case floating_point_type::double_type:
            return dispatch_optimize<double,
                                     ::maximize_tag,
                                     random_generator_type>(
              ctx, pb, p, rng, thread);
        case floating_point_type::longdouble_type:
            return dispatch_optimize<long double,
                                     ::maximize_tag,
                                     random_generator_type>(
              ctx, pb, p, rng, thread);
        }
    } else {
        switch (p.float_type) {
        case floating_point_type::float_type:
            return dispatch_optimize<float,
                                     ::minimize_tag,
                                     random_generator_type>(
              ctx, pb, p, rng, thread);
        case floating_point_type::double_type:
            return dispatch_optimize<double,
                                     ::minimize_tag,
                                     random_generator_type>(
              ctx, pb, p, rng, thread);
        case floating_point_type::longdouble_type:
            return dispatch_optimize<long double,
                                     ::minimize_tag,
                                     random_generator_type>(
              ctx, pb, p, rng, thread);
        }
    }

    return {};
}
}
}
