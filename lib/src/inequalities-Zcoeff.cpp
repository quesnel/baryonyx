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
#include <utility>
#include <vector>

#include "branch-and-bound-solver.hpp"
#include "fixed_array.hpp"
#include "itm.hpp"
#include "knapsack-dp-solver.hpp"
#include "matrix.hpp"
#include "private.hpp"
#include "utils.hpp"

#include <cassert>

namespace bx = baryonyx;

namespace {

using bx::length;

struct maximize_tag
{};

struct minimize_tag
{};

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

template<typename floatingpointT>
using AP_type = bx::SparseArray<int, floatingpointT>;

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
calculator_sort(iteratorT begin, iteratorT end, randomT& rng, bx::minimize_tag)
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
calculator_sort(iteratorT begin, iteratorT end, randomT& rng, bx::maximize_tag)
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

template<typename floatingpointT, typename randomT>
static inline bool
stop_iterating(floatingpointT value, randomT& rng, bx::minimize_tag) noexcept
{
    if (value == 0) {
        std::bernoulli_distribution d(0.5);
        return d(rng);
    }

    return value > 0;
}

template<typename floatingpointT, typename randomT>
static inline bool
stop_iterating(floatingpointT value, randomT& rng, bx::maximize_tag) noexcept
{
    if (value == 0) {
        std::bernoulli_distribution d(0.5);
        return d(rng);
    }

    return value < 0;
}

template<typename floatingpointT>
static inline bool
is_better_solution(floatingpointT lhs,
                   floatingpointT rhs,
                   bx::minimize_tag) noexcept
{
    return lhs < rhs;
}

template<typename floatingpointT>
static inline bool
is_better_solution(floatingpointT lhs,
                   floatingpointT rhs,
                   bx::maximize_tag) noexcept
{
    return lhs > rhs;
}

static inline bool
init_x(int v, bx::minimize_tag) noexcept
{
    return v <= 0;
}

static inline bool
init_x(int v, bx::maximize_tag) noexcept
{
    return v >= 0;
}

template<typename apT, typename xT, typename bT>
bool
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
int
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
void
print_missing_constraint(const std::shared_ptr<baryonyx::context>& ctx,
                         const apT& ap,
                         const xT& x,
                         const bT& b,
                         const std::vector<std::string>& names) noexcept
{
    std::vector<int> R;

    compute_missing_constraint(ap, x, b, R);
    ctx->info("Constraints remaining %d:\n", length(R));

    typename apT::const_iterator it, et;
    const auto& va = ap.A();

    for (auto k : R) {
        std::tie(it, et) = ap.row(k);
        int v = 0;

        ctx->info("%d: %d <= ", k, b(k).min);

        for (; it != et; ++it) {
            v += va[it->value] * x[it->position];

            ctx->info("%+d [%s (%d)] ",
                      va[it->value],
                      names[it->position].c_str(),
                      x[it->position]);
        }

        ctx->info(" <= %d | value: %d\n", b(k).max, v);
    }
}

template<typename floatingpointT, typename modeT, typename randomT>
struct solver
{
    using floatingpoint_type = floatingpointT;
    using mode_type = modeT;
    using random_type = randomT;

    random_type& rng;

    // Sparse matrix to store A and P values.
    AP_type<floatingpoint_type> ap;

    // Vector shared between all constraints to store the reduced cost.
    bx::fixed_array<r_data<floatingpoint_type>> R;

    // Vector for each constraint with negative coefficients.
    bx::fixed_array<bx::fixed_array<c_data>> C;

    // Vector of boolean where true informs a Z coefficient in the equation or
    // inequation.
    std::vector<bool> Z;

    // Bound vector.
    b_type b;
    const c_type<floatingpoint_type>& c;
    x_type x;
    pi_type<floatingpoint_type> pi;
    int m;
    int n;

    solver(random_type& rng_,
           int n_,
           const c_type<floatingpoint_type>& c_,
           const std::vector<bx::itm::merged_constraint>& csts)
      : rng(rng_)
      , ap(length(csts), n_)
      , C(length(csts))
      , Z(length(csts), false)
      , b(length(csts))
      , c(c_)
      , x(n_)
      , pi(length(csts))
      , m(length(csts))
      , n(n_)
    {
        {
            // Compute the number of elements in the matrix A then compute for
            // each rows and columns the number of elements to correctly
            // initialize the @c `matrix` structure.

            bx::fixed_array<int> rinit(m, 0), cinit(n, 0);
            int elem{ 0 };

            for (int i{ 0 }, e{ length(csts) }; i != e; ++i) {
                for (const auto& cst : csts[i].elements) {
                    rinit[i]++;
                    cinit[cst.variable_index]++;
                    ++elem;
                }
            }

            ap.reserve(
              elem, rinit.begin(), rinit.end(), cinit.begin(), cinit.end());
        }

        {
            // Compute the minimal bounds for each constraints, default
            // constraints are -oo <= ... <= bkmax, bkmin <= ... <= +oo and
            // bkmin <= ... <= bkmax. This code remove infinity and replace
            // with minimal or maximal value of the constraint.

            for (int i = 0, e = length(csts); i != e; ++i) {
                int lower = 0, upper = 0;

                for (const auto& cst : csts[i].elements) {
                    ap.set(i, cst.variable_index, cst.factor, 0.0);

                    if (cst.factor > 0)
                        upper += cst.factor;

                    if (cst.factor < 0)
                        lower += cst.factor;

                    if (std::abs(cst.factor) > 1)
                        Z[i] = true;
                }

                if (csts[i].min == csts[i].max) {
                    b(i).min = csts[i].min;
                    b(i).max = csts[i].max;
                } else {
                    if (csts[i].min == std::numeric_limits<int>::min()) {
                        b(i).min = lower;
                    } else {
                        if (lower < 0)
                            b(i).min = std::max(lower, csts[i].min);
                        else
                            b(i).min = csts[i].min;
                    }

                    if (csts[i].max == std::numeric_limits<int>::max()) {
                        b(i).max = upper;
                    } else {
                        b(i).max = csts[i].max;
                    }
                }
            }

            ap.sort();
        }

        {
            //
            // Compute the R vector size and the C vectors for each constraints
            // with negative coefficient.
            //

            int rsizemax = 0;
            for (int i = 0, e = length(csts); i != e; ++i) {
                int rsize = 0, csize = 0;

                for (const auto& cst : csts[i].elements) {
                    if (cst.factor < 0)
                        ++csize;
                    ++rsize;
                }

                rsizemax = std::max(rsizemax, rsize);

                if (csize <= 0) // No negative coefficient found in
                    continue;   // constraint, try the next.

                C[i] = bx::fixed_array<c_data>(csize);

                int id_in_r = 0;
                int id_in_c = 0;

                typename AP_type<floatingpoint_type>::const_iterator it, et;
                std::tie(it, et) = ap.row(i);

                for (; it != et; ++it) {
                    if (ap.A()[it->value] < 0) {
                        C[i][id_in_c].id_r = id_in_r;
                        C[i][id_in_c].id_A = it->position;
                        ++id_in_c;
                    }
                    ++id_in_r;
                }
            }

            R = bx::fixed_array<r_data<floatingpoint_type>>(rsizemax);
        }

        std::fill(ap.P().begin(), ap.P().end(), 0);
        std::fill(pi.begin(), pi.end(), 0);

        for (int i = 0, e = n; i != e; ++i)
            x(i) = init_x(c(i), mode_type());
    }

    void reinit(const x_type& best_previous, double reverse_solution)
    {
        std::fill(ap.P().begin(), ap.P().end(), 0);
        std::fill(pi.begin(), pi.end(), 0);

        if (not best_previous.empty()) {
            x = best_previous;
            std::bernoulli_distribution d(reverse_solution);

            for (int i = 0; i != n; ++i)
                if (d(rng))
                    x(i) = !x(i);
        } else {
            for (int i = 0, e = n; i != e; ++i)
                x(i) = init_x(c(i), mode_type());
        }
    }

    void print(const std::shared_ptr<bx::context>& ctx,
               const std::vector<std::string>& names,
               int print_level) const
    {
        if (print_level <= 0)
            return;

        ctx->debug("X: ");
        for (int i = 0, e = length(x); i != e; ++i)
            ctx->debug("%s=%d ", names[i].c_str(), static_cast<int>(x[i]));
        ctx->debug("\n");

        for (int k = 0, ek = m; k != ek; ++k) {
            auto ak = ap.row(k);
            int v = 0;

            for (; std::get<0>(ak) != std::get<1>(ak); ++std::get<0>(ak))
                v += ap.A()[std::get<0>(ak)->value] *
                     x(std::get<0>(ak)->position);

            bool valid = b(k).min <= v and v <= b(k).max;
            ctx->debug("C %d:%s\n", k, (valid ? "   valid" : "violated"));
        }
    }

    bx::result results(const c_type<floatingpoint_type>& original_costs,
                       const double cost_constant) const
    {
        bx::result ret;

        if (is_valid_solution(ap, x, b)) {
            ret.status = bx::result_status::success;
            double value = cost_constant;

            for (int i{ 0 }, ei{ n }; i != ei; ++i)
                value += original_costs[i] * x[i];

            ret.value = static_cast<double>(value);
        }

        ret.variable_value.resize(n, 0);

        for (int i{ 0 }, ei{ n }; i != ei; ++i)
            ret.variable_value[i] = x(i);

        ret.variables = n;
        ret.constraints = m;

        return ret;
    }

    void compute_update_row_Z_eq(int k,
                                 int bk,
                                 floatingpoint_type kappa,
                                 floatingpoint_type delta,
                                 floatingpoint_type theta,
                                 floatingpoint_type objective_amplifier)
    {
        typename AP_type<floatingpoint_type>::const_iterator it, et;
        std::tie(it, et) = ap.row(k);

        decrease_preference(it, et, theta);

        const auto& ck = C[k];
        const int c_size = length(ck);
        const int r_size = compute_reduced_costs(it, et);
        int bk_move = 0;

        //
        // Before sort and select variables, we apply the push method: for each
        // reduces cost, we had the cost multiply with an objective amplifier.
        //

        if (objective_amplifier)
            for (int i = 0; i != r_size; ++i)
                R[i].value += objective_amplifier * c[R[i].id];

        //
        // Negate reduced costs and coefficients of these variables. We need to
        // parse the row Ak[i] because we need to use r[i] not available in C.
        //

        for (int i = 0; i != c_size; ++i) {
            R[ck[i].id_r].value = -R[ck[i].id_r].value;
            ap.invert_p(k, ck[i].id_A);
            bk_move += ap.A()[ck[i].id_A];
        }

        bk += std::abs(bk_move);

        int selected =
          baryonyx::branch_and_bound_solver<modeT, floatingpoint_type>(
            ap, R, it, et, bk);

        affect_variables(k, selected, r_size, kappa, delta);

        //
        // Clean up: correct negated costs and adjust value of negated
        // variables.
        //

        for (int i = 0; i != c_size; ++i) {
            ap.invert_p(k, ck[i].id_A);
            x[ck[i].id_A] = 1 - x[ck[i].id_A];
        }
    }

    void compute_update_row_Z_ineq(int k,
                                   int bkmin,
                                   int bkmax,
                                   floatingpoint_type kappa,
                                   floatingpoint_type delta,
                                   floatingpoint_type theta,
                                   floatingpoint_type objective_amplifier)
    {
        typename AP_type<floatingpoint_type>::const_iterator it, et;
        std::tie(it, et) = ap.row(k);

        decrease_preference(it, et, theta);

        const auto& ck = C[k];
        const int c_size = (ck) ? length(ck) : 0;
        const int r_size = compute_reduced_costs(it, et);
        int bk_move = 0;

        //
        // Before sort and select variables, we apply the push method: for each
        // reduces cost, we had the cost multiply with an objective amplifier.
        //

        if (objective_amplifier)
            for (int i = 0; i != r_size; ++i)
                R[i].value += objective_amplifier * c[R[i].id];

        //
        // Negate reduced costs and coefficients of these variables. We need to
        // parse the row Ak[i] because we need to use r[i] not available in C.
        //

        for (int i = 0; i != c_size; ++i) {
            R[ck[i].id_r].value = -R[ck[i].id_r].value;
            ap.invert_p(k, ck[i].id_A);
            bk_move += ap.A()[ck[i].id_A];
        }

        bkmin += std::abs(bk_move);
        bkmax += std::abs(bk_move);

        int selected =
          baryonyx::branch_and_bound_solver<modeT, floatingpoint_type>(
            ap, R, it, et, bkmax);

        affect_variables(k, selected, r_size, kappa, delta);

        //
        // Clean up: correct negated costs and adjust value of negated
        // variables.
        //

        for (int i = 0; i != c_size; ++i) {
            ap.invert_p(k, ck[i].id_A);
            x[ck[i].id_A] = 1 - x[ck[i].id_A];
        }
    }

    void compute_update_row_01_eq(int k,
                                  int bk,
                                  floatingpoint_type kappa,
                                  floatingpoint_type delta,
                                  floatingpoint_type theta,
                                  floatingpoint_type objective_amplifier)
    {
        typename AP_type<floatingpoint_type>::const_iterator it, et;
        std::tie(it, et) = ap.row(k);

        decrease_preference(it, et, theta);

        const int r_size = compute_reduced_costs(it, et);

        //
        // Before sort and select variables, we apply the push method: for each
        // reduces cost, we had the cost multiply with an objective amplifier.
        //

        if (objective_amplifier)
            for (int i = 0; i != r_size; ++i)
                R[i].value += objective_amplifier * c[R[i].id];

        calculator_sort(R.begin(), R.begin() + r_size, rng, mode_type());

        int selected = select_variables_equality(r_size, bk);

        affect_variables(k, selected, r_size, kappa, delta);
    }

    void compute_update_row_01_ineq(int k,
                                    int bkmin,
                                    int bkmax,
                                    floatingpoint_type kappa,
                                    floatingpoint_type delta,
                                    floatingpoint_type theta,
                                    floatingpoint_type objective_amplifier)
    {
        typename AP_type<floatingpoint_type>::const_iterator it, et;
        std::tie(it, et) = ap.row(k);

        decrease_preference(it, et, theta);

        const int r_size = compute_reduced_costs(it, et);

        //
        // Before sort and select variables, we apply the push method: for each
        // reduces cost, we had the cost multiply with an objective amplifier.
        //

        if (objective_amplifier)
            for (int i = 0; i != r_size; ++i)
                R[i].value += objective_amplifier * c[R[i].id];

        calculator_sort(R.begin(), R.begin() + r_size, rng, mode_type());

        int selected = select_variables_inequality(r_size, bkmin, bkmax);

        affect_variables(k, selected, r_size, kappa, delta);
    }

    void compute_update_row_101_eq(int k,
                                   int bk,
                                   floatingpoint_type kappa,
                                   floatingpoint_type delta,
                                   floatingpoint_type theta,
                                   floatingpoint_type objective_amplifier)
    {
        typename AP_type<floatingpoint_type>::const_iterator it, et;
        std::tie(it, et) = ap.row(k);

        decrease_preference(it, et, theta);

        const auto& ck = C[k];
        const int c_size = length(ck);
        const int r_size = compute_reduced_costs(it, et);

        //
        // Before sort and select variables, we apply the push method: for each
        // reduces cost, we had the cost multiply with an objective amplifier.
        //

        if (objective_amplifier)
            for (int i = 0; i != r_size; ++i)
                R[i].value += objective_amplifier * c[R[i].id];

        //
        // Negate reduced costs and coefficients of these variables. We need to
        // parse the row Ak[i] because we need to use r[i] not available in C.
        //

        for (int i = 0; i != c_size; ++i) {
            R[ck[i].id_r].value = -R[ck[i].id_r].value;
            ap.invert_p(k, ck[i].id_A);
        }

        bk += c_size;

        calculator_sort(R.begin(), R.begin() + r_size, rng, mode_type());

        int selected = select_variables_equality(r_size, bk);

        affect_variables(k, selected, r_size, kappa, delta);

        //
        // Clean up: correct negated costs and adjust value of negated
        // variables.
        //

        for (int i = 0; i != c_size; ++i) {
            ap.invert_p(k, ck[i].id_A);
            x[ck[i].id_A] = 1 - x[ck[i].id_A];
        }
    }

    void compute_update_row_101_ineq(int k,
                                     int bkmin,
                                     int bkmax,
                                     floatingpoint_type kappa,
                                     floatingpoint_type delta,
                                     floatingpoint_type theta,
                                     floatingpoint_type objective_amplifier)
    {
        typename AP_type<floatingpoint_type>::const_iterator it, et;
        std::tie(it, et) = ap.row(k);

        decrease_preference(it, et, theta);

        const auto& ck = C[k];
        const int c_size = (ck) ? length(ck) : 0;
        const int r_size = compute_reduced_costs(it, et);

        //
        // Before sort and select variables, we apply the push method: for each
        // reduces cost, we had the cost multiply with an objective amplifier.
        //

        if (objective_amplifier)
            for (int i = 0; i != r_size; ++i)
                R[i].value += objective_amplifier * c[R[i].id];

        //
        // Negate reduced costs and coefficients of these variables. We need to
        // parse the row Ak[i] because we need to use r[i] not available in C.
        //

        for (int i = 0; i != c_size; ++i) {
            R[ck[i].id_r].value = -R[ck[i].id_r].value;
            ap.invert_p(k, ck[i].id_A);
        }

        bkmin += c_size;
        bkmax += c_size;

        calculator_sort(R.begin(), R.begin() + r_size, rng, mode_type());

        int selected = select_variables_inequality(r_size, bkmin, bkmax);

        affect_variables(k, selected, r_size, kappa, delta);

        //
        // Clean up: correct negated costs and adjust value of negated
        // variables.
        //

        for (int i = 0; i != c_size; ++i) {
            ap.invert_p(k, ck[i].id_A);
            x[ck[i].id_A] = 1 - x[ck[i].id_A];
        }
    }

    //
    // Decrease influence of local preferences. 0 will completely reset the
    // preference values for the current row. > 0 will keep former decision in
    // mind.
    //
    template<typename iteratorT>
    void decrease_preference(iteratorT begin,
                             iteratorT end,
                             floatingpoint_type theta) noexcept
    {
        for (; begin != end; ++begin)
            ap.P()[begin->value] *= theta;
    }

    //
    // Compute the reduced costs and return the size of the newly R vector.
    //
    template<typename iteratorT>
    int compute_reduced_costs(iteratorT begin, iteratorT end) noexcept
    {
        int r_size = 0;

        for (; begin != end; ++begin) {
            floatingpoint_type sum_a_pi = 0;
            floatingpoint_type sum_a_p = 0;

            typename AP_type<floatingpoint_type>::const_iterator ht, hend;
            std::tie(ht, hend) = ap.column(begin->position);

            for (; ht != hend; ++ht) {
                auto a = ap.A()[ht->value];
                sum_a_pi += a * pi[ht->position];
                sum_a_p += a * ap.P()[ht->value];
            }

            R[r_size].id = begin->position;
            R[r_size].value = c[begin->position] - sum_a_pi - sum_a_p;
            ++r_size;
        }

        return r_size;
    }

    int select_variables_equality(const int r_size, int bk)
    {
        (void)r_size;

        assert(bk <= r_size && "b(k) can not be reached, this is an "
                               "error of the preprocessing step.");

        return bk - 1;
    }

    int select_variables_inequality(const int r_size, int bkmin, int bkmax)
    {
        int i = 0;
        int selected = -1;
        int sum = 0;

        for (; i != r_size; ++i) {
            sum += 1;

            if (bkmin <= sum)
                break;
        }

        assert(bkmin <= sum && "b(0, k) can not be reached, this is an "
                               "error of the preprocessing step.");

        if (bkmin <= sum and sum <= bkmax) {
            selected = i;
            for (; i != r_size; ++i) {
                sum += 1;

                if (sum <= bkmax) {
                    if (stop_iterating(R[i].value, rng, mode_type()))
                        break;
                    ++selected;
                } else
                    break;
            }

            assert(i != r_size && "unrealizable, preprocessing error");
        }

        return selected;
    }

    //
    // The bkmin and bkmax constraint bounds are not equal and can be assigned
    // to -infinity or +infinity. We have to scan the r vector and search a
    // value j such as b(0, k) <= Sum A(k, R[j]) < b(1, k).
    //
    void affect_variables(int k,
                          int selected,
                          int r_size,
                          const floatingpointT kappa,
                          const floatingpointT delta) noexcept
    {
        if (selected < 0) {
            for (int i = 0; i != r_size; ++i) {
                x(R[i].id) = 0;
                ap.add_p(k, R[i].id, -delta);
            }
        } else if (selected + 1 >= r_size) {
            for (int i = 0; i != r_size; ++i) {
                x(R[i].id) = 1;
                ap.add_p(k, R[i].id, delta);
            }
        } else {
            pi(k) += ((R[selected].value + R[selected + 1].value) / 2.0);

            floatingpoint_type d =
              delta + ((kappa / (1.0 - kappa)) *
                       (R[selected + 1].value - R[selected].value));

            int i = 0;
            for (; i <= selected; ++i) {
                x(R[i].id) = 1;
                ap.add_p(k, R[i].id, +d);
            }

            for (; i != r_size; ++i) {
                x(R[i].id) = 0;
                ap.add_p(k, R[i].id, -d);
            }
        }
    }

    void push_and_compute_update_row(int k,
                                     floatingpoint_type kappa,
                                     floatingpoint_type delta,
                                     floatingpoint_type theta,
                                     floatingpoint_type obj_amp)
    {
        if (Z[k]) {
            if (b(k).min == b(k).max)
                compute_update_row_Z_eq(
                  k, b(k).min, kappa, delta, theta, obj_amp);
            else
                compute_update_row_Z_ineq(
                  k, b(k).min, b(k).max, kappa, delta, theta, obj_amp);
        } else if (!C[k]) {
            if (b(k).min == b(k).max)
                compute_update_row_01_eq(
                  k, b(k).min, kappa, delta, theta, obj_amp);
            else
                compute_update_row_01_ineq(
                  k, b(k).min, b(k).max, kappa, delta, theta, obj_amp);
        } else {
            if (b(k).min == b(k).max)
                compute_update_row_101_eq(
                  k, b(k).min, kappa, delta, theta, obj_amp);
            else
                compute_update_row_101_ineq(
                  k, b(k).min, b(k).max, kappa, delta, theta, obj_amp);
        }
    }

    void compute_update_row(int k,
                            floatingpoint_type kappa,
                            floatingpoint_type delta,
                            floatingpoint_type theta)
    {
        if (Z[k]) {
            if (b(k).min == b(k).max)
                compute_update_row_Z_eq(k, b(k).min, kappa, delta, theta, 0);
            else
                compute_update_row_Z_ineq(
                  k, b(k).min, b(k).max, kappa, delta, theta, 0);
        } else if (!C[k]) {
            if (b(k).min == b(k).max)
                compute_update_row_01_eq(k, b(k).min, kappa, delta, theta, 0);
            else
                compute_update_row_01_ineq(
                  k, b(k).min, b(k).max, kappa, delta, theta, 0);
        } else {
            if (b(k).min == b(k).max)
                compute_update_row_101_eq(k, b(k).min, kappa, delta, theta, 0);
            else
                compute_update_row_101_ineq(
                  k, b(k).min, b(k).max, kappa, delta, theta, 0);
        }
    }
};

// #ifndef BARYONYX_FULL_OPTIMIZATION
// template <typename floatingpointT>
// static void
// print_AP(std::shared_ptr<bx::context> ctx,
//          const AP_type<floatingpointT>& ap,
//          int k,
//          int rows,
//          int cols)
// {
//     int level = ctx->get_integer_parameter("print-level", 0l);
//     if (level <= 1)
//         return;

//     ctx->debug("P after constraint %d computation:\n", k);
//     std::vector<floatingpointT> to_show(cols);

//     for (int i{ 0 }; i != rows; ++i) {
//         std::fill(std::begin(to_show),
//                   std::end(to_show),
//                   std::numeric_limits<floatingpointT>::infinity());

//         auto its = ap.row(i);

//         for (; std::get<0>(its) != std::get<1>(its); ++std::get<0>(its))
//             to_show[std::get<0>(its)->position] =
//               ap.P()[std::get<0>(its)->value];

//         for (auto elem : to_show)
//             if (elem == std::numeric_limits<floatingpointT>::infinity())
//                 ctx->debug("          ");
//             else
//                 ctx->debug("%+.6f ", (double)elem);

//         ctx->debug("\n");
//     }
// }
// #endif

template<typename floatingpointT, typename randomT>
struct compute_none
{
    using random_type = randomT;

    std::shared_ptr<bx::context> m_ctx;
    std::vector<int> R;

    template<typename solverT>
    compute_none(std::shared_ptr<bx::context> ctx, const solverT& s, randomT&)
      : m_ctx(std::move(ctx))
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

    std::shared_ptr<bx::context> m_ctx;
    std::vector<int> R;
    int nb = 0;

    template<typename solverT>
    compute_reversing(std::shared_ptr<bx::context> ctx, solverT& s, randomT&)
      : m_ctx(std::move(ctx))
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

    std::shared_ptr<bx::context> m_ctx;
    std::vector<int> R;
    random_type& rng;

    template<typename solverT>
    compute_random(std::shared_ptr<bx::context> ctx,
                   solverT& s,
                   random_type& rng_)
      : m_ctx(std::move(ctx))
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
static void
sort(iteratorT begin, iteratorT end, compute_infeasibility_incr)
{
    std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
        return lhs.second < rhs.second;
    });
}

template<typename iteratorT>
static void
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

    std::shared_ptr<bx::context> m_ctx;
    std::vector<std::pair<int, int>> R;
    random_type& rng;

    template<typename solverT>
    compute_infeasibility(std::shared_ptr<bx::context> ctx,
                          solverT& s,
                          random_type& rng_)
      : m_ctx(std::move(ctx))
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
                     solver.x(std::get<0>(ak)->position);

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
        ::sort(R.begin(), R.end(), direction_type());

        for (auto it = R.begin(), et = R.end(); it != et; ++it)
            solver.compute_update_row(it->first, kappa, delta, theta);

        return local_compute_missing_constraint(solver);
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
    using random_type = randomT;

    std::chrono::time_point<std::chrono::steady_clock> m_begin;
    std::chrono::time_point<std::chrono::steady_clock> m_end;

    std::shared_ptr<bx::context> m_ctx;
    randomT& m_rng;
    const std::vector<std::string>& m_variable_names;
    const bx::affected_variables& m_affected_vars;

    x_type m_best_x;
    bx::result m_best;

    solver_functor(std::shared_ptr<bx::context> ctx,
                   randomT& rng,
                   const std::vector<std::string>& variable_names,
                   const bx::affected_variables& affected_vars)
      : m_ctx(std::move(ctx))
      , m_rng(rng)
      , m_variable_names(variable_names)
      , m_affected_vars(affected_vars)
    {}

    bx::result operator()(
      const std::vector<bx::itm::merged_constraint>& constraints,
      int variables,
      const c_type<floatingpointT>& original_costs,
      const c_type<floatingpointT>& norm_costs,
      double cost_constant,
      const bx::itm::parameters& p)
    {
        m_begin = std::chrono::steady_clock::now();
        m_end = m_begin;

        int i = 0;
        int pushed = -1;
        int best_remaining = -1;
        int pushing_iteration = p.pushing_iteration_limit;
        floatingpoint_type kappa = p.kappa_min;

        solver<floatingpoint_type, mode_type, random_type> slv(
          m_rng, variables, norm_costs, constraints);

        constraint_order_type compute(m_ctx, slv, m_rng);

        m_ctx->info("* solver starts:\n");

        for (;;) {
            int remaining = compute.run(slv, kappa, p.delta, p.theta);

            if (best_remaining == -1 or remaining < best_remaining) {
                best_remaining = remaining;
                m_best = slv.results(original_costs, cost_constant);
                m_best.loop = i;
                m_best.remaining_constraints = remaining;
                m_best.duration =
                  std::chrono::duration_cast<std::chrono::duration<double>>(
                    m_end - m_begin)
                    .count();

                m_ctx->info(
                  "  - constraints remaining: %d/%d at %fs (loop: %d)\n",
                  remaining,
                  m_best.constraints,
                  m_best.duration,
                  i);
            }

#ifndef BARYONYX_FULL_OPTIMIZATION
            slv.print(m_ctx, m_variable_names, p.print_level);
#endif

            if (m_best.status == bx::result_status::success) {
                ++pushing_iteration;

                if (pushed == -1)
                    m_ctx->info("  - start push system:\n");

                if (pushing_iteration >= p.pushing_iteration_limit) {
                    pushed++;
                    pushing_iteration = 0;

                    m_ctx->info(
                      "    - push %d: kappa * k: %f objective amplifier: %f\n",
                      pushed,
                      static_cast<double>(p.pushing_k_factor * kappa),
                      static_cast<double>(p.pushing_objective_amplifier));

                    remaining =
                      compute.push_and_run(slv,
                                           p.pushing_k_factor * kappa,
                                           p.delta,
                                           p.theta,
                                           p.pushing_objective_amplifier);

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
                    m_ctx->info(
                      "    - Push system limit reached. Solution found: %f\n",
                      m_best.value);
                    return m_best;
                }
            }

            if (i > p.w)
                kappa += p.kappa_step *
                         std::pow(static_cast<floatingpointT>(remaining) /
                                    static_cast<floatingpointT>(slv.m),
                                  p.alpha);

            if (++i > p.limit) {
                m_ctx->info("  - Loop limit reached: %d\n", i);
                if (pushed == -1)
                    m_best.status = bx::result_status::limit_reached;

                if (m_ctx->get_integer_parameter("print-level", 0) > 0)
                    print_missing_constraint(m_ctx,
                                             slv.ap,
                                             m_best.variable_value,
                                             slv.b,
                                             m_variable_names);

                return m_best;
            }

            if (kappa > p.kappa_max) {
                m_ctx->info("  - Kappa max reached: %+.6f\n", (double)kappa);
                if (pushed == -1)
                    m_best.status = bx::result_status::kappa_max_reached;

                if (m_ctx->get_integer_parameter("print-level", 0) > 0)
                    print_missing_constraint(m_ctx,
                                             slv.ap,
                                             m_best.variable_value,
                                             slv.b,
                                             m_variable_names);

                return m_best;
            }

            m_end = std::chrono::steady_clock::now();
            if (bx::is_time_limit(p.time_limit, m_begin, m_end)) {
                m_ctx->info(
                  "  - Time limit reached: %d %+.6f\n", i, (double)kappa);
                if (pushed == -1)
                    m_best.status = bx::result_status::time_limit_reached;

                if (m_ctx->get_integer_parameter("print-level", 0) > 0)
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

            m_ctx->info("  - Solution found: %f (i=%d t=%fs)\n",
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

template<typename floatingpointT,
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

    std::shared_ptr<bx::context> m_ctx;
    randomT m_rng;
    int m_thread_id;
    const std::vector<std::string>& m_variable_names;
    const bx::affected_variables& m_affected_vars;
    x_type m_best_x;
    bx::result m_best;

    optimize_functor(std::shared_ptr<bx::context> ctx,
                     int thread_id,
                     typename random_type::result_type seed,
                     const std::vector<std::string>& variable_names,
                     const bx::affected_variables& affected_vars)
      : m_ctx(std::move(ctx))
      , m_rng(seed)
      , m_thread_id(thread_id)
      , m_variable_names(variable_names)
      , m_affected_vars(affected_vars)
    {}

    bx::result operator()(
      const std::vector<bx::itm::merged_constraint>& constraints,
      int variables,
      const c_type<floatingpointT>& original_costs,
      const c_type<floatingpointT>& norm_costs,
      double cost_constant,
      const bx::itm::parameters& p)
    {
        m_begin = std::chrono::steady_clock::now();
        m_end = m_begin;

        int i = 0;
        int pushed = -1;
        int pushing_iteration = 0;
        floatingpoint_type kappa = p.kappa_min;

        solver<floatingpoint_type, mode_type, random_type> slv(
          m_rng, variables, norm_costs, constraints);

        constraint_order_type compute(m_ctx, slv, m_rng);

        for (; not bx::is_time_limit(p.time_limit, m_begin, m_end);
             m_end = std::chrono::steady_clock::now(), ++i) {

            int remaining = compute.run(slv, kappa, p.delta, p.theta);

            if (remaining == 0) {
                auto current = slv.results(original_costs, cost_constant);
                current.loop = i;
                current.remaining_constraints = remaining;
                if (store_if_better(current)) {
                    m_best_x = slv.x;
                    pushed = 0;
                }
            }

            if (i > p.w)
                kappa +=
                  p.kappa_step * std::pow(static_cast<double>(remaining) /
                                            static_cast<double>(slv.m),
                                          p.alpha);

            if (i >= p.limit or kappa > p.kappa_max or
                pushed > p.pushes_limit) {
                slv.reinit(m_best_x, p.reverse_solution);

                i = 0;
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

                    remaining =
                      compute.push_and_run(slv,
                                           p.pushing_k_factor * kappa,
                                           p.delta,
                                           p.theta,
                                           p.pushing_objective_amplifier);

                    if (remaining == 0) {
                        auto current =
                          slv.results(original_costs, cost_constant);
                        current.loop = i;
                        current.remaining_constraints = 0;

                        if (store_if_better(current))
                            m_best_x = slv.x;
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
 * input vector is too small or with infinity value, the c is unchanged.
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

template<typename floatingpointT>
static inline c_type<floatingpointT>
make_objective_function(const bx::objective_function& obj, int n)
{
    c_type<floatingpointT> ret(n, 0);

    for (const auto& elem : obj.elements)
        ret(elem.variable_index) += elem.factor;

    return ret;
}

template<typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT>
static bx::result
solve(std::shared_ptr<bx::context> ctx,
      bx::problem& pb,
      const bx::itm::parameters& p)
{
    ctx->info("Solver initializing\n");

    bx::result ret;
    auto affected_vars = std::move(pb.affected_vars);

    auto constraints{ bx::itm::make_merged_constraints(ctx, pb, p) };
    if (not constraints.empty() and not pb.vars.values.empty()) {
        randomT rng(ctx->get_integer_parameter(
          "seed",
          std::chrono::system_clock::now().time_since_epoch().count()));

        auto variables = bx::numeric_cast<int>(pb.vars.values.size());
        auto cost =
          make_objective_function<floatingpointT>(pb.objective, variables);
        auto norm_costs = normalize_costs(ctx, p.norm, cost, rng);
        auto cost_constant = pb.objective.value;
        auto names = std::move(pb.vars.names);

        bx::clear(pb);

        solver_functor<floatingpointT, modeT, constraintOrderT, randomT> slv(
          ctx, rng, names, affected_vars);

        ret = slv(constraints, variables, cost, norm_costs, cost_constant, p);

        ret.method = "inequalities_Zcoeff solver";
        ret.variable_name = std::move(names);
    } else {
        ret.status = bx::result_status::success;
    }
    ret.affected_vars = std::move(affected_vars);

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
         int thread)
{
    bx::Expects(thread >= 1, "optimize: bad thread number");

    ctx->info("Optimizer initializing\n");

    bx::result ret;
    auto affected_vars = std::move(pb.affected_vars);

    auto constraints{ bx::itm::make_merged_constraints(ctx, pb, p) };
    if (not constraints.empty() and not pb.vars.values.empty()) {

        randomT rng(ctx->get_integer_parameter(
          "seed",
          std::chrono::system_clock::now().time_since_epoch().count()));

        auto variables = bx::numeric_cast<int>(pb.vars.values.size());
        auto cost =
          make_objective_function<floatingpointT>(pb.objective, variables);
        auto norm_costs = normalize_costs(ctx, p.norm, cost, rng);
        auto cost_constant = pb.objective.value;
        auto names = std::move(pb.vars.names);

        bx::clear(pb);

        std::vector<std::thread> pool(thread);
        pool.clear();
        std::vector<std::future<bx::result>> results(thread);
        results.clear();

        if (thread == 1)
            ctx->info("optimizer starts with one thread\n");
        else
            ctx->info("Optimizer starts with %d threads\n", thread);

        std::uniform_int_distribution<typename randomT::result_type> dst(
          std::numeric_limits<typename randomT::result_type>::min(),
          std::numeric_limits<typename randomT::result_type>::max());

        for (int i{ 0 }; i != thread; ++i) {
            auto seed = bx::numeric_cast<typename randomT::result_type>(dst(rng));

            std::packaged_task<bx::result()> task(std::bind(
              optimize_functor<floatingpointT,
                               modeT,
                               constraintOrderT,
                               randomT>(ctx, i, seed, names, affected_vars),
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
            if (current.status == bx::result_status::success) {
                if (ret.status != bx::result_status::success or
                    is_better_solution(current.value, ret.value, modeT()))
                    ret = current;
            }
        }

        ret.method = "inequalities_Zcoeff optimizer";
        ret.variable_name = std::move(names);
    } else {
        ret.status = bx::result_status::success;
    }
    ret.affected_vars = std::move(affected_vars);

    return ret;
}

template<typename realT, typename modeT, typename randomT>
static bx::result
dispatch_solve(std::shared_ptr<bx::context> ctx,
               bx::problem& pb,
               const bx::itm::parameters& p)
{
    switch (p.order) {
    case bx::itm::constraint_order::none:
        return ::solve<realT, modeT, ::compute_none<realT, randomT>, randomT>(
          ctx, pb, p);
    case bx::itm::constraint_order::reversing:
        return ::
          solve<realT, modeT, ::compute_reversing<realT, randomT>, randomT>(
            ctx, pb, p);
    case bx::itm::constraint_order::random_sorting:
        return ::solve<realT, modeT, ::compute_random<realT, randomT>, randomT>(
          ctx, pb, p);
    case bx::itm::constraint_order::infeasibility_decr:
        return ::solve<
          realT,
          modeT,
          ::compute_infeasibility<realT, randomT, ::compute_infeasibility_decr>, randomT>(
          ctx, pb, p);
    case bx::itm::constraint_order::infeasibility_incr:
        return ::solve<
          realT,
          modeT,
          ::compute_infeasibility<realT, randomT, ::compute_infeasibility_incr>, randomT>(
          ctx, pb, p);
    }

    return {};
}

template<typename realT, typename modeT, typename randomT>
static bx::result
dispatch_optimize(std::shared_ptr<bx::context> ctx,
                  bx::problem& pb,
                  const bx::itm::parameters& p,
                  int thread)
{
    switch (p.order) {
    case bx::itm::constraint_order::none:
        return ::
          optimize<realT, modeT, ::compute_none<realT, randomT>, randomT>(
            ctx, pb, p, thread);
    case bx::itm::constraint_order::reversing:
        return ::
          optimize<realT, modeT, ::compute_reversing<realT, randomT>, randomT>(
            ctx, pb, p, thread);
    case bx::itm::constraint_order::random_sorting:
        return ::optimize<realT, modeT, ::compute_random<realT, randomT>, randomT>(
          ctx, pb, p, thread);
    case bx::itm::constraint_order::infeasibility_decr:
        return ::optimize<
          realT,
          modeT,
          ::compute_infeasibility<realT, randomT, ::compute_infeasibility_decr>, randomT>(
          ctx, pb, p, thread);
    case bx::itm::constraint_order::infeasibility_incr:
        return ::optimize<
          realT,
          modeT,
          ::compute_infeasibility<realT, randomT, ::compute_infeasibility_incr>, randomT>(
          ctx, pb, p, thread);
    }

    return {};
}

} // anonymous namespace

namespace baryonyx {
namespace itm {

result
inequalities_Zcoeff_wedelin_solve(
  const std::shared_ptr<baryonyx::context>& ctx,
  problem& pb)
{
    ctx->info("inequalities_Zcoeff_wedelin_solve\n");
    parameters p(ctx);

    using random_type = std::default_random_engine;
    // using random_type = std::mt19937;

    if (pb.type == baryonyx::objective_function_type::maximize) {
        switch (p.float_type) {
        case floating_point_type::float_type:
            return dispatch_solve<float, maximize_tag, random_type>(
              ctx, pb, p);
        case floating_point_type::double_type:
            return dispatch_solve<double, maximize_tag, random_type>(
              ctx, pb, p);
        case floating_point_type::longdouble_type:
            return dispatch_solve<long double, maximize_tag, random_type>(
              ctx, pb, p);
        }
    } else {
        switch (p.float_type) {
        case floating_point_type::float_type:
            return dispatch_solve<float, minimize_tag, random_type>(
              ctx, pb, p);
        case floating_point_type::double_type:
            return dispatch_solve<double, minimize_tag, random_type>(
              ctx, pb, p);
        case floating_point_type::longdouble_type:
            return dispatch_solve<long double, minimize_tag, random_type>(
              ctx, pb, p);
        }
    }

    return {};
}

result
inequalities_Zcoeff_wedelin_optimize(
  const std::shared_ptr<baryonyx::context>& ctx,
  problem& pb,
  int thread)
{
    ctx->info("inequalities_Zcoeff_wedelin_optimize\n");
    parameters p(ctx);

    using random_type = std::default_random_engine;
    // using random_type = std::mt19937;

    if (pb.type == baryonyx::objective_function_type::maximize) {
        switch (p.float_type) {
        case floating_point_type::float_type:
            return dispatch_optimize<float, maximize_tag, random_type>(
              ctx, pb, p, thread);
        case floating_point_type::double_type:
            return dispatch_optimize<double, maximize_tag, random_type>(
              ctx, pb, p, thread);
        case floating_point_type::longdouble_type:
            return dispatch_optimize<long double, maximize_tag, random_type>(
              ctx, pb, p, thread);
        }
    } else {
        switch (p.float_type) {
        case floating_point_type::float_type:
            return dispatch_optimize<float, minimize_tag, random_type>(
              ctx, pb, p, thread);
        case floating_point_type::double_type:
            return dispatch_optimize<double, minimize_tag, random_type>(
              ctx, pb, p, thread);
        case floating_point_type::longdouble_type:
            return dispatch_optimize<long double, minimize_tag, random_type>(
              ctx, pb, p, thread);
        }
    }

    return {};
}
}
}
