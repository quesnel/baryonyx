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

template<typename floatingpointT>
using AP_type = bx::SparseArray<int, floatingpointT>;

using b_type = baryonyx::fixed_array<bound>;

template<typename floatingpointT>
using c_type = baryonyx::fixed_array<floatingpointT>;
using x_type = baryonyx::fixed_array<std::int8_t>;

template<typename floatingpointT>
using pi_type = baryonyx::fixed_array<floatingpointT>;

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
           const std::vector<bx::itm::merged_constraint>& csts,
           bx::itm::init_policy_type init_type,
           double init_random)
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
                    ap.set(i,
                           cst.variable_index,
                           cst.factor,
                           static_cast<floatingpoint_type>(0.0));

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

        x_type empty;
        reinit(empty, init_type, init_random);
    }

    void reinit(const x_type& best_previous,
                bx::itm::init_policy_type type,
                double init_random)
    {
        std::fill(ap.P().begin(), ap.P().end(), 0);
        std::fill(pi.begin(), pi.end(), 0);

        //
        // Default, we randomly change the init policy using using the
        // bernoulli distribution with p = 0.1.
        //

        {
            std::bernoulli_distribution d(0.1);

            if (d(rng)) {
                std::uniform_int_distribution<int> di(0, 2);
                auto ret = di(rng);
                type = (ret == 0)
                         ? bx::itm::init_policy_type::best
                         : (ret == 1) ? bx::itm::init_policy_type::random
                                      : bx::itm::init_policy_type::best;
            }
        }

        //
        // If no solution was found previously and the policy type is best, we
        // randomly replace the best policy init type with random or bastert
        // policy solution using the berouilli distribution with p = 0.5.
        //

        if (best_previous.empty() and
            type == bx::itm::init_policy_type::best) {
            std::bernoulli_distribution d(0.5);
            if (type == bx::itm::init_policy_type::best) {
                type = d(rng) ? bx::itm::init_policy_type::random
                              : bx::itm::init_policy_type::bastert;
            }
        }

        init_random = bx::clamp(init_random, 0.0, 1.0);

        std::bernoulli_distribution d(init_random);

        switch (type) {
        case bx::itm::init_policy_type::bastert:
            if (init_random == 0.0 or init_random == 1.0) {
                bool value_if_cost_0 = init_random == 1.0;

                for (int i = 0, e = n; i != e; ++i)
                    x(i) = init_x(c(i), value_if_cost_0, mode_type());
            } else {
                for (int i = 0, e = n; i != e; ++i)
                    x(i) = init_x(c(i), d(rng), mode_type());
            }
            break;
        case bx::itm::init_policy_type::random:
            for (int i = 0; i != n; ++i)
                x(i) = d(rng);
            break;
        case bx::itm::init_policy_type::best:
            for (int i = 0; i != n; ++i)
                x(i) = (d(rng)) ? best_previous(i) : d(rng);
            break;
        }
    }

    void print(const std::shared_ptr<bx::context>& ctx,
               const std::vector<std::string>& names,
               int print_level) const
    {
        if (print_level <= 0)
            return;

        debug(ctx, "  - X: {} to {}\n", 0, length(x));
        for (int i = 0, e = length(x); i != e; ++i)
            debug(ctx,
                  "    - {} {}={}/c_i:{}\n",
                  i,
                  names[i],
                  static_cast<int>(x[i]),
                  c[i]);
        debug(ctx, "\n");

        for (int k = 0, ek = m; k != ek; ++k) {
            auto ak = ap.row(k);
            int v = 0;

            for (; std::get<0>(ak) != std::get<1>(ak); ++std::get<0>(ak))
                v += ap.A()[std::get<0>(ak)->value] *
                     x(std::get<0>(ak)->position);

            bool valid = b(k).min <= v and v <= b(k).max;
            debug(ctx,
                  "C {}:{} (Lmult: {})\n",
                  k,
                  (valid ? "   valid" : "violated"),
                  pi[k]);
        }
    }

    bx::result results(const c_type<floatingpoint_type>& original_costs,
                       const double cost_constant) const
    {
        bx::result ret;

        if (is_valid_solution(ap, x, b)) {
            ret.status = bx::result_status::success;
            double value = static_cast<double>(cost_constant);

            for (int i{ 0 }, ei{ n }; i != ei; ++i)
                value += static_cast<double>(original_costs[i] * x[i]);

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
            x[ck[i].id_A] = !x[ck[i].id_A]; // 1 - x[ck[i].id_A];
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
            x[ck[i].id_A] = !x[ck[i].id_A]; // 1 - x[ck[i].id_A];
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
            x[ck[i].id_A] = !x[ck[i].id_A]; // 1 - x[ck[i].id_A];
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
            x[ck[i].id_A] = !x[ck[i].id_A]; // 1 - x[ck[i].id_A];
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
                auto a = static_cast<floatingpoint_type>(ap.A()[ht->value]);
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
            pi(k) += ((R[selected].value + R[selected + 1].value) /
                      static_cast<floatingpoint_type>(2.0));

            floatingpoint_type d =
              delta +
              ((kappa / (static_cast<floatingpoint_type>(1.0) - kappa)) *
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
                compute_update_row_Z_eq(k,
                                        b(k).min,
                                        kappa,
                                        delta,
                                        theta,
                                        static_cast<floatingpoint_type>(0));
            else
                compute_update_row_Z_ineq(k,
                                          b(k).min,
                                          b(k).max,
                                          kappa,
                                          delta,
                                          theta,
                                          static_cast<floatingpoint_type>(0));
        } else if (!C[k]) {
            if (b(k).min == b(k).max)
                compute_update_row_01_eq(k,
                                         b(k).min,
                                         kappa,
                                         delta,
                                         theta,
                                         static_cast<floatingpoint_type>(0));
            else
                compute_update_row_01_ineq(k,
                                           b(k).min,
                                           b(k).max,
                                           kappa,
                                           delta,
                                           theta,
                                           static_cast<floatingpoint_type>(0));
        } else {
            if (b(k).min == b(k).max)
                compute_update_row_101_eq(k,
                                          b(k).min,
                                          kappa,
                                          delta,
                                          theta,
                                          static_cast<floatingpoint_type>(0));
            else
                compute_update_row_101_ineq(
                  k,
                  b(k).min,
                  b(k).max,
                  kappa,
                  delta,
                  theta,
                  static_cast<floatingpoint_type>(0));
        }
    }
};

}
}
