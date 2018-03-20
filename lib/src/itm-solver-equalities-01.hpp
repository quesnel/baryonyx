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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_EQUALITIES_01COEFF_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_EQUALITIES_01COEFF_HPP

#include "itm-solver-common.hpp"
#include "sparse-matrix.hpp"

namespace baryonyx {
namespace itm {

template<typename floatingpointT, typename modeT, typename randomT>
struct solver_equalities_01coeff
{
    using floatingpoint_type = floatingpointT;
    using mode_type = modeT;
    using random_type = randomT;

    using AP_type = sparse_matrix<int>;
    using b_type = baryonyx::fixed_array<int>;
    using c_type = baryonyx::fixed_array<floatingpointT>;
    using pi_type = baryonyx::fixed_array<floatingpointT>;
    using P_type = fixed_array<floatingpointT>;

    random_type& rng;

    // Sparse matrix to store A and P values.
    AP_type ap;
    fixed_array<floatingpointT> P;

    // Vector shared between all constraints to store the reduced cost.
    fixed_array<r_data<floatingpoint_type>> R;

    // Bound vector.
    b_type b;
    const c_type& c;
    x_type x;
    pi_type pi;
    int m;
    int n;

    solver_equalities_01coeff(random_type& rng_,
                              int n_,
                              const c_type& c_,
                              const std::vector<itm::merged_constraint>& csts,
                              itm::init_policy_type init_type,
                              double init_random)
      : rng(rng_)
      , ap(csts, length(csts), n_)
      , P(element_number(csts), 0)
      , b(length(csts))
      , c(c_)
      , x(n_)
      , pi(length(csts))
      , m(length(csts))
      , n(n_)
    {
        {
            // Compute the minimal bounds for each constraints, default
            // constraints are -oo <= ... <= bkmax, bkmin <= ... <= +oo and
            // bkmin <= ... <= bkmax. This code remove infinity and replace
            // with minimal or maximal value of the constraint.

            for (int i = 0; i != m; ++i) {
                int lower = 0, upper = 0;

                for (const auto& cst : csts[i].elements) {
                    if (cst.factor > 0)
                        upper += cst.factor;

                    if (cst.factor < 0)
                        lower += cst.factor;
                }

                assert(lower == 0);
                assert(upper == length(csts[i].elements));
                assert(csts[i].min == csts[i].max);
                assert(lower == 0 and upper == length(csts[i].elements) and
                       csts[i].min == csts[i].max and
                       "Preprocessor select error");

                b(i) = csts[i].min;
            }
        }

        {
            //
            // Compute the R vector size and the C vectors for each constraints
            // with negative coefficient.
            //

            int rsizemax = length(csts[0].elements);
            for (int i = 1; i != m; ++i)
                if (rsizemax < length(csts[i].elements))
                    rsizemax = length(csts[i].elements);

            R = fixed_array<r_data<floatingpoint_type>>(rsizemax);
        }

        x_type empty;
        reinit(empty, init_type, init_random);
    }

    int factor(int /*value*/) const noexcept
    {
        return 1;
    }

    int bound_min(int constraint) const noexcept
    {
        return b[constraint];
    }

    int bound_max(int constraint) const noexcept
    {
        return b[constraint];
    }

    int bound_init(int constraint) const
    {
        return b[constraint];
    }

    floatingpointT compute_sum_A_pi(int variable) const
    {
        floatingpointT ret{ 0 };

        AP_type::const_col_iterator ht, hend;
        std::tie(ht, hend) = ap.column(variable);

        for (; ht != hend; ++ht)
            ret += pi[ht->row];

        return ret;
    }

    void reinit(const x_type& best_previous,
                itm::init_policy_type type,
                double init_random)
    {
        std::fill(
          P.begin(), P.end(), static_cast<typename P_type::value_type>(0));
        std::fill(
          pi.begin(), pi.end(), static_cast<typename pi_type::value_type>(0));

        //
        // Default, we randomly change the init policy using using the
        // bernoulli distribution with p = 0.1.
        //

        {
            std::bernoulli_distribution d(0.1);

            if (d(rng)) {
                std::uniform_int_distribution<int> di(0, 2);
                auto ret = di(rng);
                type = (ret == 0) ? itm::init_policy_type::best
                                  : (ret == 1) ? itm::init_policy_type::random
                                               : itm::init_policy_type::best;
            }
        }

        //
        // If no solution was found previously and the policy type is best, we
        // randomly replace the best policy init type with random or bastert
        // policy solution using the berouilli distribution with p = 0.5.
        //

        if (best_previous.empty() and type == itm::init_policy_type::best) {
            std::bernoulli_distribution d(0.5);
            if (type == itm::init_policy_type::best) {
                type = d(rng) ? itm::init_policy_type::random
                              : itm::init_policy_type::bastert;
            }
        }

        init_random = clamp(init_random, 0.0, 1.0);

        std::bernoulli_distribution d(init_random);

        switch (type) {
        case itm::init_policy_type::bastert:
            if (init_random == 0.0 or init_random == 1.0) {
                bool value_if_cost_0 = init_random == 1.0;

                for (int i = 0, e = n; i != e; ++i)
                    x[i] = init_x(c(i), value_if_cost_0, mode_type());
            } else {
                for (int i = 0, e = n; i != e; ++i)
                    x[i] = init_x(c(i), d(rng), mode_type());
            }
            break;
        case itm::init_policy_type::random:
            for (int i = 0; i != n; ++i)
                x[i] = d(rng);
            break;
        case itm::init_policy_type::best:
            for (int i = 0; i != n; ++i) {
                x[i] = (d(rng)) ? (best_previous[i]) : d(rng);
            }
            break;
        }
    }

    void print(const context_ptr& ctx,
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
            typename AP_type::const_row_iterator it, et;

            std::tie(it, et) = ap.row(k);
            int v = 0;

            for (; it != et; ++it)
                v += x[it->column];

            bool valid = b(k) == v;
            debug(ctx,
                  "C {}:{} (Lmult: {})\n",
                  k,
                  (valid ? "   valid" : "violated"),
                  pi[k]);
        }
    }

    bool is_valid_solution() const
    {
        typename AP_type::const_row_iterator it, et;

        for (int k = 0; k != m; ++k) {
            std::tie(it, et) = ap.row(k);
            int v = 0;

            for (; it != et; ++it)
                v += x[it->column];

            if (b[k] != v)
                return false;
        }

        return true;
    }

    template<typename Container>
    int compute_violated_constraints(Container& c) const
    {
        typename AP_type::const_row_iterator it, et;

        c.clear();

        for (int k = 0; k != m; ++k) {
            std::tie(it, et) = ap.row(k);
            int v = 0;

            for (; it != et; ++it)
                v += x[it->column];

            if (b[k] != v)
                c.emplace_back(k);
        }

        return length(c);
    }

    double results(const c_type& original_costs,
                   const double cost_constant) const
    {
        assert(is_valid_solution());

        double value = static_cast<double>(cost_constant);

        for (int i{ 0 }, ei{ n }; i != ei; ++i)
            value += static_cast<double>(original_costs[i] * x[i]);

        return value;
    }

    typename AP_type::row_iterator ap_value(typename AP_type::row_iterator it,
                                            int id_in_r)
    {
        return it + id_in_r;
    }

    //
    // Decrease influence of local preferences. 0 will completely reset the
    // preference values for the current row. > 0 will keep former decision
    // in mind.
    //
    template<typename iteratorT>
    void decrease_preference(iteratorT begin,
                             iteratorT end,
                             floatingpoint_type theta) noexcept
    {
        for (; begin != end; ++begin)
            P[begin->value] *= theta;
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

            typename AP_type::const_col_iterator ht, hend;
            std::tie(ht, hend) = ap.column(begin->column);

            for (; ht != hend; ++ht) {
                sum_a_pi += pi[ht->row];
                sum_a_p += P[ht->value];
            }

            R[r_size].id = r_size;
            R[r_size].value = c[begin->column] - sum_a_pi - sum_a_p;
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

    //
    // The bkmin and bkmax constraint bounds are not equal and can be
    // assigned to -infinity or +infinity. We have to scan the r vector and
    // search a value j such as b(0, k) <= Sum A(k, R[j]) < b(1, k).
    //
    template<typename Iterator>
    void affect_variables(Iterator it,
                          int k,
                          int selected,
                          int r_size,
                          const floatingpointT kappa,
                          const floatingpointT delta) noexcept
    {
        if (selected < 0) {
            for (int i = 0; i != r_size; ++i) {
                auto var = ap_value(it, R[i].id);

                x[var->column] = 0;
                P[var->value] -= delta;
            }
        } else if (selected + 1 >= r_size) {
            for (int i = 0; i != r_size; ++i) {
                auto var = ap_value(it, R[i].id);

                x[var->column] = 1;
                P[var->value] += delta;
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
                auto var = ap_value(it, R[i].id);

                x[var->column] = 1;
                P[var->value] += d;
            }

            for (; i != r_size; ++i) {
                auto var = ap_value(it, R[i].id);

                x[var->column] = 0;
                P[var->value] -= d;
            }
        }
    }

    void push_and_compute_update_row(int k,
                                     floatingpoint_type kappa,
                                     floatingpoint_type delta,
                                     floatingpoint_type theta,
                                     floatingpoint_type objective_amplifier)
    {
        typename AP_type::row_iterator it, et;
        std::tie(it, et) = ap.row(k);

        decrease_preference(it, et, theta);
        const int r_size = compute_reduced_costs(it, et);

        if (objective_amplifier)
            for (int i = 0; i != r_size; ++i)
                R[i].value += objective_amplifier * c[R[i].id];

        calculator_sort(R.begin(), R.begin() + r_size, rng, mode_type());
        int selected = select_variables_equality(r_size, b[k]);
        affect_variables(it, k, selected, r_size, kappa, delta);
    }

    void compute_update_row(int k,
                            floatingpoint_type kappa,
                            floatingpoint_type delta,
                            floatingpoint_type theta)
    {
        typename AP_type::row_iterator it, et;
        std::tie(it, et) = ap.row(k);

        decrease_preference(it, et, theta);
        const int r_size = compute_reduced_costs(it, et);
        calculator_sort(R.begin(), R.begin() + r_size, rng, mode_type());
        int selected = select_variables_equality(r_size, b[k]);
        affect_variables(it, k, selected, r_size, kappa, delta);
    }
};

} // namespace itm
} // namespace baryonyx

#endif
