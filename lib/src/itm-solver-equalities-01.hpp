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

#include <memory>

namespace baryonyx {
namespace itm {

template<typename floatingpointT, typename modeT, typename randomT>
struct solver_equalities_01coeff
{
    using floatingpoint_type = floatingpointT;
    using mode_type = modeT;
    using random_type = randomT;

    random_type& rng;

    sparse_matrix<int> ap;
    std::vector<bool> x;
    std::unique_ptr<floatingpointT[]> P;
    std::unique_ptr<r_data<floatingpoint_type>[]> R;
    std::unique_ptr<int[]> b;
    std::unique_ptr<floatingpointT[]> pi;

    const std::unique_ptr<floatingpointT[]>& c;

    int m;
    int n;

    solver_equalities_01coeff(random_type& rng_,
                              int m_,
                              int n_,
                              const std::unique_ptr<floatingpointT[]>& c_,
                              const std::vector<itm::merged_constraint>& csts,
                              itm::init_policy_type init_type,
                              double init_random)
      : rng(rng_)
      , ap(csts, m_, n_)
      , x(n_)
      , P(std::make_unique<floatingpointT[]>(ap.size()))
      , R(std::make_unique<r_data<floatingpoint_type>[]>(
          compute_reduced_costs_vector_size(csts)))
      , b(std::make_unique<int[]>(m_))
      , pi(std::make_unique<floatingpointT[]>(m_))
      , c(c_)
      , m(m_)
      , n(n_)
    {
        //
        // Assigns bound for each constraints.
        //

        for (int i = 0; i != m; ++i) {
            for (const auto& cst : csts[i].elements)
                Ensures(cst.factor == 1);

            Ensures(csts[i].min == csts[i].max);

            b[i] = csts[i].min;
        }

        x_type empty;
        init_solver(*this, empty, init_type, init_random);
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

        sparse_matrix<int>::const_col_iterator ht, hend;
        std::tie(ht, hend) = ap.column(variable);

        for (; ht != hend; ++ht)
            ret += pi[ht->row];

        return ret;
    }

    bool is_valid_solution() const
    {
        for (int k = 0; k != m; ++k) {
            sparse_matrix<int>::const_row_iterator it, et;

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
        typename sparse_matrix<int>::const_row_iterator it, et;

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

    double results(const std::unique_ptr<floatingpointT[]>& original_costs,
                   const double cost_constant) const
    {
        assert(is_valid_solution());

        double value = static_cast<double>(cost_constant);

        for (int i{ 0 }, ei{ n }; i != ei; ++i)
            value += static_cast<double>(original_costs[i] * x[i]);

        return value;
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

            sparse_matrix<int>::const_col_iterator ht, hend;
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

        assert(bk <= r_size && "b[k] can not be reached, this is an "
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
                auto var = it + R[i].id;

                x[var->column] = false;
                P[var->value] -= delta;
            }
        } else if (selected + 1 >= r_size) {
            pi[k] += R[selected].value;

            for (int i = 0; i != r_size; ++i) {
                auto var = it + R[i].id;

                x[var->column] = true;
                P[var->value] += delta;
            }
        } else {
            pi[k] += ((R[selected].value + R[selected + 1].value) /
                      static_cast<floatingpoint_type>(2.0));

            floatingpoint_type d =
              delta +
              ((kappa / (static_cast<floatingpoint_type>(1.0) - kappa)) *
               (R[selected + 1].value - R[selected].value));

            int i = 0;
            for (; i <= selected; ++i) {
                auto var = it + R[i].id;

                x[var->column] = true;
                P[var->value] += d;
            }

            for (; i != r_size; ++i) {
                auto var = it + R[i].id;

                x[var->column] = false;
                P[var->value] -= d;
            }
        }
    }

    template<typename Iterator>
    void push_and_compute_update_row(Iterator first,
                                     Iterator last,
                                     floatingpoint_type kappa,
                                     floatingpoint_type delta,
                                     floatingpoint_type theta,
                                     floatingpoint_type objective_amplifier)
    {
        for (; first != last; ++first) {
            auto k = constraint(first);

            sparse_matrix<int>::row_iterator it, et;
            std::tie(it, et) = ap.row(k);

            decrease_preference(it, et, theta);
            const int r_size = compute_reduced_costs(it, et);

            if (objective_amplifier)
                for (int i = 0; i != r_size; ++i)
                    R[i].value += objective_amplifier * c[R[i].id];

            calculator_sort(R.get(), R.get() + r_size, rng, mode_type());
            int selected = select_variables_equality(r_size, b[k]);
            affect_variables(it, k, selected, r_size, kappa, delta);
        }
    }

    template<typename Iterator>
    void compute_update_row(Iterator first,
                            Iterator last,
                            floatingpoint_type kappa,
                            floatingpoint_type delta,
                            floatingpoint_type theta)
    {
        for (; first != last; ++first) {
            auto k = constraint(first);

            sparse_matrix<int>::row_iterator it, et;
            std::tie(it, et) = ap.row(k);

            decrease_preference(it, et, theta);
            const int r_size = compute_reduced_costs(it, et);
            calculator_sort(R.get(), R.get() + r_size, rng, mode_type());
            int selected = select_variables_equality(r_size, b[k]);
            affect_variables(it, k, selected, r_size, kappa, delta);
        }
    }
};

} // namespace itm
} // namespace baryonyx

#endif
