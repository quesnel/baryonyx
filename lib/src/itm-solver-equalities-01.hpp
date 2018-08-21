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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_SOLVER_EQUALITIES_01_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_SOLVER_EQUALITIES_01_HPP

#include "itm-common.hpp"

namespace baryonyx {
namespace itm {

template<typename Float, typename Mode, typename Random>
struct solver_equalities_01coeff
{
    using mode_type = Mode;
    using float_type = Float;

    Random& rng;

    sparse_matrix<int> ap;
    std::vector<bool> x;
    std::unique_ptr<Float[]> P;
    std::unique_ptr<r_data<Float>[]> R;
    std::unique_ptr<int[]> b;
    std::unique_ptr<Float[]> pi;

    const std::unique_ptr<Float[]>& c;

    int m;
    int n;

    solver_equalities_01coeff(Random& rng_,
                              int m_,
                              int n_,
                              const std::unique_ptr<Float[]>& c_,
                              const std::vector<merged_constraint>& csts,
                              solver_parameters::init_policy_type init_type,
                              double init_random);

    int factor(int /*value*/) const noexcept;

    int bound_min(int constraint) const noexcept;

    int bound_max(int constraint) const noexcept;

    int bound_init(int constraint) const;

    Float compute_sum_A_pi(int variable) const;

    bool is_valid_solution() const;

    int compute_violated_constraints(std::vector<int>& container) const;

    double results(const std::unique_ptr<Float[]>& original_costs,
                   const double cost_constant) const;

    //
    // Decrease influence of local preferences. 0 will completely reset the
    // preference values for the current row. > 0 will keep former decision
    // in mind.
    //
    void decrease_preference(sparse_matrix<int>::row_iterator begin,
                             sparse_matrix<int>::row_iterator end,
                             Float theta) noexcept;

    //
    // Compute the reduced costs and return the size of the newly R vector.
    //
    int compute_reduced_costs(sparse_matrix<int>::row_iterator begin,
                              sparse_matrix<int>::row_iterator end) noexcept;

    int select_variables_equality(const int r_size, int bk);

    //
    // The bkmin and bkmax constraint bounds are not equal and can be
    // assigned to -infinity or +infinity. We have to scan the r vector and
    // search a value j such as b(0, k) <= Sum A(k, R[j]) < b(1, k).
    //
    void affect_variables(sparse_matrix<int>::row_iterator it,
                          int k,
                          int selected,
                          int r_size,
                          const Float kappa,
                          const Float delta) noexcept;

    void push_and_compute_update_row(std::vector<int>::iterator first,
                                     std::vector<int>::iterator last,
                                     Float kappa,
                                     Float delta,
                                     Float theta,
                                     Float objective_amplifier);

    void compute_update_row(std::vector<int>::iterator first,
                            std::vector<int>::iterator last,
                            Float kappa,
                            Float delta,
                            Float theta);
};

template<typename Float, typename Mode, typename Random>
solver_equalities_01coeff<Float, Mode, Random>::solver_equalities_01coeff(
  Random& rng_,
  int m_,
  int n_,
  const std::unique_ptr<Float[]>& c_,
  const std::vector<merged_constraint>& csts,
  solver_parameters::init_policy_type init_type,
  double init_random)
  : rng(rng_)
  , ap(csts, m_, n_)
  , x(n_)
  , P(std::make_unique<Float[]>(ap.size()))
  , R(std::make_unique<r_data<Float>[]>(
      compute_reduced_costs_vector_size(csts)))
  , b(std::make_unique<int[]>(m_))
  , pi(std::make_unique<Float[]>(m_))
  , c(c_)
  , m(m_)
  , n(n_)
{
    //
    // Assigns bound for each constraints.
    //

    for (int i = 0; i != m; ++i) {
        for (const auto& cst : csts[i].elements)
            bx_ensures(cst.factor == 1);

        bx_ensures(csts[i].min == csts[i].max);

        b[i] = csts[i].min;
    }

    x_type empty;
    init_solver(*this, empty, init_type, init_random);
}

template<typename Float, typename Mode, typename Random>
int
solver_equalities_01coeff<Float, Mode, Random>::factor(int /*value*/) const
  noexcept
{
    return 1;
}

template<typename Float, typename Mode, typename Random>
int
solver_equalities_01coeff<Float, Mode, Random>::bound_min(int constraint) const
  noexcept
{
    return b[constraint];
}

template<typename Float, typename Mode, typename Random>
int
solver_equalities_01coeff<Float, Mode, Random>::bound_max(int constraint) const
  noexcept
{
    return b[constraint];
}

template<typename Float, typename Mode, typename Random>
int
solver_equalities_01coeff<Float, Mode, Random>::bound_init(
  int constraint) const
{
    return b[constraint];
}

template<typename Float, typename Mode, typename Random>
Float
solver_equalities_01coeff<Float, Mode, Random>::compute_sum_A_pi(
  int variable) const
{
    Float ret{ 0 };

    sparse_matrix<int>::const_col_iterator ht, hend;
    std::tie(ht, hend) = ap.column(variable);

    for (; ht != hend; ++ht)
        ret += pi[ht->row];

    return ret;
}

template<typename Float, typename Mode, typename Random>
bool
solver_equalities_01coeff<Float, Mode, Random>::is_valid_solution() const
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

template<typename Float, typename Mode, typename Random>
int
solver_equalities_01coeff<Float, Mode, Random>::compute_violated_constraints(
  std::vector<int>& container) const
{
    typename sparse_matrix<int>::const_row_iterator it, et;

    container.clear();

    for (int k = 0; k != m; ++k) {
        std::tie(it, et) = ap.row(k);
        int v = 0;

        for (; it != et; ++it)
            v += x[it->column];

        if (b[k] != v)
            container.emplace_back(k);
    }

    return length(container);
}

template<typename Float, typename Mode, typename Random>
double
solver_equalities_01coeff<Float, Mode, Random>::results(
  const std::unique_ptr<Float[]>& original_costs,
  const double cost_constant) const
{
    bx_expects(is_valid_solution());

    auto value = static_cast<double>(cost_constant);

    for (int i{ 0 }, ei{ n }; i != ei; ++i)
        value += static_cast<double>(original_costs[i] * x[i]);

    return value;
}

//
// Decrease influence of local preferences. 0 will completely reset the
// preference values for the current row. > 0 will keep former decision
// in mind.
//
template<typename Float, typename Mode, typename Random>
void
solver_equalities_01coeff<Float, Mode, Random>::decrease_preference(
  sparse_matrix<int>::row_iterator begin,
  sparse_matrix<int>::row_iterator end,
  Float theta) noexcept
{
    for (; begin != end; ++begin)
        P[begin->value] *= theta;
}

//
// Compute the reduced costs and return the size of the newly R vector.
//
template<typename Float, typename Mode, typename Random>
int
solver_equalities_01coeff<Float, Mode, Random>::compute_reduced_costs(
  sparse_matrix<int>::row_iterator begin,
  sparse_matrix<int>::row_iterator end) noexcept
{
    int r_size = 0;

    for (; begin != end; ++begin) {
        Float sum_a_pi = 0;
        Float sum_a_p = 0;

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

template<typename Float, typename Mode, typename Random>
int
solver_equalities_01coeff<Float, Mode, Random>::select_variables_equality(
  const int r_size,
  int bk)
{
    bk = std::min(bk, r_size);

    return bk - 1;
}

//
// The bkmin and bkmax constraint bounds are not equal and can be
// assigned to -infinity or +infinity. We have to scan the r vector and
// search a value j such as b(0, k) <= Sum A(k, R[j]) < b(1, k).
//
template<typename Float, typename Mode, typename Random>
void
solver_equalities_01coeff<Float, Mode, Random>::affect_variables(
  sparse_matrix<int>::row_iterator it,
  int k,
  int selected,
  int r_size,
  const Float kappa,
  const Float delta) noexcept
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
                  static_cast<Float>(2.0));

        Float d = delta + ((kappa / (static_cast<Float>(1.0) - kappa)) *
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

template<typename Float, typename Mode, typename Random>
void
solver_equalities_01coeff<Float, Mode, Random>::push_and_compute_update_row(
  std::vector<int>::iterator first,
  std::vector<int>::iterator last,
  Float kappa,
  Float delta,
  Float theta,
  Float objective_amplifier)
{
    for (; first != last; ++first) {
        auto k = *first;

        sparse_matrix<int>::row_iterator it, et;
        std::tie(it, et) = ap.row(k);

        decrease_preference(it, et, theta);
        const int r_size = compute_reduced_costs(it, et);

        if (objective_amplifier)
            for (int i = 0; i != r_size; ++i)
                R[i].value += objective_amplifier * c[R[i].id];

        calculator_sort(R.get(), R.get() + r_size, rng, Mode());
        int selected = select_variables_equality(r_size, b[k]);
        affect_variables(it, k, selected, r_size, kappa, delta);
    }
}

template<typename Float, typename Mode, typename Random>
void
solver_equalities_01coeff<Float, Mode, Random>::compute_update_row(
  std::vector<int>::iterator first,
  std::vector<int>::iterator last,
  Float kappa,
  Float delta,
  Float theta)
{
    for (; first != last; ++first) {
        auto k = *first;

        sparse_matrix<int>::row_iterator it, et;
        std::tie(it, et) = ap.row(k);

        decrease_preference(it, et, theta);
        const int r_size = compute_reduced_costs(it, et);
        calculator_sort(R.get(), R.get() + r_size, rng, Mode());
        int selected = select_variables_equality(r_size, b[k]);
        affect_variables(it, k, selected, r_size, kappa, delta);
    }
}

#if 0
extern template struct solver_equalities_01coeff<float,
                                                 minimize_tag,
                                                 std::default_random_engine>;
extern template struct solver_equalities_01coeff<double,
                                                 minimize_tag,
                                                 std::default_random_engine>;
extern template struct solver_equalities_01coeff<long double,
                                                 minimize_tag,
                                                 std::default_random_engine>;
extern template struct solver_equalities_01coeff<float,
                                                 maximize_tag,
                                                 std::default_random_engine>;
extern template struct solver_equalities_01coeff<double,
                                                 maximize_tag,
                                                 std::default_random_engine>;
extern template struct solver_equalities_01coeff<long double,
                                                 maximize_tag,
                                                 std::default_random_engine>;
#endif

} // namespace itm
} // namespace baryonyx

#endif
