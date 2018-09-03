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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_SOLVER_INEQUALITIES_101_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_SOLVER_INEQUALITIES_101_HPP

#include "itm-common.hpp"
#include "sparse-vector.hpp"

namespace baryonyx {
namespace itm {

template<typename Float, typename Mode, typename Random>
struct solver_inequalities_101coeff
{
    using mode_type = Mode;
    using float_type = Float;

    Random& rng;

    sparse_matrix<int> ap;
    std::unique_ptr<Float[]> P;
    std::unique_ptr<int[]> A;
    std::unique_ptr<r_data<Float>[]> R;

    sparse_vector<c_data> C;

    std::unique_ptr<bound[]> b;
    std::unique_ptr<Float[]> pi;

    const std::unique_ptr<Float[]>& c;
    int m;
    int n;

    solver_inequalities_101coeff(Random& rng_,
                                 int m_,
                                 int n_,
                                 const std::unique_ptr<Float[]>& c_,
                                 const std::vector<merged_constraint>& csts);

    int factor(int value) const noexcept;

    int bound_min(int constraint) const noexcept;

    int bound_max(int constraint) const noexcept;

    int bound_init(int constraint) const;

    int bound_init(int constraint, minimize_tag) const;

    int bound_init(int constraint, maximize_tag) const;

    Float compute_sum_A_pi(int variable) const;

    template<typename Xtype>
    bool is_valid_solution(const Xtype& x) const;

    template<typename Xtype>
    int compute_violated_constraints(const Xtype& x,
                                     std::vector<int>& container) const;

    template<typename Xtype>
    double results(const Xtype& x,
                   const std::unique_ptr<Float[]>& original_costs,
                   const double cost_constant) const;

    template<typename Xtype>
    void compute_update_row_01_eq(Xtype& x,
                                  int k,
                                  int bk,
                                  Float kappa,
                                  Float delta,
                                  Float theta,
                                  Float objective_amplifier);

    template<typename Xtype>
    void compute_update_row_01_ineq(Xtype& x,
                                    int k,
                                    int bkmin,
                                    int bkmax,
                                    Float kappa,
                                    Float delta,
                                    Float theta,
                                    Float objective_amplifier);

    template<typename Xtype>
    void compute_update_row_101_eq(Xtype& x,
                                   int k,
                                   int bk,
                                   Float kappa,
                                   Float delta,
                                   Float theta,
                                   Float objective_amplifier);

    template<typename Xtype>
    void compute_update_row_101_ineq(Xtype& x,
                                     int k,
                                     int bkmin,
                                     int bkmax,
                                     Float kappa,
                                     Float delta,
                                     Float theta,
                                     Float objective_amplifier);

    //
    // Decrease influence of local preferences. 0 will completely reset the
    // preference values for the current row. > 0 will keep former decision in
    // mind.
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
    int select_variables_inequality(const int r_size, int bkmin, int bkmax);
    //
    // The bkmin and bkmax constraint bounds are not equal and can be assigned
    // to -infinity or +infinity. We have to scan the r vector and search a
    // value j such as b(0, k) <= Sum A(k, R[j]) < b(1, k).
    //
    template<typename Xtype>
    void affect_variables(Xtype& x,
                          sparse_matrix<int>::row_iterator it,
                          int k,
                          int selected,
                          int r_size,
                          const Float kappa,
                          const Float delta) noexcept;

    template<typename Xtype>
    void push_and_compute_update_row(Xtype& x,
                                     std::vector<int>::iterator first,
                                     std::vector<int>::iterator last,
                                     Float kappa,
                                     Float delta,
                                     Float theta,
                                     Float obj_amp);

    template<typename Xtype>
    void compute_update_row(Xtype& x,
                            std::vector<int>::iterator first,
                            std::vector<int>::iterator last,
                            Float kappa,
                            Float delta,
                            Float theta);
};

template<typename Float, typename Mode, typename Random>
solver_inequalities_101coeff<Float, Mode, Random>::
  solver_inequalities_101coeff(Random& rng_,
                               int m_,
                               int n_,
                               const std::unique_ptr<Float[]>& c_,
                               const std::vector<merged_constraint>& csts)
  : rng(rng_)
  , ap(csts, m_, n_)
  , P(std::make_unique<Float[]>(ap.size()))
  , A(std::make_unique<int[]>(ap.size()))
  , R(std::make_unique<r_data<Float>[]>(
      compute_reduced_costs_vector_size(csts)))
  , C(csts)
  , b(std::make_unique<bound[]>(m_))
  , pi(std::make_unique<Float[]>(m_))
  , c(c_)
  , m(m_)
  , n(n_)
{
    int id = 0;
    for (int i = 0, e = length(csts); i != e; ++i) {
        int lower = 0, upper = 0;

        for (const auto& cst : csts[i].elements) {
            bx_ensures(std::abs(cst.factor) == 1);
            A[id++] = cst.factor;

            if (cst.factor > 0)
                upper++;
            else
                lower++;
        }

        if (csts[i].min == csts[i].max) {
            b[i].min = csts[i].min;
            b[i].max = csts[i].max;
        } else {
            b[i].min = std::max(-lower, csts[i].min);
            b[i].max = std::min(upper, csts[i].max);
        }
    }
}

template<typename Float, typename Mode, typename Random>
int
solver_inequalities_101coeff<Float, Mode, Random>::factor(int value) const
  noexcept
{
    return A[value];
}

template<typename Float, typename Mode, typename Random>
int
solver_inequalities_101coeff<Float, Mode, Random>::bound_min(
  int constraint) const noexcept
{
    return b[constraint].min;
}

template<typename Float, typename Mode, typename Random>
int
solver_inequalities_101coeff<Float, Mode, Random>::bound_max(
  int constraint) const noexcept
{
    return b[constraint].max;
}

template<typename Float, typename Mode, typename Random>
int
solver_inequalities_101coeff<Float, Mode, Random>::bound_init(
  int constraint) const
{
    return bound_init(constraint, Mode());
}

template<typename Float, typename Mode, typename Random>
int
solver_inequalities_101coeff<Float, Mode, Random>::bound_init(
  int constraint,
  minimize_tag) const
{
    return b[constraint].min;
}

template<typename Float, typename Mode, typename Random>
int
solver_inequalities_101coeff<Float, Mode, Random>::bound_init(
  int constraint,
  maximize_tag) const
{
    return b[constraint].max;
}

template<typename Float, typename Mode, typename Random>
Float
solver_inequalities_101coeff<Float, Mode, Random>::compute_sum_A_pi(
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
template<typename Xtype>
bool
solver_inequalities_101coeff<Float, Mode, Random>::is_valid_solution(
  const Xtype& x) const
{
    for (int k = 0; k != m; ++k) {
        typename sparse_matrix<int>::const_row_iterator it, et;

        std::tie(it, et) = ap.row(k);
        int v = 0;

        for (; it != et; ++it)
            v += A[it->value] * x[it->column];

        if (!(b[k].min <= v && v <= b[k].max))
            return false;
    }

    return true;
}

template<typename Float, typename Mode, typename Random>
template<typename Xtype>
int
solver_inequalities_101coeff<Float, Mode, Random>::
  compute_violated_constraints(const Xtype& x,
                               std::vector<int>& container) const
{
    typename sparse_matrix<int>::const_row_iterator it, et;

    container.clear();

    for (int k = 0; k != m; ++k) {
        std::tie(it, et) = ap.row(k);
        int v = 0;

        for (; it != et; ++it)
            v += A[it->value] * x[it->column];

        if (!(b[k].min <= v && v <= b[k].max))
            container.emplace_back(k);
    }

    return length(container);
}

template<typename Float, typename Mode, typename Random>
template<typename Xtype>
double
solver_inequalities_101coeff<Float, Mode, Random>::results(
  const Xtype& x,
  const std::unique_ptr<Float[]>& original_costs,
  const double cost_constant) const
{
    bx_expects(is_valid_solution(x));

    auto value = static_cast<double>(cost_constant);

    for (int i{ 0 }, ei{ n }; i != ei; ++i)
        value += static_cast<double>(original_costs[i] * x[i]);

    return value;
}

template<typename Float, typename Mode, typename Random>
template<typename Xtype>
void
solver_inequalities_101coeff<Float, Mode, Random>::compute_update_row_01_eq(
  Xtype& x,
  int k,
  int bk,
  Float kappa,
  Float delta,
  Float theta,
  Float objective_amplifier)
{
    typename sparse_matrix<int>::row_iterator it, et;
    std::tie(it, et) = ap.row(k);

    decrease_preference(it, et, theta);

    const int r_size = compute_reduced_costs(it, et);

    //
    // Before sort and select variables, we apply the push method: for each
    // reduces cost, we had the cost multiply with an objective amplifier.
    //

    if (objective_amplifier)
        for (int i = 0; i != r_size; ++i)
            R[i].value += objective_amplifier * c[(it + R[i].id)->column];

    calculator_sort(R.get(), R.get() + r_size, rng, Mode());

    int selected = select_variables_equality(r_size, bk);

    affect_variables(x, it, k, selected, r_size, kappa, delta);
}

template<typename Float, typename Mode, typename Random>
template<typename Xtype>
void
solver_inequalities_101coeff<Float, Mode, Random>::compute_update_row_01_ineq(
  Xtype& x,
  int k,
  int bkmin,
  int bkmax,
  Float kappa,
  Float delta,
  Float theta,
  Float objective_amplifier)
{
    typename sparse_matrix<int>::row_iterator it, et;
    std::tie(it, et) = ap.row(k);

    decrease_preference(it, et, theta);

    const int r_size = compute_reduced_costs(it, et);

    //
    // Before sort and select variables, we apply the push method: for each
    // reduces cost, we had the cost multiply with an objective amplifier.
    //

    if (objective_amplifier)
        for (int i = 0; i != r_size; ++i)
            R[i].value += objective_amplifier * c[(it + R[i].id)->column];

    calculator_sort(R.get(), R.get() + r_size, rng, Mode());

    int selected = select_variables_inequality(r_size, bkmin, bkmax);

    affect_variables(x, it, k, selected, r_size, kappa, delta);
}

template<typename Float, typename Mode, typename Random>
template<typename Xtype>
void
solver_inequalities_101coeff<Float, Mode, Random>::compute_update_row_101_eq(
  Xtype& x,
  int k,
  int bk,
  Float kappa,
  Float delta,
  Float theta,
  Float objective_amplifier)
{
    typename sparse_matrix<int>::row_iterator it, et;
    std::tie(it, et) = ap.row(k);

    decrease_preference(it, et, theta);

    const int r_size = compute_reduced_costs(it, et);

    //
    // Before sort and select variables, we apply the push method: for each
    // reduces cost, we had the cost multiply with an objective amplifier.
    //

    if (objective_amplifier)
        for (int i = 0; i != r_size; ++i)
            R[i].value += objective_amplifier * c[(it + R[i].id)->column];

    //
    // Negate reduced costs and coefficients of these variables. We need to
    // parse the row Ak[i] because we need to use r[i] not available in C.
    //

    typename sparse_vector<c_data>::iterator c_begin, c_end;
    std::tie(c_begin, c_end) = C.range(k);

    for (auto c_it = c_begin; c_it != c_end; ++c_it) {
        R[c_it->id_r].value = -R[c_it->id_r].value;
        auto var = it + c_it->id_r;
        P[var->value] = -P[var->value];
    }

    bk += static_cast<int>(std::distance(c_begin, c_end));

    calculator_sort(R.get(), R.get() + r_size, rng, Mode());

    int selected = select_variables_equality(r_size, bk);

    affect_variables(x, it, k, selected, r_size, kappa, delta);

    //
    // Clean up: correct negated costs and adjust value of negated
    // variables.
    //

    for (auto c_it = c_begin; c_it != c_end; ++c_it) {
        auto var = it + c_it->id_r;
        P[var->value] = -P[var->value];
        x.invert(var->column);
    }
}

template<typename Float, typename Mode, typename Random>
template<typename Xtype>
void
solver_inequalities_101coeff<Float, Mode, Random>::compute_update_row_101_ineq(
  Xtype& x,
  int k,
  int bkmin,
  int bkmax,
  Float kappa,
  Float delta,
  Float theta,
  Float objective_amplifier)
{
    typename sparse_matrix<int>::row_iterator it, et;
    std::tie(it, et) = ap.row(k);

    decrease_preference(it, et, theta);

    typename sparse_vector<c_data>::iterator c_begin, c_end;
    std::tie(c_begin, c_end) = C.range(k);

    const int r_size = compute_reduced_costs(it, et);

    //
    // Before sort and select variables, we apply the push method: for each
    // reduces cost, we had the cost multiply with an objective amplifier.
    //

    if (objective_amplifier)
        for (int i = 0; i != r_size; ++i)
            R[i].value += objective_amplifier * c[(it + R[i].id)->column];

    //
    // Negate reduced costs and coefficients of these variables. We need to
    // parse the row Ak[i] because we need to use r[i] not available in C.
    //

    for (auto c_it = c_begin; c_it != c_end; ++c_it) {
        R[c_it->id_r].value = -R[c_it->id_r].value;
        auto var = it + c_it->id_r;
        P[var->value] = -P[var->value];
    }

    bkmin += static_cast<int>(std::distance(c_begin, c_end));
    bkmax += static_cast<int>(std::distance(c_begin, c_end));

    calculator_sort(R.get(), R.get() + r_size, rng, Mode());

    int selected = select_variables_inequality(r_size, bkmin, bkmax);

    affect_variables(x, it, k, selected, r_size, kappa, delta);

    //
    // Clean up: correct negated costs and adjust value of negated
    // variables.
    //

    for (auto c_it = c_begin; c_it != c_end; ++c_it) {
        auto var = it + c_it->id_r;
        P[var->value] = -P[var->value];
        x.invert(var->column);
    }
}

//
// Decrease influence of local preferences. 0 will completely reset the
// preference values for the current row. > 0 will keep former decision in
// mind.
//
template<typename Float, typename Mode, typename Random>
void
solver_inequalities_101coeff<Float, Mode, Random>::decrease_preference(
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
solver_inequalities_101coeff<Float, Mode, Random>::compute_reduced_costs(
  sparse_matrix<int>::row_iterator begin,
  sparse_matrix<int>::row_iterator end) noexcept
{
    int r_size = 0;

    for (; begin != end; ++begin) {
        Float sum_a_pi = 0;
        Float sum_a_p = 0;

        typename sparse_matrix<int>::const_col_iterator ht, hend;
        std::tie(ht, hend) = ap.column(begin->column);

        for (; ht != hend; ++ht) {
            auto a = static_cast<Float>(A[ht->value]);

            sum_a_pi += a * pi[ht->row];
            sum_a_p += a * P[ht->value];
        }

        R[r_size].id = r_size;
        R[r_size].value = c[begin->column] - sum_a_pi - sum_a_p;
        ++r_size;
    }

    return r_size;
}

template<typename Float, typename Mode, typename Random>
int
solver_inequalities_101coeff<Float, Mode, Random>::select_variables_equality(
  const int r_size,
  int bk)
{
    bk = std::min(bk, r_size);

    return bk - 1;
}

template<typename Float, typename Mode, typename Random>
int
solver_inequalities_101coeff<Float, Mode, Random>::select_variables_inequality(
  const int r_size,
  int bkmin,
  int bkmax)
{
    bkmin = std::min(bkmin, r_size);
    bkmax = std::min(bkmax, r_size);

    for (int i = bkmin; i <= bkmax; ++i)
        if (stop_iterating(R[i].value, rng, Mode()))
            return i - 1;

    return bkmax - 1;
}

//
// The bkmin and bkmax constraint bounds are not equal and can be assigned
// to -infinity or +infinity. We have to scan the r vector and search a
// value j such as b(0, k) <= Sum A(k, R[j]) < b(1, k).
//
template<typename Float, typename Mode, typename Random>
template<typename Xtype>
void
solver_inequalities_101coeff<Float, Mode, Random>::affect_variables(
  Xtype& x,
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

            x.set(var->column, false);
            P[var->value] -= delta;
        }
    } else if (selected + 1 >= r_size) {
        pi[k] += R[selected].value;

        for (int i = 0; i != r_size; ++i) {
            auto var = it + R[i].id;

            x.set(var->column, true);
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

            x.set(var->column, true);
            P[var->value] += d;
        }

        for (; i != r_size; ++i) {
            auto var = it + R[i].id;

            x.set(var->column, false);
            P[var->value] -= d;
        }
    }
}

template<typename Float, typename Mode, typename Random>
template<typename Xtype>
void
solver_inequalities_101coeff<Float, Mode, Random>::push_and_compute_update_row(
  Xtype& x,
  std::vector<int>::iterator first,
  std::vector<int>::iterator last,
  Float kappa,
  Float delta,
  Float theta,
  Float obj_amp)
{
    for (; first != last; ++first) {
        auto k = *first;

        if (C.empty(k)) {
            if (b[k].min == b[k].max)
                compute_update_row_01_eq(
                  x, k, b[k].min, kappa, delta, theta, obj_amp);
            else
                compute_update_row_01_ineq(
                  x, k, b[k].min, b[k].max, kappa, delta, theta, obj_amp);
        } else {
            if (b[k].min == b[k].max)
                compute_update_row_101_eq(
                  x, k, b[k].min, kappa, delta, theta, obj_amp);
            else
                compute_update_row_101_ineq(
                  x, k, b[k].min, b[k].max, kappa, delta, theta, obj_amp);
        }
    }
}

template<typename Float, typename Mode, typename Random>
template<typename Xtype>
void
solver_inequalities_101coeff<Float, Mode, Random>::compute_update_row(
  Xtype& x,
  std::vector<int>::iterator first,
  std::vector<int>::iterator last,
  Float kappa,
  Float delta,
  Float theta)
{
    for (; first != last; ++first) {
        auto k = *first;

        if (C.empty(k)) {
            if (b[k].min == b[k].max)
                compute_update_row_01_eq(
                  x, k, b[k].min, kappa, delta, theta, static_cast<Float>(0));
            else
                compute_update_row_01_ineq(x,
                                           k,
                                           b[k].min,
                                           b[k].max,
                                           kappa,
                                           delta,
                                           theta,
                                           static_cast<Float>(0));
        } else {
            if (b[k].min == b[k].max)
                compute_update_row_101_eq(
                  x, k, b[k].min, kappa, delta, theta, static_cast<Float>(0));
            else
                compute_update_row_101_ineq(x,
                                            k,
                                            b[k].min,
                                            b[k].max,
                                            kappa,
                                            delta,
                                            theta,
                                            static_cast<Float>(0));
        }
    }
}

#if 0
extern template struct solver_inequalities_101coeff<
  float,
  minimize_tag,
  std::default_random_engine>;
extern template struct solver_inequalities_101coeff<
  double,
  minimize_tag,
  std::default_random_engine>;
extern template struct solver_inequalities_101coeff<
  long double,
  minimize_tag,
  std::default_random_engine>;
extern template struct solver_inequalities_101coeff<
  float,
  maximize_tag,
  std::default_random_engine>;
extern template struct solver_inequalities_101coeff<
  double,
  maximize_tag,
  std::default_random_engine>;
extern template struct solver_inequalities_101coeff<
  long double,
  maximize_tag,
  std::default_random_engine>;
#endif

} // namespace itm
} // namespace baryonyx

#endif
