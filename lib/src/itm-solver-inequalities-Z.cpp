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

#include "itm-common.hpp"
#include "itm-optimizer-common.hpp"
#include "itm-solver-common.hpp"
#include "knapsack-dp-solver.hpp"

namespace baryonyx {
namespace itm {

template<typename Float, typename Mode, typename Random>
struct solver_inequalities_Zcoeff
{
    using mode_type = Mode;
    using float_type = Float;

    Random& rng;

    sparse_matrix<int> ap;
    std::unique_ptr<Float[]> P;
    std::unique_ptr<int[]> A;
    std::unique_ptr<r_data<Float>[]> R;
    fixed_array<fixed_array<c_data>> C;
    std::vector<bool> Z;

    std::unique_ptr<bound[]> b;
    std::unique_ptr<Float[]> pi;

    const std::unique_ptr<Float[]>& c;
    int m;
    int n;

    solver_inequalities_Zcoeff(Random& rng_,
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
      , C(m_)
      , Z(m_, false)
      , b(std::make_unique<bound[]>(m_))
      , pi(std::make_unique<Float[]>(m_))
      , c(c_)
      , m(m_)
      , n(n_)
    {
        int id = 0;
        for (int i = 0, e = length(csts); i != e; ++i) {
            int lower_size = 0, upper_size = 0;
            int lower = 0, upper = 0;

            for (const auto& cst : csts[i].elements) {
                bx_ensures(cst.factor);
                A[id++] = cst.factor;

                if (cst.factor > 0) {
                    upper_size++;
                    upper += cst.factor;
                } else {
                    lower_size++;
                    lower += cst.factor;
                }

                if (std::abs(cst.factor) > 1)
                    Z[i] = true;
            }

            if (csts[i].min == csts[i].max) {
                b[i].min = csts[i].min;
                b[i].max = csts[i].max;
            } else {
                b[i].min = std::max(-lower, csts[i].min);
                b[i].max = std::min(upper, csts[i].max);
            }

            if (lower_size > 0) {
                C[i] = fixed_array<c_data>(lower_size);

                int id_in_r = 0;
                int id_in_c = 0;

                typename sparse_matrix<int>::const_row_iterator it, et;
                std::tie(it, et) = ap.row(i);

                for (; it != et; ++it) {
                    if (A[it->value] < 0) {
                        C[i][id_in_c].id_r = id_in_r;
                        ++id_in_c;
                    }
                    ++id_in_r;
                }
            }
        }
    }

    int factor(int value) const noexcept
    {
        return A[value];
    }

    int bound_min(int constraint) const noexcept
    {
        return b[constraint].min;
    }

    int bound_max(int constraint) const noexcept
    {
        return b[constraint].max;
    }

    int bound_init(int constraint) const
    {
        return bound_init(constraint, Mode());
    }

    int bound_init(int constraint, minimize_tag) const
    {
        return b[constraint].min;
    }

    int bound_init(int constraint, maximize_tag) const
    {
        return b[constraint].max;
    }

    Float compute_sum_A_pi(int variable) const
    {
        Float ret{ 0 };

        sparse_matrix<int>::const_col_iterator ht, hend;
        std::tie(ht, hend) = ap.column(variable);

        for (; ht != hend; ++ht)
            ret += std::abs(static_cast<Float>(A[ht->value])) * pi[ht->row];

        return ret;
    }

    template<typename Xtype>
    bool is_valid_solution(const Xtype& x) const
    {
        for (int k = 0, ek = m; k != ek; ++k) {
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

    template<typename Xtype>
    int compute_violated_constraints(const Xtype& x,
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

    template<typename Xtype>
    double results(const Xtype& x,
                   const std::unique_ptr<Float[]>& original_costs,
                   const double cost_constant) const
    {
        bx_expects(is_valid_solution(x));

        auto value = static_cast<double>(cost_constant);

        for (int i{ 0 }, ei{ n }; i != ei; ++i)
            value += static_cast<double>(original_costs[i] * x[i]);

        return value;
    }

    template<typename Xtype>
    void compute_update_row_Z_eq(Xtype& x,
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
                R[i].value += objective_amplifier * c[(it + R[i].id)->column];

        //
        // Negate reduced costs and coefficients of these variables. We need to
        // parse the row Ak[i] because we need to use r[i] not available in C.
        //

        for (int i = 0; i != c_size; ++i) {
            R[ck[i].id_r].value = -R[ck[i].id_r].value;

            auto var = it + ck[i].id_r;

            P[var->value] = -P[var->value];
            bk_move += A[var->value];
        }

        bk += std::abs(bk_move);

        calculator_sort(R.get(), R.get() + r_size, rng, Mode());

        //
        //
        //

        int selected = select_variables_equality_Z(it, r_size, bk);

        affect_variables(x, it, k, selected, r_size, kappa, delta);

        //
        // Clean up: correct negated costs and adjust value of negated
        // variables.
        //

        for (int i = 0; i != c_size; ++i) {
            auto var = it + ck[i].id_r;

            P[var->value] = -P[var->value];
            x.invert(var->column);
        }
    }

    template<typename Xtype>
    void compute_update_row_Z_ineq(Xtype& x,
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
                R[i].value += objective_amplifier * c[(it + R[i].id)->column];

        //
        // Negate reduced costs and coefficients of these variables. We need to
        // parse the row Ak[i] because we need to use r[i] not available in C.
        //

        for (int i = 0; i != c_size; ++i) {
            R[ck[i].id_r].value = -R[ck[i].id_r].value;

            auto var = it + ck[i].id_r;

            P[var->value] = -P[var->value];
            bk_move += A[var->value];
        }

        bkmin += std::abs(bk_move);
        bkmax += std::abs(bk_move);

        calculator_sort(R.get(), R.get() + r_size, rng, Mode());

        //
        //
        //

        int selected = select_variables_inequality_Z(it, r_size, bkmin, bkmax);

        affect_variables(x, it, k, selected, r_size, kappa, delta);

        //
        // Clean up: correct negated costs and adjust value of negated
        // variables.
        //

        for (int i = 0; i != c_size; ++i) {
            auto var = it + ck[i].id_r;

            P[var->value] = -P[var->value];
            x.invert(var->column);
        }
    }

    template<typename Xtype>
    void compute_update_row_01_eq(Xtype& x,
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

    template<typename Xtype>
    void compute_update_row_01_ineq(Xtype& x,
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

    template<typename Xtype>
    void compute_update_row_101_eq(Xtype& x,
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

        const auto& ck = C[k];
        const int c_size = length(ck);
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

        for (int i = 0; i != c_size; ++i) {
            R[ck[i].id_r].value = -R[ck[i].id_r].value;

            auto var = it + ck[i].id_r;

            P[var->value] = -P[var->value];
        }

        bk += c_size;

        calculator_sort(R.get(), R.get() + r_size, rng, Mode());

        int selected = select_variables_equality(r_size, bk);

        affect_variables(x, it, k, selected, r_size, kappa, delta);

        //
        // Clean up: correct negated costs and adjust value of negated
        // variables.
        //

        for (int i = 0; i != c_size; ++i) {
            auto var = it + ck[i].id_r;

            P[var->value] = -P[var->value];
            x.invert(var->column);
        }
    }

    template<typename Xtype>
    void compute_update_row_101_ineq(Xtype& x,
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

        const auto& ck = C[k];
        const int c_size = (ck) ? length(ck) : 0;
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

        for (int i = 0; i != c_size; ++i) {
            R[ck[i].id_r].value = -R[ck[i].id_r].value;

            auto var = it + ck[i].id_r;

            P[var->value] = -P[var->value];
        }

        bkmin += c_size;
        bkmax += c_size;

        calculator_sort(R.get(), R.get() + r_size, rng, Mode());

        int selected = select_variables_inequality(r_size, bkmin, bkmax);

        affect_variables(x, it, k, selected, r_size, kappa, delta);

        //
        // Clean up: correct negated costs and adjust value of negated
        // variables.
        //

        for (int i = 0; i != c_size; ++i) {
            auto var = it + ck[i].id_r;

            P[var->value] = -P[var->value];
            x.invert(var->column);
        }
    }

    //
    // Decrease influence of local preferences. 0 will completely reset the
    // preference values for the current row. > 0 will keep former decision in
    // mind.
    //
    void decrease_preference(sparse_matrix<int>::row_iterator begin,
                             sparse_matrix<int>::row_iterator end,
                             Float theta) noexcept
    {
        for (; begin != end; ++begin)
            P[begin->value] *= theta;
    }

    //
    // Compute the reduced costs and return the size of the newly R vector.
    //
    int compute_reduced_costs(sparse_matrix<int>::row_iterator begin,
                              sparse_matrix<int>::row_iterator end)
    {
        int r_size = 0;

        for (; begin != end; ++begin) {
            Float sum_a_pi = 0;
            Float sum_a_p = 0;

            typename sparse_matrix<int>::const_col_iterator ht, hend;
            std::tie(ht, hend) = ap.column(begin->column);

            for (; ht != hend; ++ht) {
                auto f = A[ht->value];
                auto a = static_cast<Float>(f);

                sum_a_pi += a * pi[ht->row];
                sum_a_p += a * P[ht->value];
            }

            // R[r_size].id = begin->column;
            R[r_size].id = r_size;
            R[r_size].value = c[begin->column] - sum_a_pi - sum_a_p;
            ++r_size;
        }

        return r_size;
    }

    int select_variables_equality(const int r_size, int bk)
    {
        bk = std::min(bk, r_size);

        return bk - 1;
    }

    int select_variables_inequality(const int r_size, int bkmin, int bkmax)
    {
        bkmin = std::min(bkmin, r_size);
        bkmax = std::min(bkmax, r_size);

        for (int i = bkmin; i <= bkmax; ++i)
            if (stop_iterating(R[i].value, rng, Mode()))
                return i - 1;

        return bkmax - 1;
    }

    int select_variables_equality_Z(sparse_matrix<int>::row_iterator it,
                                    const int r_size,
                                    int bk)
    {
        int sum = A[(it + R[0].id)->value];
        if (sum == bk)
            return 0;

        for (int i = 1; i != r_size; ++i) {
            sum += A[(it + R[i].id)->value];

            if (sum == bk)
                return i;

            if (sum > bk)
                break;
        }

        //
        // If the previous selection failed, we try a dynamic programming to
        // solve this knapsack 01 problem.
        //

        return knapsack_dp_solver<Mode, Float>(A, R, it, r_size, bk);
    }

    int select_variables_inequality_Z(sparse_matrix<int>::row_iterator it,
                                      const int r_size,
                                      int bkmin,
                                      int bkmax)
    {
        int sum = A[(it + R[0].id)->value];
        if (bkmin <= sum && sum <= bkmax)
            return 0;

        for (int i = 1; i != r_size; ++i) {
            sum += A[(it + R[i].id)->value];

            if (bkmin <= sum && sum <= bkmax)
                return i;

            if (sum > bkmax)
                break;
        }

        //
        // If the previous selection failed, we try a dynamic programming to
        // solve this knapsack 01 problem.
        //

        return knapsack_dp_solver<Mode, Float>(A, R, it, r_size, bkmax);
    }

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

    template<typename Xtype, typename Iterator>
    void push_and_compute_update_row(Xtype& x,
                                     Iterator first,
                                     Iterator last,
                                     Float kappa,
                                     Float delta,
                                     Float theta,
                                     Float obj_amp)
    {
        for (; first != last; ++first) {
            auto k = *first;

            if (Z[k]) {
                if (b[k].min == b[k].max)
                    compute_update_row_Z_eq(
                      x, k, b[k].min, kappa, delta, theta, obj_amp);
                else
                    compute_update_row_Z_ineq(
                      x, k, b[k].min, b[k].max, kappa, delta, theta, obj_amp);
            } else if (!C[k]) {
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

    template<typename Xtype, typename Iterator>
    void compute_update_row(Xtype& x,
                            Iterator first,
                            Iterator last,
                            Float kappa,
                            Float delta,
                            Float theta)
    {
        for (; first != last; ++first) {
            auto k = *first;
            if (Z[k]) {
                if (b[k].min == b[k].max)
                    compute_update_row_Z_eq(x,
                                            k,
                                            b[k].min,
                                            kappa,
                                            delta,
                                            theta,
                                            static_cast<Float>(0));
                else
                    compute_update_row_Z_ineq(x,
                                              k,
                                              b[k].min,
                                              b[k].max,
                                              kappa,
                                              delta,
                                              theta,
                                              static_cast<Float>(0));
            } else if (!C[k]) {
                if (b[k].min == b[k].max)
                    compute_update_row_01_eq(x,
                                             k,
                                             b[k].min,
                                             kappa,
                                             delta,
                                             theta,
                                             static_cast<Float>(0));
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
                    compute_update_row_101_eq(x,
                                              k,
                                              b[k].min,
                                              kappa,
                                              delta,
                                              theta,
                                              static_cast<Float>(0));
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
};

template<typename Solver,
         typename Float,
         typename Mode,
         typename Order,
         typename Random>
static result
solve_or_optimize(const context_ptr& ctx,
                  const problem& pb,
                  bool is_optimization)
{
    return is_optimization
             ? optimize_problem<Solver, Float, Mode, Order, Random>(ctx, pb)
             : solve_problem<Solver, Float, Mode, Order, Random>(ctx, pb);
}

template<typename Float, typename Mode, typename Random>
static result
select_order(const context_ptr& ctx, const problem& pb, bool is_optimization)
{
    const auto c = static_cast<int>(ctx->parameters.order);

    if (c == 0)
        return solve_or_optimize<
          solver_inequalities_Zcoeff<Float, Mode, Random>,
          Float,
          Mode,
          constraint_sel<Float, Random, 0>,
          Random>(ctx, pb, is_optimization);
    else if (c == 1)
        return solve_or_optimize<
          solver_inequalities_Zcoeff<Float, Mode, Random>,
          Float,
          Mode,
          constraint_sel<Float, Random, 1>,
          Random>(ctx, pb, is_optimization);
    else if (c == 2)
        return solve_or_optimize<
          solver_inequalities_Zcoeff<Float, Mode, Random>,
          Float,
          Mode,
          constraint_sel<Float, Random, 2>,
          Random>(ctx, pb, is_optimization);
    else if (c == 3)
        return solve_or_optimize<
          solver_inequalities_Zcoeff<Float, Mode, Random>,
          Float,
          Mode,
          constraint_sel<Float, Random, 3>,
          Random>(ctx, pb, is_optimization);
    else
        return solve_or_optimize<
          solver_inequalities_Zcoeff<Float, Mode, Random>,
          Float,
          Mode,
          constraint_sel<Float, Random, 4>,
          Random>(ctx, pb, is_optimization);
}

template<typename Float, typename Mode>
static result
select_random(const context_ptr& ctx, const problem& pb, bool is_optimization)
{
    return select_order<Float, Mode, std::default_random_engine>(
      ctx, pb, is_optimization);
}

template<typename Float>
static result
select_mode(const context_ptr& ctx, const problem& pb, bool is_optimization)
{
    const auto m = static_cast<int>(pb.type);

    return m == 0
             ? select_random<Float, mode_sel<0>>(ctx, pb, is_optimization)
             : select_random<Float, mode_sel<1>>(ctx, pb, is_optimization);
}

static result
select_float(const context_ptr& ctx, const problem& pb, bool is_optimization)
{
    const auto f = static_cast<int>(ctx->parameters.float_type);

    if (f == 0)
        return select_mode<float_sel<0>>(ctx, pb, is_optimization);
    else if (f == 1)
        return select_mode<float_sel<1>>(ctx, pb, is_optimization);
    else
        return select_mode<float_sel<2>>(ctx, pb, is_optimization);
}

result
solve_inequalities_Z(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "  - solve_inequalities_Z\n");
    return select_float(ctx, pb, false);
}

result
optimize_inequalities_Z(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "  - optimize_inequalities_Z\n");
    return select_float(ctx, pb, true);
}

} // namespace itm
} // namespace baryonyx
