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

#include "itm-common.hpp"
#include "itm-optimizer-common.hpp"
#include "itm-solver-common.hpp"

#include <algorithm>

namespace baryonyx {
namespace itm {

template<typename Float, typename Mode, typename Cost>
struct solver_random_inequalities_101coeff
{
    using mode_type = Mode;
    using float_type = Float;
    using cost_type = Cost;

    random_engine& rng;

    struct rc_data
    {
        Float value;
        int id;
        int a;

        constexpr bool is_negative() const
        {
            return a < 0;
        }

        constexpr int factor() const noexcept
        {
            return a;
        }
    };

    struct bound_factor
    {
        int min;
        int max;
        int negative_factor;
    };

    sparse_matrix<int> ap;
    std::unique_ptr<int[]> A;
    std::unique_ptr<rc_data[]> R;
    std::unique_ptr<bound_factor[]> b;
    std::unique_ptr<Float[]> pi;
    std::unique_ptr<Float[]> P;
    const cost_type& c;
    std::bernoulli_distribution dist;
    int m;
    int n;

    solver_random_inequalities_101coeff(
      random_engine& rng_,
      int m_,
      int n_,
      const cost_type& c_,
      const std::vector<merged_constraint>& csts)
      : rng(rng_)
      , ap(csts, m_, n_)
      , A(std::make_unique<int[]>(ap.size()))
      , R(std::make_unique<rc_data[]>(compute_reduced_costs_vector_size(csts)))
      , b(std::make_unique<bound_factor[]>(m_))
      , pi(std::make_unique<Float[]>(m_))
      , P(std::make_unique<Float[]>(ap.size()))
      , c(c_)
      , dist(0.5)
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

            bx_ensures(b[i].min <= b[i].max);
        }
    }

    void reset() const noexcept
    {
        std::fill_n(P.get(), ap.length(), Float{ 0 });
        std::fill_n(pi.get(), m, Float{ 0 });
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

    void decrease_preference(sparse_matrix<int>::row_iterator begin,
                             sparse_matrix<int>::row_iterator end,
                             Float theta) noexcept
    {
        for (; begin != end; ++begin)
            P[begin->value] *= theta;
    }

    template<typename Xtype, typename Iterator>
    bool push_and_compute_update_row(Xtype& x,
                                     Iterator first,
                                     Iterator last,
                                     Float kappa,
                                     Float delta,
                                     Float theta,
                                     Float /*obj_amp*/)
    {
        return compute_update_row(x, first, last, kappa, delta, theta);
    }

    template<typename Xtype, typename Iterator>
    bool compute_update_row(Xtype& x,
                            Iterator first,
                            Iterator last,
                            Float kappa,
                            Float delta,
                            Float theta)
    {
        for (; first != last; ++first) {
            auto k = constraint(first);
            bx_expects(k < m);

            auto [row_begin, row_end] = ap.row(k);
            decrease_preference(row_begin, row_end, theta);

            int r_size = 0;
            for (auto it = row_begin; it != row_end; ++it, ++r_size) {
                Float sum_a_p = 0;

                auto [col_begin, col_end] = ap.column(it->column);
                for (; col_begin != col_end; ++col_begin)
                    sum_a_p += static_cast<Float>(A[col_begin->value]) *
                               P[col_begin->value];

                R[r_size].id = r_size;
                R[r_size].a = A[it->value];
                R[r_size].value = c(row_begin->column, x) - sum_a_p;
            }

            calculator_sort<Mode>(R.get(), R.get() + r_size, rng);

            int value = 0;
            bool valid = b[k].min <= 0;
            int i = -1;

            // if (!valid) {
            //     do {
            //         ++i;
            //         value += R[i].a;
            //         valid = b[k].min <= value;
            //     } while (!valid && i + 1 < r_size);
            // }

            // valid = b[k].min <= value && value <= b[k].max;
            // while (i + 1 < r_size && valid) {
            //     ++i;
            //     value += R[i].a;
            //     valid = b[k].min <= value && value <= b[k].max;
            // }

            // std::uniform_real_distribution<Float> real_dist(Float{ 0 },
            //                                                 Float{ 1 });

            // auto selected = i;
            // auto d = delta;

            // if (selected < 0) {
            //     d += (kappa / (one - kappa)) * (R[0].value / two);
            //     d *= real_dist(rng);

            //     for (i = 0; i != r_size; ++i) {
            //         auto var = row_begin + R[i].id;

            //         if (R[i].is_negative()) {
            //             x.set(var->column, 1);
            //             P[var->value] += d;
            //         } else {
            //             x.set(var->column, 0);
            //             P[var->value] -= d;
            //         }
            //     }
            // } else if (selected + 1 >= r_size) {
            //     d += (kappa / (one - kappa)) * (R[selected].value * middle);
            //     d *= real_dist(rng);

            //     for (i = 0; i != r_size; ++i) {
            //         auto var = row_begin + R[i].id;

            //         if (R[i].is_negative()) {
            //             x.set(var->column, 0);
            //             P[var->value] -= d;
            //         } else {
            //             x.set(var->column, 1);
            //             P[var->value] += d;
            //         }
            //     }
            // } else {
            //     d += (kappa / (one - kappa)) *
            //          (R[selected + 1].value - R[selected].value);
            //     d *= real_dist(rng);

            //     for (i = 0; i <= selected; ++i) {
            //         auto var = row_begin + R[i].id;

            //         if (R[i].is_negative()) {
            //             x.set(var->column, 0);
            //             P[var->value] -= d;
            //         } else {
            //             x.set(var->column, 1);
            //             P[var->value] += d;
            //         }
            //     }

            //     for (; i != r_size; ++i) {
            //         auto var = row_begin + R[i].id;

            //         if (R[i].is_negative()) {
            //             x.set(var->column, 1);
            //             P[var->value] += d;
            //         } else {
            //             x.set(var->column, 0);
            //             P[var->value] -= d;
            //         }
            //     }

            constexpr Float one{ 1 };
            constexpr Float two{ 2 };

            if (!valid) {
                do {
                    ++i;
                    value += R[i].a;
                    valid = b[k].min <= value;
                    auto var = row_begin + R[i].id;
                    x.set(var->column);
                    P[var->value] +=
                      delta + (kappa / (one - kappa)) * (R[i].value / two);
                } while (!valid && i + 1 < r_size);
            }

            valid = b[k].min <= value && value <= b[k].max;
            while (i + 1 < r_size && valid) {
                ++i;
                value += R[i].a;
                valid = b[k].min <= value && value <= b[k].max;
                auto var = row_begin + R[i].id;

                if (valid)
                    valid = stop_iterating<Mode>(R[i].value, rng);

                if (valid) {
                    x.set(var->column);
                    P[var->value] +=
                      delta + (kappa / (one - kappa)) * (R[i].value / two);
                } else {
                    x.unset(var->column);
                    P[var->value] -=
                      delta + (kappa / (one - kappa)) * (R[i].value / two);
                }
            }

            while (i + 1 < r_size) {
                ++i;
                auto var = row_begin + R[i].id;
                x.unset(var->column);
                P[var->value] -=
                  delta + (kappa / (one - kappa)) * (R[i].value / two);
            }

            bx_expects(is_valid_constraint(*this, k, x));
        }

        return false;
    }
};

template<typename Float, typename Mode, typename Cost>
static result
solve_or_optimize(const context_ptr& ctx,
                  const problem& pb,
                  bool is_optimization)
{
    using Solver = solver_random_inequalities_101coeff<Float, Mode, Cost>;

    return is_optimization
             ? optimize_problem<Solver, Float, Mode, Cost>(ctx, pb)
             : solve_problem<Solver, Float, Mode, Cost>(ctx, pb);
}

template<typename Float, typename Mode>
static result
select_cost(const context_ptr& ctx, const problem& pb, bool is_optimization)
{
    return pb.objective.qelements.empty()
             ? solve_or_optimize<Float,
                                 Mode,
                                 baryonyx::itm::default_cost_type<Float>>(
                 ctx, pb, is_optimization)
             : solve_or_optimize<Float,
                                 Mode,
                                 baryonyx::itm::quadratic_cost_type<Float>>(
                 ctx, pb, is_optimization);
}

template<typename Float>
static result
select_mode(const context_ptr& ctx, const problem& pb, bool is_optimization)
{
    const auto m = static_cast<int>(pb.type);

    return m == 0 ? select_cost<Float, mode_sel<0>>(ctx, pb, is_optimization)
                  : select_cost<Float, mode_sel<1>>(ctx, pb, is_optimization);
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
solve_random_equalities_01(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "  - random::solver_equalities_01coeff\n");
    return select_float(ctx, pb, false);
}

result
optimize_random_equalities_01(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "  - random::optimizer_equalities_01coeff\n");
    return select_float(ctx, pb, false);
}

result
solve_random_inequalities_01(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "  - random::solver_inequalities_01coeff\n");
    return select_float(ctx, pb, true);
}

result
optimize_random_inequalities_01(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "  - random::optimizer_inequalities_01coeff\n");
    return select_float(ctx, pb, true);
}

result
solve_random_equalities_101(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "  - random::solver_equalities_101coeff\n");
    return select_float(ctx, pb, false);
}

result
optimize_random_equalities_101(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "  - random::optimizer_equalities_101coeff\n");
    return select_float(ctx, pb, false);
}

result
solve_random_inequalities_101(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "  - random::solver_inequalities_101coeff\n");
    return select_float(ctx, pb, true);
}

result
optimize_random_inequalities_101(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "  - random::optimizer_inequalities_101coeff\n");
    return select_float(ctx, pb, true);
}

} // namespace itm
} // namespace baryonyx
