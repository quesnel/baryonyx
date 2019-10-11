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

#include "bit-array.hpp"
#include "itm-common.hpp"
#include "itm-optimizer-common.hpp"
#include "itm-solver-common.hpp"

namespace baryonyx {
namespace itm {

template<typename Float, typename Mode, typename Cost>
struct solver_inequalities_101coeff
{
    using mode_type = Mode;
    using float_type = Float;
    using cost_type = Cost;

    random_engine& rng;

    struct rc_data
    {
        Float value;
        int f;
        int id;

        constexpr bool is_negative_factor() const noexcept
        {
            return f < 0;
        }

        int factor() const noexcept
        {
            return f;
        }
    };

    struct rc_size
    {
        int r_size;
        int c_size;
    };

    struct bound_factor
    {
        int min;
        int max;
        int negative_factor;
    };

    sparse_matrix<int> ap;
    std::unique_ptr<Float[]> P;
    value_bit_array<-1, 1> A;
    std::unique_ptr<rc_data[]> R;
    std::unique_ptr<bound_factor[]> b;
    std::unique_ptr<Float[]> pi;

    const cost_type& c;
    int m;
    int n;

    solver_inequalities_101coeff(random_engine& rng_,
                                 int m_,
                                 int n_,
                                 const cost_type& c_,
                                 const std::vector<merged_constraint>& csts)
      : rng(rng_)
      , ap(csts, m_, n_)
      , P(std::make_unique<Float[]>(ap.size()))
      , A(static_cast<int>(ap.size()))
      , R(std::make_unique<rc_data[]>(compute_reduced_costs_vector_size(csts)))
      , b(std::make_unique<bound_factor[]>(m_))
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

                if (cst.factor > 0)
                    A.set(id++);
                else
                    A.unset(id++);

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

            b[i].negative_factor = lower;

            bx_ensures(b[i].min <= b[i].max);
            bx_ensures(upper + lower == length(csts[i].elements));
        }
    }

    void reset() noexcept
    {
        std::fill_n(P.get(), ap.length(), Float{ 0 });
        std::fill_n(pi.get(), m, Float{ 0 });
    }

    int factor(int value) const noexcept
    {
        return A.get(value);
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
            ret += pi[ht->row];

        return ret;
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
    template<typename Xtype>
    rc_size compute_reduced_costs(sparse_matrix<int>::row_iterator begin,
                                  sparse_matrix<int>::row_iterator end,
                                  const Xtype& x) noexcept
    {
        int r_size = 0;
        int c_size = 0;

        for (; begin != end; ++begin) {
            Float sum_a_pi = 0;
            Float sum_a_p = 0;

            auto ht = ap.column(begin->column);

            for (; std::get<0>(ht) != std::get<1>(ht); ++std::get<0>(ht)) {
                auto a = static_cast<Float>(A.get(std::get<0>(ht)->value));

                sum_a_pi += a * pi[std::get<0>(ht)->row];
                sum_a_p += a * P[std::get<0>(ht)->value];
            }

            R[r_size].id = r_size;
            R[r_size].value = c(begin->column, x) - sum_a_pi - sum_a_p;
            R[r_size].f = A.get(begin->value);

            if (R[r_size].is_negative_factor()) {
                R[r_size].value = -R[r_size].value;
                ++c_size;
            }

            ++r_size;
        }

        return { r_size, c_size };
    }

    int select_variables(const rc_size& sizes, int bkmin, int bkmax)
    {
        if (bkmin == bkmax)
            return std::min(bkmin + sizes.c_size, sizes.r_size) - 1;

        bkmin += sizes.c_size;
        bkmax = std::min(bkmax + sizes.c_size, sizes.r_size);

        for (int i = bkmin; i <= bkmax; ++i)
            if (stop_iterating<Mode>(R[i].value, rng))
                return i - 1;

        return bkmax - 1;
    }

    template<typename Xtype, typename Iterator>
    bool push_and_compute_update_row(Xtype& x,
                                     Iterator first,
                                     Iterator last,
                                     Float kappa,
                                     Float delta,
                                     Float theta,
                                     Float obj_amp)
    {
        auto at_least_one_pi_changed{ false };

        for (; first != last; ++first) {
            auto k = constraint(first);

            bx_expects(k < m);

            const auto it = ap.row(k);
            decrease_preference(std::get<0>(it), std::get<1>(it), theta);

            const auto sizes =
              compute_reduced_costs(std::get<0>(it), std::get<1>(it), x);

            for (int i = 0; i != sizes.r_size; ++i)
                R[i].value +=
                  obj_amp * c((std::get<0>(it) + R[i].id)->column, x);

            calculator_sort<Mode>(R.get(), R.get() + sizes.r_size, rng);
            int selected = select_variables(sizes, b[k].min, b[k].max);

            auto pi_change = affect(*this,
                                    x,
                                    std::get<0>(it),
                                    k,
                                    selected,
                                    sizes.r_size,
                                    kappa,
                                    delta);

            at_least_one_pi_changed = at_least_one_pi_changed || pi_change;
        }

        return at_least_one_pi_changed;
    }

    template<typename Xtype, typename Iterator>
    bool compute_update_row(Xtype& x,
                            Iterator first,
                            Iterator last,
                            Float kappa,
                            Float delta,
                            Float theta)
    {
        auto at_least_one_pi_changed{ false };

        for (; first != last; ++first) {
            auto k = constraint(first);

            const auto it = ap.row(k);
            decrease_preference(std::get<0>(it), std::get<1>(it), theta);

            const auto sizes =
              compute_reduced_costs(std::get<0>(it), std::get<1>(it), x);

            calculator_sort<Mode>(R.get(), R.get() + sizes.r_size, rng);
            int selected = select_variables(sizes, b[k].min, b[k].max);

            auto pi_change = affect(*this,
                                    x,
                                    std::get<0>(it),
                                    k,
                                    selected,
                                    sizes.r_size,
                                    kappa,
                                    delta);

            at_least_one_pi_changed = at_least_one_pi_changed || pi_change;
        }

        return at_least_one_pi_changed;
    }
};

template<typename Solver,
         typename Float,
         typename Mode,
         typename Order,
         typename Cost>
static result
solve_or_optimize(const context_ptr& ctx,
                  const problem& pb,
                  bool is_optimization)
{
    return is_optimization
             ? optimize_problem<Solver, Float, Mode, Order, Cost>(ctx, pb)
             : solve_problem<Solver, Float, Mode, Order, Cost>(ctx, pb);
}

template<typename Float, typename Mode, typename Cost>
static result
select_order(const context_ptr& ctx, const problem& pb, bool is_optimization)
{
    switch (ctx->parameters.order) {
    case solver_parameters::constraint_order::none:
        return solve_or_optimize<
          solver_inequalities_101coeff<Float, Mode, Cost>,
          Float,
          Mode,
          constraint_sel<Float, 0>,
          Cost>(ctx, pb, is_optimization);
    case solver_parameters::constraint_order::reversing:
        return solve_or_optimize<
          solver_inequalities_101coeff<Float, Mode, Cost>,
          Float,
          Mode,
          constraint_sel<Float, 1>,
          Cost>(ctx, pb, is_optimization);
    case solver_parameters::constraint_order::random_sorting:
        return solve_or_optimize<
          solver_inequalities_101coeff<Float, Mode, Cost>,
          Float,
          Mode,
          constraint_sel<Float, 2>,
          Cost>(ctx, pb, is_optimization);
    case solver_parameters::constraint_order::infeasibility_decr:
        return solve_or_optimize<
          solver_inequalities_101coeff<Float, Mode, Cost>,
          Float,
          Mode,
          constraint_sel<Float, 3>,
          Cost>(ctx, pb, is_optimization);
    case solver_parameters::constraint_order::infeasibility_incr:
        return solve_or_optimize<
          solver_inequalities_101coeff<Float, Mode, Cost>,
          Float,
          Mode,
          constraint_sel<Float, 4>,
          Cost>(ctx, pb, is_optimization);
    case solver_parameters::constraint_order::lagrangian_decr:
        return solve_or_optimize<
          solver_inequalities_101coeff<Float, Mode, Cost>,
          Float,
          Mode,
          constraint_sel<Float, 5>,
          Cost>(ctx, pb, is_optimization);
    case solver_parameters::constraint_order::lagrangian_incr:
        return solve_or_optimize<
          solver_inequalities_101coeff<Float, Mode, Cost>,
          Float,
          Mode,
          constraint_sel<Float, 6>,
          Cost>(ctx, pb, is_optimization);
    case solver_parameters::constraint_order::pi_sign_change:
        return solve_or_optimize<
          solver_inequalities_101coeff<Float, Mode, Cost>,
          Float,
          Mode,
          constraint_sel<Float, 7>,
          Cost>(ctx, pb, is_optimization);
    default:
        bx_reach();
    }
}

template<typename Float, typename Mode>
static result
select_cost(const context_ptr& ctx, const problem& pb, bool is_optimization)
{
    return pb.objective.qelements.empty()
             ? select_order<Float,
                            Mode,
                            baryonyx::itm::default_cost_type<Float>>(
                 ctx, pb, is_optimization)
             : select_order<Float,
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
solve_inequalities_101(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "  - solve_inequalities_101\n");
    return select_float(ctx, pb, false);
}

result
optimize_inequalities_101(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "  - optimize_inequalities_101\n");
    return select_float(ctx, pb, true);
}

} // namespace itm
} // namespace baryonyx
