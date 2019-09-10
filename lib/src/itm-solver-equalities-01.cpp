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

namespace baryonyx {
namespace itm {

template<typename Float, typename Mode, typename Random, typename Cost>
struct solver_equalities_01coeff
{
    using mode_type = Mode;
    using float_type = Float;
    using cost_type = Cost;

    Random& rng;

    struct rc_data
    {
        Float value; ///< Reduced cost value.
        int id;      ///< Index in ap.row() vector.

        constexpr bool is_negative_factor() const noexcept
        {
            return false;
        }

        constexpr int factor() const noexcept
        {
            return 1;
        }
    };

    using bound_factor = int;

    sparse_matrix<int> ap;
    std::unique_ptr<Float[]> P;
    std::unique_ptr<rc_data[]> R;
    std::unique_ptr<bound_factor[]> b;
    std::unique_ptr<Float[]> pi;

    const cost_type& c;
    int m;
    int n;

    solver_equalities_01coeff(Random& rng_,
                              int m_,
                              int n_,
                              const cost_type& c_,
                              const std::vector<merged_constraint>& csts)
      : rng(rng_)
      , ap(csts, m_, n_)
      , P(std::make_unique<Float[]>(ap.size()))
      , R(std::make_unique<rc_data[]>(compute_reduced_costs_vector_size(csts)))
      , b(std::make_unique<bound_factor[]>(m_))
      , pi(std::make_unique<Float[]>(m_))
      , c(c_)
      , m(m_)
      , n(n_)
    {
        for (int i = 0; i != m; ++i) {
#if !defined(BARYONYX_FULL_OPTIMIZATION)
            // mscv 15.9.6 fail to build this line:
            // for ([[maybe_unused]] const auto& cst : csts[i].elements)
            //    bx_ensures(cst.factor == 1);
            for (const auto& cst : csts[i].elements)
                bx_ensures(cst.factor == 1);
#endif

            bx_ensures(csts[i].min == csts[i].max);

            b[i] = csts[i].min;
        }
    }

    void reset() noexcept
    {
        std::fill_n(P.get(), ap.length(), Float{ 0 });
        std::fill_n(pi.get(), m, Float{ 0 });
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

    Float compute_sum_A_pi(int variable) const
    {
        Float ret{ 0 };

        sparse_matrix<int>::const_col_iterator ht, hend;
        std::tie(ht, hend) = ap.column(variable);

        for (; ht != hend; ++ht)
            ret += pi[ht->row];

        return ret;
    }

    void decrease_preference(sparse_matrix<int>::row_iterator begin,
                             sparse_matrix<int>::row_iterator end,
                             Float theta) noexcept
    {
        for (; begin != end; ++begin)
            P[begin->value] *= theta;
    }

    template<typename Xtype>
    int compute_reduced_costs(sparse_matrix<int>::row_iterator begin,
                              sparse_matrix<int>::row_iterator end,
                              const Xtype& x) noexcept
    {
        int r_size = 0;

        for (; begin != end; ++begin) {
            Float sum_a_pi = 0;
            Float sum_a_p = 0;

            auto ht = ap.column(begin->column);

            for (; std::get<0>(ht) != std::get<1>(ht); ++std::get<0>(ht)) {
                sum_a_pi += pi[std::get<0>(ht)->row];
                sum_a_p += P[std::get<0>(ht)->value];
            }

            R[r_size].id = r_size;
            R[r_size].value = c(begin->column, x) - sum_a_pi - sum_a_p;

            ++r_size;
        }

        return r_size;
    }

    int select_variables(const int r_size, const int bk)
    {
        // We need to affect all the bk variables even if the reduced cost sign
        // change.
        // For example: v + x + y + z = 2
        // bk = 2, r_size = 4
        // This function returns 1 i.e, we select the 0 and 1 variables.

        return std::min(bk, r_size) - 1;
    }

    template<typename Xtype, typename Iterator>
    bool push_and_compute_update_row(Xtype& x,
                                     Iterator first,
                                     Iterator last,
                                     Float kappa,
                                     Float delta,
                                     Float theta,
                                     Float objective_amplifier)
    {
        auto at_least_one_pi_changed{ false };

        for (; first != last; ++first) {
            auto k = constraint(first);

            const auto it = ap.row(k);
            decrease_preference(std::get<0>(it), std::get<1>(it), theta);

            const auto r_size =
              compute_reduced_costs(std::get<0>(it), std::get<1>(it), x);

            // Before sort and select variables, we apply the push method: for
            // each reduces cost, we had the cost multiply with an objective
            // amplifier.

            for (int i = 0; i != r_size; ++i)
                R[i].value += objective_amplifier * c(R[i].id, x);

            calculator_sort<Mode>(R.get(), R.get() + r_size, rng);
            int selected = select_variables(r_size, b[k]);

            auto pi_change = affect(
              *this, x, std::get<0>(it), k, selected, r_size, kappa, delta);

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

            const auto r_size =
              compute_reduced_costs(std::get<0>(it), std::get<1>(it), x);

            calculator_sort<Mode>(R.get(), R.get() + r_size, rng);
            int selected = select_variables(r_size, b[k]);

            auto pi_change = affect(
              *this, x, std::get<0>(it), k, selected, r_size, kappa, delta);

            at_least_one_pi_changed = at_least_one_pi_changed || pi_change;
        }

        return at_least_one_pi_changed;
    }
};

template<typename Solver,
         typename Float,
         typename Mode,
         typename Order,
         typename Random,
         typename Cost>
static result
solve_or_optimize(const context_ptr& ctx,
                  const problem& pb,
                  bool is_optimization)
{
    return is_optimization
             ? optimize_problem<Solver, Float, Mode, Order, Random, Cost>(ctx,
                                                                          pb)
             : solve_problem<Solver, Float, Mode, Order, Random, Cost>(ctx,
                                                                       pb);
}

template<typename Float, typename Mode, typename Random, typename Cost>
static result
select_order(const context_ptr& ctx, const problem& pb, bool is_optimization)
{
    switch (ctx->parameters.order) {
    case solver_parameters::constraint_order::none:
        return solve_or_optimize<
          solver_equalities_01coeff<Float, Mode, Random, Cost>,
          Float,
          Mode,
          constraint_sel<Float, Random, 0>,
          Random,
          Cost>(ctx, pb, is_optimization);
    case solver_parameters::constraint_order::reversing:
        return solve_or_optimize<
          solver_equalities_01coeff<Float, Mode, Random, Cost>,
          Float,
          Mode,
          constraint_sel<Float, Random, 1>,
          Random,
          Cost>(ctx, pb, is_optimization);
    case solver_parameters::constraint_order::random_sorting:
        return solve_or_optimize<
          solver_equalities_01coeff<Float, Mode, Random, Cost>,
          Float,
          Mode,
          constraint_sel<Float, Random, 2>,
          Random,
          Cost>(ctx, pb, is_optimization);
    case solver_parameters::constraint_order::infeasibility_decr:
        return solve_or_optimize<
          solver_equalities_01coeff<Float, Mode, Random, Cost>,
          Float,
          Mode,
          constraint_sel<Float, Random, 3>,
          Random,
          Cost>(ctx, pb, is_optimization);
    case solver_parameters::constraint_order::infeasibility_incr:
        return solve_or_optimize<
          solver_equalities_01coeff<Float, Mode, Random, Cost>,
          Float,
          Mode,
          constraint_sel<Float, Random, 4>,
          Random,
          Cost>(ctx, pb, is_optimization);
    case solver_parameters::constraint_order::lagrangian_decr:
        return solve_or_optimize<
          solver_equalities_01coeff<Float, Mode, Random, Cost>,
          Float,
          Mode,
          constraint_sel<Float, Random, 5>,
          Random,
          Cost>(ctx, pb, is_optimization);
    case solver_parameters::constraint_order::lagrangian_incr:
        return solve_or_optimize<
          solver_equalities_01coeff<Float, Mode, Random, Cost>,
          Float,
          Mode,
          constraint_sel<Float, Random, 6>,
          Random,
          Cost>(ctx, pb, is_optimization);
    case solver_parameters::constraint_order::pi_sign_change:
        return solve_or_optimize<
          solver_equalities_01coeff<Float, Mode, Random, Cost>,
          Float,
          Mode,
          constraint_sel<Float, Random, 7>,
          Random,
          Cost>(ctx, pb, is_optimization);
    default:
        bx_reach();
    }
}

template<typename Float, typename Mode, typename Random>
static result
select_cost(const context_ptr& ctx, const problem& pb, bool is_optimization)
{
    return pb.objective.qelements.empty()
             ? select_order<Float,
                            Mode,
                            Random,
                            baryonyx::itm::default_cost_type<Float>>(
                 ctx, pb, is_optimization)
             : select_order<Float,
                            Mode,
                            Random,
                            baryonyx::itm::quadratic_cost_type<Float>>(
                 ctx, pb, is_optimization);
}

template<typename Float, typename Mode>
static result
select_random(const context_ptr& ctx, const problem& pb, bool is_optimization)
{
    return select_cost<Float, Mode, std::default_random_engine>(
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
solve_equalities_01(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "  - solve_equalities_01\n");
    return select_float(ctx, pb, false);
}

result
optimize_equalities_01(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "  - optimize_equalities_01\n");
    return select_float(ctx, pb, true);
}

} // namespace itm
} // namespace baryonyx
