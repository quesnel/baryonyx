/* Copyright (C) 2016-2021 INRAE
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

template<typename Mode, typename Cost, bool debug>
struct solver_inequalities_01coeff : debug_logger<debug>
{
    using logger = debug_logger<debug>;
    using mode_type = Mode;
    using cost_type = Cost;

    random_engine& rng;

    struct rc_data
    {
        real value;
        int id;

        constexpr bool is_negative_factor() const noexcept
        {
            return false;
        }

        constexpr int factor() const noexcept
        {
            return 1;
        }
    };

    struct bound_factor
    {
        int min;
        int max;
    };

    sparse_matrix<int> ap;
    std::unique_ptr<real[]> P;
    std::unique_ptr<rc_data[]> R;
    std::unique_ptr<bound_factor[]> b;
    std::unique_ptr<real[]> pi;

    const cost_type& c;
    int m;
    int n;

    solver_inequalities_01coeff(random_engine& rng_,
                                int m_,
                                int n_,
                                const cost_type& c_,
                                const std::vector<merged_constraint>& csts)
      : logger("solver_inequalities_01coeff")
      , rng(rng_)
      , ap(csts, m_, n_)
      , P(std::make_unique<real[]>(ap.size()))
      , R(std::make_unique<rc_data[]>(compute_reduced_costs_vector_size(csts)))
      , b(std::make_unique<bound_factor[]>(m_))
      , pi(std::make_unique<real[]>(m_))
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

            if (csts[i].min == csts[i].max) {
                b[i].min = csts[i].min;
                b[i].max = csts[i].max;
            } else {
                b[i].min = std::max(0, csts[i].min);
                b[i].max = std::min(length(csts[i].elements), csts[i].max);
            }

            bx_ensures(b[i].min <= b[i].max);
        }
    }

    void reset() noexcept
    {
        std::fill_n(P.get(), ap.length(), real{ 0 });
        std::fill_n(pi.get(), m, real{ 0 });
    }

    int factor(int /*value*/) const noexcept
    {
        return 1;
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

    real compute_sum_A_pi(int variable) const
    {
        real ret{ 0 };

        sparse_matrix<int>::const_col_iterator ht, hend;
        std::tie(ht, hend) = ap.column(variable);

        for (; ht != hend; ++ht)
            ret += pi[ht->row];

        return ret;
    }

    void decrease_preference(sparse_matrix<int>::row_iterator begin,
                             sparse_matrix<int>::row_iterator end,
                             real theta) noexcept
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
            real sum_a_pi = 0;
            real sum_a_p = 0;

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

    int select_variables(const int r_size, int bkmin, int bkmax)
    {
        if (bkmin == bkmax)
            return std::min(bkmin, r_size) - 1;

        bkmin = std::min(bkmin, r_size);
        bkmax = std::min(bkmax, r_size);

        for (int i = bkmin; i <= bkmax; ++i)
            if (stop_iterating<Mode>(R[i].value, rng))
                return i - 1;

        return bkmax - 1;
    }

    template<typename Xtype, typename Iterator>
    bool push_and_compute_update_row(Xtype& x,
                                     Iterator first,
                                     Iterator last,
                                     real kappa,
                                     real delta,
                                     real theta,
                                     real objective_amplifier)
    {
        auto at_least_one_pi_changed{ false };

        logger::log("push-update-row {} {} {}\n", kappa, delta, theta);

        for (; first != last; ++first) {
            auto k = constraint(first);

            bx_expects(k < m);

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
            int selected = select_variables(r_size, b[k].min, b[k].max);

            logger::log("constraints {}: {} <= ", k, b[k].min);

            for (int i = 0; i < r_size; ++i)
                logger::log("({} {}) ", R[i].value, R[i].id);

            logger::log("<= {} => Selected: {}\n", b[k].max, selected);

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
                            real kappa,
                            real delta,
                            real theta)
    {
        auto at_least_one_pi_changed{ false };

        logger::log("update-row {} {} {}\n", kappa, delta, theta);

        for (; first != last; ++first) {
            auto k = constraint(first);

            const auto it = ap.row(k);
            decrease_preference(std::get<0>(it), std::get<1>(it), theta);

            const auto r_size =
              compute_reduced_costs(std::get<0>(it), std::get<1>(it), x);

            calculator_sort<Mode>(R.get(), R.get() + r_size, rng);
            int selected = select_variables(r_size, b[k].min, b[k].max);

            logger::log("constraints {}: {} <= ", k, b[k].min);

            for (int i = 0; i < r_size; ++i)
                logger::log("({} {}) ", R[i].value, R[i].id);

            logger::log("<= {} => Selected: {}\n", b[k].max, selected);

            auto pi_change = affect(
              *this, x, std::get<0>(it), k, selected, r_size, kappa, delta);

            at_least_one_pi_changed = at_least_one_pi_changed || pi_change;
        }

        return at_least_one_pi_changed;
    }
};

template<typename Mode, typename Cost>
static result
solve_or_optimize(const context& ctx, const problem& pb, bool is_optimization)
{
    if (ctx.parameters.debug) {
        using Solver = solver_inequalities_01coeff<Mode, Cost, true>;

        return is_optimization
                 ? optimize_problem<Solver, Mode, Cost>(ctx, pb)
                 : solve_problem<Solver, Mode, Cost>(ctx, pb);
    } else {
        using Solver = solver_inequalities_01coeff<Mode, Cost, false>;

        return is_optimization
                 ? optimize_problem<Solver, Mode, Cost>(ctx, pb)
                 : solve_problem<Solver, Mode, Cost>(ctx, pb);
    }
}

template<typename Mode>
static result
select_cost(const context& ctx, const problem& pb, bool is_optimization)
{
    return pb.objective.qelements.empty()
             ? solve_or_optimize<Mode, baryonyx::itm::default_cost_type>(
                 ctx, pb, is_optimization)
             : solve_or_optimize<Mode,
                                 baryonyx::itm::quadratic_cost_type>(
                 ctx, pb, is_optimization);
}

static result
select_mode(const context& ctx, const problem& pb, bool is_optimization)
{
    const auto m = static_cast<int>(pb.type);

    return m == 0 ? select_cost<mode_sel<0>>(ctx, pb, is_optimization)
                  : select_cost<mode_sel<1>>(ctx, pb, is_optimization);
}

result
solve_inequalities_01(const context& ctx, const problem& pb)
{
    info(ctx, "  - solve_inequalities_01\n");
    return select_mode(ctx, pb, false);
}

result
optimize_inequalities_01(const context& ctx, const problem& pb)
{
    info(ctx, "  - optimize_inequalities_01\n");
    return select_mode(ctx, pb, true);
}

} // namespace itm
} // namespace baryonyx
