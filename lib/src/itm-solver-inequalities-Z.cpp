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

#include "branch-and-bound-solver.hpp"
#include "exhaustive-solver.hpp"
#include "itm-common.hpp"
#include "itm-optimizer-common.hpp"
#include "itm-solver-common.hpp"

namespace baryonyx {
namespace itm {

template<typename Mode, typename Cost, bool debug>
struct solver_inequalities_Zcoeff : debug_logger<debug>
{
    using logger = debug_logger<debug>;
    using mode_type = Mode;
    using cost_type = Cost;

    static inline const int maximum_factor_exhaustive_solver = 10;

    random_engine& rng;

    struct rc_data
    {
        real value;
        int id;
        int f;

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

    enum class subsolver_type : char
    {
        branch_and_bound,
        exhaustive,
        unconstrained_branch_and_bound,
        unconstrained_exhaustive,
        linear,
    };

    sparse_matrix<int> ap;
    std::unique_ptr<real[]> P;
    std::unique_ptr<int[]> A;
    std::unique_ptr<rc_data[]> R;
    std::unique_ptr<subsolver_type[]> Z;
    std::unique_ptr<bound_factor[]> b;
    std::unique_ptr<real[]> pi;

    branch_and_bound_solver<Mode> bb;
    exhaustive_solver<Mode> ex;

    const cost_type& c;
    int m;
    int n;

    solver_inequalities_Zcoeff(random_engine& rng_,
                               int m_,
                               int n_,
                               const cost_type& c_,
                               const std::vector<merged_constraint>& csts)
      : logger("solver_inequalities_Zcoeff")
      , rng(rng_)
      , ap(csts, m_, n_)
      , P(std::make_unique<real[]>(ap.size()))
      , A(std::make_unique<int[]>(ap.size()))
      , R(std::make_unique<rc_data[]>(compute_reduced_costs_vector_size(csts)))
      , Z(std::make_unique<subsolver_type[]>(m_))
      , b(std::make_unique<bound_factor[]>(m_))
      , pi(std::make_unique<real[]>(m_))
      , c(c_)
      , m(m_)
      , n(n_)
    {
        // Count the maximum function elements in the constraint where a least
        // one coefficient is in Z. To be use with the branch-and-bound and
        // exhaustive solvers.
        std::size_t z_variables_max = 0;

        // Count the maximum constraint number where at leat one coefficient is
        // in Z. To be used with the exhaustive solver.
        std::size_t z_constraint_exhaustive = 0;

        int id = 0;
        for (int i = 0, e = length(csts); i != e; ++i) {
            int lower = 0, upper = 0;
            std::size_t local_z_variables_max = 0;

            Z[i] =
              subsolver_type::linear; // Default, Z solver use the classical
                                      // Bastert selection for 101 problem.

            for (const auto& cst : csts[i].elements) {
                bx_ensures(cst.factor);
                A[id++] = cst.factor;
                ++local_z_variables_max;

                if (cst.factor > 0)
                    upper += cst.factor;
                else
                    lower += cst.factor;

                if (cst.factor < -1 || cst.factor > 1) {
                    Z[i] = subsolver_type::branch_and_bound;
                }
            }

            if (Z[i] == subsolver_type::branch_and_bound &&
                local_z_variables_max < maximum_factor_exhaustive_solver)
                Z[i] = subsolver_type::exhaustive;

            if (csts[i].min == csts[i].max) {
                b[i].min = csts[i].min;
                b[i].max = csts[i].max;
            } else {
                if (lower >= csts[i].min || upper <= csts[i].max) {
                    if (Z[i] == subsolver_type::branch_and_bound)
                        Z[i] = subsolver_type::unconstrained_branch_and_bound;
                    else if (Z[i] == subsolver_type::exhaustive)
                        Z[i] = subsolver_type::unconstrained_exhaustive;
                }

                b[i].min = std::max(lower, csts[i].min);
                b[i].max = std::min(upper, csts[i].max);
            }

            z_variables_max = std::max(z_variables_max, local_z_variables_max);

            if (Z[i] == subsolver_type::exhaustive ||
                Z[i] == subsolver_type::unconstrained_exhaustive) {
                z_constraint_exhaustive++;
                ex.build_constraints(i, csts[i].elements, b[i].min, b[i].max);
            }

            logger::log("Is Z: {} ({}) with {} <= {}\n",
                        static_cast<int>(Z[i]),
                        local_z_variables_max,
                        b[i].min,
                        b[i].max);

            bx_ensures(b[i].min <= b[i].max);
        }

        if (z_constraint_exhaustive > 0)
            ex.reserve(z_variables_max, z_constraint_exhaustive);

        bb.reserve(z_variables_max);
    }

    void reset() noexcept
    {
        logger::log("reset-solver\n");

        std::fill_n(P.get(), ap.length(), real{ 0 });
        std::fill_n(pi.get(), m, real{ 0 });
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

    real compute_sum_A_pi(int variable) const
    {
        real ret{ 0 };

        sparse_matrix<int>::const_col_iterator ht, hend;
        std::tie(ht, hend) = ap.column(variable);

        for (; ht != hend; ++ht)
            ret += std::abs(static_cast<real>(A[ht->value])) * pi[ht->row];

        return ret;
    }

    void decrease_preference(sparse_matrix<int>::row_iterator begin,
                             sparse_matrix<int>::row_iterator end,
                             real theta) noexcept
    {
        for (; begin != end; ++begin)
            P[begin->value] *= theta;
    }

    //
    // Compute the reduced costs and return the size of the newly R vector.
    //
    template<typename Xtype>
    int compute_reduced_costs(sparse_matrix<int>::row_iterator begin,
                              sparse_matrix<int>::row_iterator end,
                              const Xtype& x) noexcept
    {
        // to_log(
        //  debug_os, "  compute-reduced-cost {}\n", std::distance(begin,
        //  end));

        int r_size = 0;

        for (; begin != end; ++begin) {
            real sum_a_pi_p = 0;

            for (auto [first, last] = ap.column(begin->column); first != last;
                 ++first) {
                auto a = std::abs(static_cast<real>(A[first->value]));
                sum_a_pi_p += a * (pi[first->row] + P[first->value]);
            }

            R[r_size].id = r_size;
            R[r_size].f = A[begin->value];
            R[r_size].value = c(begin->column, x) - sum_a_pi_p;

            // to_log(debug_os,
            //       4u,
            //       "Compute: {} = {} - {} - {}\n",
            //       r_size,
            //       c(begin->column, x),
            //       sum_pi,
            //       sum_p);

            // to_log(debug_os, 4u, "{}x{}\n", R[r_size].f, R[r_size].value);

            ++r_size;
        }

        // to_log(debug_os, "\n");

        return r_size;
    }

    template<typename Xtype>
    real local_compute_reduced_cost(int variable, const Xtype& x) noexcept
    {
        real sum_a_pi_p = 0;

        for (auto [ht, hte] = ap.column(variable); ht != hte; ++ht) {
            auto a = std::abs(static_cast<real>(A[ht->value]));
            sum_a_pi_p += a * (pi[ht->row] + P[ht->value]);
        }

        return c(variable, x) - sum_a_pi_p;
    }

    int select_variables_101(const int r_size, int bkmin, int bkmax)
    {
        int sum = 0;
        int best = -2;

        for (int i = -1; i < r_size; ++i) {
            if (bkmin <= sum && sum <= bkmax)
                best = i;

            if (best != -2 && i - 1 < r_size &&
                stop_iterating<Mode>(R[i + 1].value))
                break;

            sum += R[i + 1].f;
        }

        return best;
    }

    int select_variables(const int r_size, int bkmin, int bkmax)
    {
        int sum = 0;
        int best = -2;

        for (int i = -1; i < r_size; ++i) {
            if (bkmin <= sum && sum <= bkmax)
                best = i;

            if (best != -2 && i - 1 < r_size &&
                stop_iterating<Mode>(R[i + 1].value))
                break;

            sum += R[i + 1].f;
        }

        return best;
    }

    template<typename Xtype, typename Iterator>
    bool local_affect(Xtype& x,
                      Iterator it,
                      int k,
                      int selected,
                      int r_size,
                      const real kappa,
                      const real delta)
    {
        constexpr real one{ 1 };
        constexpr real two{ 2 };
        constexpr real middle{ (two + one) / two };

        const auto old_pi = pi[k];

        auto d = (kappa / (one - kappa)) + delta;

        if (selected < 0) {
            pi[k] += R[0].value / two;

            logger::log("    selected: {}/{} ({})/2 = pi {}\n",
                        selected,
                        r_size,
                        R[0].value,
                        pi[k]);

            for (int i = 0; i != r_size; ++i) {
                auto var = it + R[i].id;

                x.unset(var->column);
                P[var->value] -= d;

                auto repair = local_compute_reduced_cost(var->column, x);
                if (repair <= 0)
                    P[var->value] = P[var->value] + repair - d;
            }
        } else if (selected + 1 >= r_size) {
            pi[k] += R[selected].value * middle;

            logger::log("    selected: {}/{} ({})/2 = pi {}\n",
                        selected,
                        r_size,
                        R[selected].value,
                        pi[k]);

            for (int i = 0; i != r_size; ++i) {
                auto var = it + R[i].id;

                x.set(var->column);
                P[var->value] += d;

                auto repair = local_compute_reduced_cost(var->column, x);
                if (repair >= 0)
                    P[var->value] = P[var->value] - repair + d;
            }
        } else {
            pi[k] += ((R[selected].value + R[selected + 1].value) / two);

            logger::log("    selected: {}/{} ({}x{})/2 = pi {}\n",
                        selected,
                        r_size,
                        R[selected].value,
                        R[selected + 1].value,
                        pi[k]);

            int i = 0;
            for (; i <= selected; ++i) {
                auto var = it + R[i].id;

                x.set(var->column);
                P[var->value] += d;

                auto repair = local_compute_reduced_cost(var->column, x);
                if (repair >= 0)
                    P[var->value] = P[var->value] - repair + d;
            }

            for (; i != r_size; ++i) {
                auto var = it + R[i].id;

                x.unset(var->column);
                P[var->value] -= d;

                auto repair = local_compute_reduced_cost(var->column, x);
                if (repair <= 0)
                    P[var->value] = P[var->value] + repair - d;
            }
        }

        // TODO job: develops is_valid_constraint for all the solvers
        bx_expects(is_valid_constraint(*this, k, x));

        return is_signbit_change(old_pi, pi[k]);
    }

    template<typename Xtype, typename Iterator>
    bool push_and_compute_update_row(Xtype& x,
                                     Iterator first,
                                     Iterator last,
                                     real kappa,
                                     real delta,
                                     real theta,
                                     real obj_amp)
    {
        auto at_least_one_pi_changed{ false };

        logger::log("push-update-row {} {} {}\n", kappa, delta, theta);

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
                R[i].value +=
                  obj_amp * c((std::get<0>(it) + R[i].id)->column, x);

            real pi_change;
            int selected;

            switch (Z[k]) {
            case subsolver_type::branch_and_bound:
                calculator_sort<Mode>(R.get(), R.get() + r_size, rng);
                selected = bb.solve(R, r_size, b[k].min, b[k].max);
                break;
            case subsolver_type::exhaustive:
                selected = ex.solve(k, R, r_size);
                break;
            case subsolver_type::unconstrained_branch_and_bound:
                calculator_sort<Mode>(R.get(), R.get() + r_size, rng);
                selected = select_variables(r_size, b[k].min, b[k].max);
                if (selected == -2)
                    selected = bb.solve(R, r_size, b[k].min, b[k].max);
                break;
            case subsolver_type::unconstrained_exhaustive:
                calculator_sort<Mode>(R.get(), R.get() + r_size, rng);
                selected = select_variables(r_size, b[k].min, b[k].max);
                if (selected == -2)
                    selected = ex.solve(k, R, r_size);
                break;
            case subsolver_type::linear:
            default:
                calculator_sort<Mode>(R.get(), R.get() + r_size, rng);
                selected = select_variables_101(r_size, b[k].min, b[k].max);
                break;
            }

            pi_change = local_affect(
              x, std::get<0>(it), k, selected, r_size, kappa, delta);
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

            real pi_change;
            int selected = -2;

            switch (Z[k]) {
            case subsolver_type::branch_and_bound:
                calculator_sort<Mode>(R.get(), R.get() + r_size, rng);
                selected = bb.solve(R, r_size, b[k].min, b[k].max);
                break;
            case subsolver_type::exhaustive:
                selected = ex.solve(k, R, r_size);
                break;
            case subsolver_type::unconstrained_branch_and_bound:
                calculator_sort<Mode>(R.get(), R.get() + r_size, rng);
                selected = select_variables(r_size, b[k].min, b[k].max);
                if (selected == -2)
                    selected = bb.solve(R, r_size, b[k].min, b[k].max);
                break;
            case subsolver_type::unconstrained_exhaustive:
                calculator_sort<Mode>(R.get(), R.get() + r_size, rng);
                selected = select_variables(r_size, b[k].min, b[k].max);
                if (selected == -2)
                    selected = ex.solve(k, R, r_size);
                break;
            case subsolver_type::linear:
                calculator_sort<Mode>(R.get(), R.get() + r_size, rng);
                selected = select_variables_101(r_size, b[k].min, b[k].max);
                break;
            }

            pi_change = local_affect(
              x, std::get<0>(it), k, selected, r_size, kappa, delta);
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
        using Solver = solver_inequalities_Zcoeff<Mode, Cost, true>;

        return is_optimization
                 ? optimize_problem<Solver, Mode, Cost>(ctx, pb)
                 : solve_problem<Solver, Mode, Cost>(ctx, pb);
    } else {
        using Solver = solver_inequalities_Zcoeff<Mode, Cost, false>;

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
             ? solve_or_optimize<Mode,
                                 baryonyx::itm::default_cost_type>(
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
solve_inequalities_Z(const context& ctx, const problem& pb)
{
    info(ctx, "  - solve_inequalities_Z\n");
    return select_mode(ctx, pb, false);
}

result
optimize_inequalities_Z(const context& ctx, const problem& pb)
{
    info(ctx, "  - optimize_inequalities_Z\n");
    return select_mode(ctx, pb, true);
}

} // namespace itm
} // namespace baryonyx
