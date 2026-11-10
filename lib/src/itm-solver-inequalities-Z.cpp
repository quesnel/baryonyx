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

// Rows are handled for y = 1 - x on negative coefficients (Bastert et al.
// 2010, Alg. 3): weights are |a|, reduced costs of negative coefficients are
// negated and the bounds are shifted, so every row reads
// lo <= sum |a| y <= hi with R holding the weights |a|. Reduced costs use the
// signed coefficients.
template<typename Mode, typename Cost, bool debug>
struct solver_inequalities_Zcoeff : debug_logger<debug>
{
    using logger = debug_logger<debug>;
    using mode_type = Mode;
    using cost_type = Cost;

    static inline const int maximum_factor_exhaustive_solver = 10;

    // q = sign * r makes the attractive variables have q < 0 in both modes.
    static constexpr real sign =
      std::is_same_v<Mode, minimize_tag> ? One : -One;

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

    struct bound_factor
    {
        int min;
        int max;
        int negative_factor;
    };

    enum class subsolver_type : char
    {
        linear,
        exhaustive,
        branch_and_bound,
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
        std::size_t z_variables_max = 0;
        std::size_t z_constraint_exhaustive = 0;
        std::vector<function_element> y_elements;

        int id = 0;
        for (int i = 0, e = length(csts); i != e; ++i) {
            int lower = 0, upper = 0;
            Z[i] = subsolver_type::linear;

            for (const auto& cst : csts[i].elements) {
                bx_ensures(cst.factor);
                A[id++] = cst.factor;

                if (cst.factor > 0)
                    upper += cst.factor;
                else
                    lower += cst.factor;

                if (cst.factor < -1 || cst.factor > 1)
                    Z[i] = subsolver_type::branch_and_bound;
            }

            const auto row_size = csts[i].elements.size();
            if (Z[i] == subsolver_type::branch_and_bound &&
                row_size < maximum_factor_exhaustive_solver)
                Z[i] = subsolver_type::exhaustive;

            b[i].min = std::max(lower, csts[i].min);
            b[i].max = std::min(upper, csts[i].max);
            b[i].negative_factor = -lower;

            z_variables_max = std::max(z_variables_max, row_size);

            if (Z[i] == subsolver_type::exhaustive) {
                y_elements.clear();
                for (const auto& cst : csts[i].elements)
                    y_elements.emplace_back(std::abs(cst.factor),
                                            cst.variable_index);

                ++z_constraint_exhaustive;
                ex.build_constraints(i,
                                     y_elements,
                                     b[i].min + b[i].negative_factor,
                                     b[i].max + b[i].negative_factor);
            }

            bx_ensures(b[i].min <= b[i].max);
        }

        if (z_constraint_exhaustive > 0)
            ex.reserve(z_variables_max, z_constraint_exhaustive);

        bb.reserve(z_variables_max);
    }

    void reset() noexcept
    {
        std::fill_n(P.get(), ap.length(), Zero);
        std::fill_n(pi.get(), m, Zero);
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

        for (auto [ht, hend] = ap.column(variable); ht != hend; ++ht)
            ret += static_cast<real>(ht->factor) * pi[ht->row];

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
    real reduced_cost(int variable, const Xtype& x) const noexcept
    {
        real sum_a_pi = 0;
        real sum_a_p = 0;

        for (auto [ht, hte] = ap.column(variable); ht != hte; ++ht) {
            const auto a = static_cast<real>(ht->factor);
            sum_a_pi += a * pi[ht->row];
            sum_a_p += a * P[ht->value];
        }

        return c(variable, x) - sum_a_pi - sum_a_p;
    }

    template<typename Xtype>
    int compute_reduced_costs(sparse_matrix<int>::row_iterator begin,
                              sparse_matrix<int>::row_iterator end,
                              const Xtype& x,
                              real obj_amp) noexcept
    {
        int r_size = 0;

        for (; begin != end; ++begin, ++r_size) {
            auto value = reduced_cost(begin->column, x);
            if (obj_amp != Zero)
                value += obj_amp * c(begin->column, x);

            R[r_size].id = r_size;
            R[r_size].f = std::abs(A[begin->value]);
            R[r_size].value = A[begin->value] < 0 ? -value : value;
        }

        return r_size;
    }

    static real weight(const rc_data& r) noexcept
    {
        return static_cast<real>(r.f);
    }

    static real ratio(const rc_data& r) noexcept
    {
        return sign * r.value / weight(r);
    }

    // Greedy scan by increasing reduced cost per unit of weight (Bastert et
    // al. 2010, §3.2 and Alg. 4): items are taken until the lower bound is
    // met, then while they stay attractive and fit. Returns the number of
    // selected items or -1 when the scan cannot satisfy the row.
    int scan(int r_size, int lo, int hi, real phi) noexcept
    {
        std::sort(R.get(), R.get() + r_size, [](const auto& l, const auto& r) {
            return ratio(l) < ratio(r);
        });

        int total = 0;
        int count = 0;
        while (count != r_size) {
            const int w = R[count].f;
            if (total + w > hi)
                break;
            if (total >= lo && ratio(R[count]) + phi >= Zero)
                break;

            total += w;
            ++count;
        }

        return lo <= total && total <= hi ? count : -1;
    }

    // Returns the number of selected items, moved to the front of R, or -1
    // when the row cannot be satisfied.
    int select(int k, int r_size, bool& knapsack) noexcept
    {
        const int lo = b[k].min + b[k].negative_factor;
        const int hi = b[k].max + b[k].negative_factor;

        knapsack = false;
        if (const int count = scan(r_size, lo, hi, sign * pi[k]); count >= 0)
            return count;

        knapsack = true;
        int selected = -1;

        switch (Z[k]) {
        case subsolver_type::exhaustive:
            std::sort(
              R.get(), R.get() + r_size, [](const auto& l, const auto& r) {
                  return l.id < r.id;
              });
            selected = ex.solve(k, R, r_size);
            break;
        case subsolver_type::branch_and_bound:
            selected = bb.solve(R, r_size, lo, hi);
            if (!bb.has_solution)
                return -1;
            break;
        case subsolver_type::linear:
            return -1;
        }

        return std::min(selected + 1, r_size);
    }

    // pi[k] is set in the middle of the multipliers separating the selected
    // items from the others by reduced cost per unit of weight, with the
    // dual sign of the bounds that can be tight. Preferences follow Bastert
    // et al. 2010, Alg. 4, whose repair is applied after the greedy scan too:
    // with weights above one the scan alone leaves the solver cycling.
    template<typename Xtype, typename Iterator>
    bool affect(Xtype& x,
                Iterator it,
                int k,
                int r_size,
                int count,
                bool knapsack,
                real kappa,
                real delta)
    {
        constexpr real inf = std::numeric_limits<real>::infinity();

        const real old_pi = pi[k];
        const real old_phi = sign * old_pi;

        int total = 0;
        real lo = -inf;
        real hi = +inf;
        for (int i = 0; i != r_size; ++i) {
            const real t = ratio(R[i]) + old_phi;
            if (i < count) {
                total += R[i].f;
                lo = std::max(lo, t);
            } else {
                hi = std::min(hi, t);
            }
        }

        const real gap =
          !knapsack && count > 0 && count < r_size ? hi - lo : One;

        const int lower = b[k].min + b[k].negative_factor;
        const int upper = b[k].max + b[k].negative_factor;
        if (total > lower)
            hi = std::min(hi, Zero);
        if (total < upper)
            lo = std::max(lo, Zero);

        real phi = Zero;
        if (std::isfinite(lo) && std::isfinite(hi))
            phi = (lo + hi) / Two;
        else if (std::isfinite(lo))
            phi = lo;
        else if (std::isfinite(hi))
            phi = hi;

        pi[k] = sign * phi;

        const real d = sign * (delta + (kappa / (One - kappa)) * gap);

        for (int i = 0; i != r_size; ++i) {
            const auto var = it + R[i].id;
            const bool selected = i < count;

            if (selected != (A[var->value] < 0))
                x.set(var->column);
            else
                x.unset(var->column);

            P[var->value] += selected ? d : -d;
        }

        for (int i = 0; i != r_size; ++i) {
            const auto var = it + R[i].id;
            const real r_x = reduced_cost(var->column, x);
            const real q = sign * (A[var->value] < 0 ? -r_x : r_x);
            const real w = weight(R[i]);

            if (i < count && q >= Zero)
                P[var->value] += sign * (q + delta) / w;
            else if (i >= count && q <= Zero)
                P[var->value] += sign * (q - delta) / w;
        }

        bx_expects(is_valid_constraint(*this, k, x));

        return is_signbit_change(old_pi, pi[k]);
    }

    template<typename Xtype, typename Iterator>
    bool update_rows(Xtype& x,
                     Iterator first,
                     Iterator last,
                     real kappa,
                     real delta,
                     real theta,
                     real obj_amp)
    {
        bool pi_changed = false;

        for (; first != last; ++first) {
            const auto k = constraint(first);
            const auto [row_begin, row_end] = ap.row(k);

            decrease_preference(row_begin, row_end, theta);
            const int r_size =
              compute_reduced_costs(row_begin, row_end, x, obj_amp);

            bool knapsack = false;
            const int count = select(k, r_size, knapsack);
            if (count < 0)
                continue;

            pi_changed =
              affect(x, row_begin, k, r_size, count, knapsack, kappa, delta) ||
              pi_changed;
        }

        return pi_changed;
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
        return update_rows(x, first, last, kappa, delta, theta, obj_amp);
    }

    template<typename Xtype, typename Iterator>
    bool compute_update_row(Xtype& x,
                            Iterator first,
                            Iterator last,
                            real kappa,
                            real delta,
                            real theta)
    {
        return update_rows(x, first, last, kappa, delta, theta, Zero);
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
