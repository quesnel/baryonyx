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

#include <algorithm>
#include <cmath>
#include <limits>

namespace baryonyx {
namespace itm {

// Rows are updated by batches of rows sharing no variable: inside a batch the
// row computations read nothing another row of the batch writes (pi, P and x
// of its own variables), so they can run concurrently without locks while
// producing exactly the result of the sequential Wedelin sweep.
template<typename Mode, typename Cost, bool debug>
struct solver_inequalities_101coeff_buffered : debug_logger<debug>
{
    using logger = debug_logger<debug>;
    using mode_type = Mode;
    using cost_type = Cost;

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

    struct row_update
    {
        real pi;
        real d;
        int k;
        int count;
    };

    // Reduced costs are handled as q = sign * r so that attractive variables
    // have q < 0 for both minimization and maximization.
    static constexpr real sign =
      std::is_same_v<Mode, minimize_tag> ? One : -One;

    sparse_matrix<> ap;
    std::unique_ptr<real[]> P;
    std::unique_ptr<std::int8_t[]> A;
    std::unique_ptr<rc_data[]> R;
    std::unique_ptr<bound_factor[]> b;
    std::unique_ptr<real[]> pi;
    std::vector<row_update> batch;
    std::vector<int> variable_batch;
    std::vector<int> pending;
    std::vector<int> deferred;
    int batch_id = 0;

    const cost_type& c;
    int m;
    int n;

    solver_inequalities_101coeff_buffered(
      random_engine& rng_,
      int m_,
      int n_,
      const cost_type& c_,
      const std::vector<merged_constraint>& csts)
      : logger("solver_inequalities_101coeff_buffered")
      , rng(rng_)
      , ap(csts, m_, n_)
      , P(std::make_unique<real[]>(ap.size()))
      , A(std::make_unique<std::int8_t[]>(ap.size()))
      , R(std::make_unique<rc_data[]>(ap.size()))
      , b(std::make_unique<bound_factor[]>(m_))
      , pi(std::make_unique<real[]>(m_))
      , variable_batch(static_cast<std::size_t>(n_), -1)
      , c(c_)
      , m(m_)
      , n(n_)
    {
        batch.reserve(static_cast<std::size_t>(m));
        pending.reserve(static_cast<std::size_t>(m));
        deferred.reserve(static_cast<std::size_t>(m));

        int id = 0;
        for (int i = 0, e = length(csts); i != e; ++i) {
            int lower = 0, upper = 0;

            for (const auto& cst : csts[i].elements) {
                bx_ensures(std::abs(cst.factor) == 1);
                A[id++] = static_cast<std::int8_t>(cst.factor);

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

    int bound_init(int constraint) const noexcept
    {
        return bound_init(constraint, Mode());
    }

    int bound_init(int constraint, minimize_tag) const noexcept
    {
        return b[constraint].min;
    }

    int bound_init(int constraint, maximize_tag) const noexcept
    {
        return b[constraint].max;
    }

    real compute_sum_A_pi(int variable) const noexcept
    {
        real ret{ 0 };

        auto h = ap.column(variable);

        for (; std::get<0>(h) != std::get<1>(h); ++std::get<0>(h))
            ret += pi[std::get<0>(h)->row];

        return ret;
    }

    // Reduced costs are stored for y = 1 - x on negative coefficients
    // (Bastert et al. 2010, Alg. 3) and the row multiplier pi[k] is placed in
    // the middle of the multipliers supporting the selection, with the dual
    // sign of the bounds that can be tight.
    template<typename Xtype>
    row_update compute_row(const Xtype& x,
                           int k,
                           real kappa,
                           real delta,
                           real theta,
                           real obj_amp) noexcept
    {
        constexpr real inf = std::numeric_limits<real>::infinity();

        auto [row_first, row_last] = ap.row(k);
        const int r_size = static_cast<int>(row_last - row_first);
        rc_data* r = &R[row_first->value];

        for (auto it = row_first; it != row_last; ++it)
            P[it->value] *= theta;

        int i = 0;
        for (auto it = row_first; it != row_last; ++it, ++i) {
            real sum_a_pi = Zero;
            real sum_a_p = Zero;

            for (auto [cb, ce] = ap.column(it->column); cb != ce; ++cb) {
                const auto a = static_cast<real>(A[cb->value]);
                sum_a_pi += a * pi[cb->row];
                sum_a_p += a * P[cb->value];
            }

            r[i].id = i;
            r[i].f = A[it->value];
            r[i].value = c(it->column, x) - sum_a_pi - sum_a_p;

            if (r[i].is_negative_factor())
                r[i].value = -r[i].value;

            if (obj_amp != Zero) {
                const auto push = obj_amp * c(it->column, x);
                r[i].value += r[i].is_negative_factor() ? -push : push;
            }
        }

        calculator_sort<Mode>(r, r + r_size, rng);

        const real old_phi = sign * pi[k];
        const auto q = [&](int j) { return sign * r[j].value + old_phi; };

        int attractive = 0;
        while (attractive != r_size && q(attractive) < Zero)
            ++attractive;

        const int bmin =
          std::clamp(b[k].min + b[k].negative_factor, 0, r_size);
        const int bmax =
          std::clamp(b[k].max + b[k].negative_factor, bmin, r_size);
        const int count = std::clamp(attractive, bmin, bmax);

        real lo = count > 0 ? q(count - 1) : -inf;
        real hi = count < r_size ? q(count) : +inf;
        if (count > bmin)
            hi = std::min(hi, Zero);
        if (count < bmax)
            lo = std::max(lo, Zero);

        real phi = Zero;
        if (std::isfinite(lo) && std::isfinite(hi))
            phi = (lo + hi) / Two;
        else if (std::isfinite(lo))
            phi = lo;
        else if (std::isfinite(hi))
            phi = hi;

        const real gap =
          count > 0 && count < r_size ? q(count) - q(count - 1) : Zero;

        return { sign * phi,
                 delta + (kappa / (One - kappa)) * gap,
                 k,
                 count };
    }

    template<typename Xtype>
    void apply(const row_update& u, Xtype& x) noexcept
    {
        pi[u.k] = u.pi;

        auto [row_first, row_last] = ap.row(u.k);
        const int r_size = static_cast<int>(row_last - row_first);
        const rc_data* r = &R[row_first->value];
        const real d = sign * u.d;

        for (int i = 0; i != r_size; ++i) {
            const auto var = row_first + r[i].id;
            const bool selected = i < u.count;

            if (selected != r[i].is_negative_factor())
                x.set(var->column);
            else
                x.unset(var->column);

            P[var->value] += selected ? d : -d;
        }

        bx_expects(is_valid_constraint(*this, u.k, x));
    }

    bool shares_variable_with_batch(int k) const noexcept
    {
        for (auto [it, et] = ap.row(k); it != et; ++it)
            if (variable_batch[static_cast<std::size_t>(it->column)] ==
                batch_id)
                return true;

        return false;
    }

    void add_to_batch(int k) noexcept
    {
        for (auto [it, et] = ap.row(k); it != et; ++it)
            variable_batch[static_cast<std::size_t>(it->column)] = batch_id;
    }

    template<typename Xtype>
    bool flush_batch(Xtype& x) noexcept
    {
        bool pi_sign_changed = false;

        for (const auto& u : batch) {
            pi_sign_changed =
              pi_sign_changed || is_signbit_change(pi[u.k], u.pi);
            apply(u, x);
        }

        batch.clear();
        ++batch_id;

        return pi_sign_changed;
    }

    // Rows conflicting with the current batch are deferred to the next one
    // (greedy colouring) so that batches stay large enough to be computed in
    // parallel.
    template<typename Xtype, typename Iterator>
    bool update_rows(Xtype& x,
                     Iterator first,
                     Iterator last,
                     real kappa,
                     real delta,
                     real theta,
                     real obj_amp) noexcept
    {
        bool pi_sign_changed = false;

        pending.clear();
        for (; first != last; ++first)
            pending.emplace_back(constraint(first));

        while (!pending.empty()) {
            deferred.clear();

            for (const int k : pending) {
                if (shares_variable_with_batch(k)) {
                    deferred.emplace_back(k);
                } else {
                    add_to_batch(k);
                    batch.emplace_back(
                      compute_row(x, k, kappa, delta, theta, obj_amp));
                }
            }

            pi_sign_changed = flush_batch(x) || pi_sign_changed;
            pending.swap(deferred);
        }

        return pi_sign_changed;
    }

    template<typename Xtype, typename Iterator>
    bool push_and_compute_update_row(Xtype& x,
                                     Iterator first,
                                     Iterator last,
                                     real kappa,
                                     real delta,
                                     real theta,
                                     real obj_amp) noexcept
    {
        logger::log("push-update-row {} {} {}\n", kappa, delta, theta);

        return update_rows(x, first, last, kappa, delta, theta, obj_amp);
    }

    template<typename Xtype, typename Iterator>
    bool compute_update_row(Xtype& x,
                            Iterator first,
                            Iterator last,
                            real kappa,
                            real delta,
                            real theta) noexcept
    {
        logger::log("update-row {} {} {}\n", kappa, delta, theta);

        return update_rows(x, first, last, kappa, delta, theta, Zero);
    }
};

template<typename Mode, typename Cost>
static result
solve_or_optimize(const context& ctx, const problem& pb, bool is_optimization)
{
    if (ctx.parameters.debug) {
        using Solver = solver_inequalities_101coeff_buffered<Mode, Cost, true>;

        return is_optimization ? optimize_problem<Solver, Mode, Cost>(ctx, pb)
                               : solve_problem<Solver, Mode, Cost>(ctx, pb);
    } else {
        using Solver =
          solver_inequalities_101coeff_buffered<Mode, Cost, false>;

        return is_optimization ? optimize_problem<Solver, Mode, Cost>(ctx, pb)
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
             : solve_or_optimize<Mode, baryonyx::itm::quadratic_cost_type>(
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
solve_inequalities_101_buffered(const context& ctx, const problem& pb)
{
    info(ctx, "  - solve_inequalities_101_buffered\n");
    return select_mode(ctx, pb, false);
}

result
optimize_inequalities_101_buffered(const context& ctx, const problem& pb)
{
    info(ctx, "  - optimize_inequalities_101_buffered\n");
    return select_mode(ctx, pb, true);
}

} // namespace itm
} // namespace baryonyx
