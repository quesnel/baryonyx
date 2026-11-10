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
#include "parallel-for.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace baryonyx {
namespace itm {

// Rows are updated by batches of rows sharing no variable, given by a greedy
// coloring of the row conflict graph computed once: inside a batch the row
// computations read nothing another row of the batch writes (pi, P and x of
// its own variables), so they run concurrently without locks. Each row
// draws its ties from its own random engine, seeded from the batch, so the
// result does not depend on the number of threads.
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
    std::vector<int> batch_rows;
    std::vector<int> row_color;
    std::vector<int> color_start;
    std::vector<int> rows_by_color;
    std::vector<long> row_work;
    std::unique_ptr<parallel_for> pool;

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
      , c(c_)
      , m(m_)
      , n(n_)
    {
        batch.reserve(static_cast<std::size_t>(m));
        batch_rows.reserve(static_cast<std::size_t>(m));
        rows_by_color.resize(static_cast<std::size_t>(m));
        color_rows();

        row_work.resize(static_cast<std::size_t>(m));
        for (int k = 0; k != m; ++k) {
            long work = 0;
            for (auto [it, et] = ap.row(k); it != et; ++it) {
                auto [cb, ce] = ap.column(it->column);
                work += static_cast<long>(ce - cb);
            }
            row_work[static_cast<std::size_t>(k)] = work;
        }

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

    void set_thread_number(int thread_number)
    {
        pool = thread_number > 1
                 ? std::make_unique<parallel_for>(thread_number)
                 : nullptr;
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
                           real obj_amp,
                           random_engine& row_rng) noexcept
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

        calculator_sort<Mode>(r, r + r_size, row_rng);

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

    void color_rows()
    {
        row_color.assign(static_cast<std::size_t>(m), -1);
        std::vector<int> forbidden_by;

        for (int k = 0; k != m; ++k) {
            for (auto [it, et] = ap.row(k); it != et; ++it)
                for (auto [cb, ce] = ap.column(it->column); cb != ce; ++cb)
                    if (const int color =
                          row_color[static_cast<std::size_t>(cb->row)];
                        color >= 0)
                        forbidden_by[static_cast<std::size_t>(color)] = k;

            int color = 0;
            while (color != length(forbidden_by) &&
                   forbidden_by[static_cast<std::size_t>(color)] == k)
                ++color;

            if (color == length(forbidden_by))
                forbidden_by.emplace_back(-1);

            row_color[static_cast<std::size_t>(k)] = color;
        }

        color_start.assign(forbidden_by.size() + 1u, 0);
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

        return pi_sign_changed;
    }

    template<typename Xtype>
    void compute_batch(const Xtype& x,
                       real kappa,
                       real delta,
                       real theta,
                       real obj_amp)
    {
        batch.resize(batch_rows.size());
        const auto seed = rng();

        auto compute = [&](int i) {
            const int k = batch_rows[static_cast<std::size_t>(i)];
            random_engine row_rng(
              seed + static_cast<random_engine::result_type>(k));

            batch[static_cast<std::size_t>(i)] =
              compute_row(x, k, kappa, delta, theta, obj_amp, row_rng);
        };

        const int size = static_cast<int>(batch_rows.size());

        // Below about a thousand multiply-adds the batch is cheaper than
        // handing it to the workers.
        constexpr long parallel_min_work = 1000;

        long work = 0;
        if (pool)
            for (const int k : batch_rows)
                work += row_work[static_cast<std::size_t>(k)];

        if (pool && work >= parallel_min_work)
            pool->run(size, compute);
        else
            for (int i = 0; i != size; ++i)
                compute(i);
    }

    // The requested rows are bucketed by color, keeping their order inside a
    // color, and each color is one batch.
    template<typename Xtype, typename Iterator>
    bool update_rows(Xtype& x,
                     Iterator first,
                     Iterator last,
                     real kappa,
                     real delta,
                     real theta,
                     real obj_amp)
    {
        bool pi_sign_changed = false;

        std::fill(color_start.begin(), color_start.end(), 0);
        for (auto it = first; it != last; ++it)
            ++color_start[static_cast<std::size_t>(
              row_color[static_cast<std::size_t>(constraint(it))]) + 1u];

        for (std::size_t c = 1; c < color_start.size(); ++c)
            color_start[c] += color_start[c - 1];

        const auto colors = color_start.size() - 1u;
        batch_rows.assign(color_start.begin(), color_start.end() - 1);
        for (auto it = first; it != last; ++it) {
            const int k = constraint(it);
            const auto color =
              static_cast<std::size_t>(row_color[static_cast<std::size_t>(k)]);
            rows_by_color[static_cast<std::size_t>(batch_rows[color]++)] = k;
        }

        for (std::size_t color = 0; color != colors; ++color) {
            const auto begin = rows_by_color.begin() + color_start[color];
            const auto end = rows_by_color.begin() + color_start[color + 1];
            if (begin == end)
                continue;

            batch_rows.assign(begin, end);
            compute_batch(x, kappa, delta, theta, obj_amp);
            pi_sign_changed = flush_batch(x) || pi_sign_changed;
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
    return select_cost<minimize_tag>(ctx, pb, is_optimization);
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
