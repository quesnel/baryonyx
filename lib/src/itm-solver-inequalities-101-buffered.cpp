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

#include <cstring>

namespace baryonyx {
namespace itm {

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

    struct computation_buffer
    {
        computation_buffer(const real pi_,
                           const real d_,
                           const int k_,
                           const int selected_,
                           bool use_loop_)
          : pi(pi_)
          , d(d_)
          , k(k_)
          , selected(selected_)
          , use_loop(use_loop_)
        {}

        real pi;      // pi update
        real d;       // delta difference to apply to P matrix
        int k;        // Constraint id
        int selected; // -1: unselect all, >= R_size: select all

        // true if select_variable use the loop, false if select_variable
        // function use the negative coefficient function.
        bool use_loop;
    };

    sparse_matrix<int> ap;
    std::unique_ptr<real[]> P;
    std::unique_ptr<std::int8_t[]> A;
    std::unique_ptr<rc_data[]> R;
    std::unique_ptr<bound_factor[]> b;
    std::unique_ptr<real[]> pi;
    std::vector<computation_buffer> buffer;

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
        buffer.reserve(static_cast<size_t>(m));

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
        std::memset(P.get(), 0, sizeof(real) * ap.length());
        std::memset(pi.get(), 0, sizeof(real) * m);
        std::memset(R.get(), 0, sizeof(computation_buffer) * m);
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

    template<typename ConstraintIterator>
    void decrease_preference(ConstraintIterator first,
                             ConstraintIterator last,
                             const real theta) noexcept
    {
        for (; first != last; ++first) {
            const int k = constraint(first);
            auto [row_first, row_last] = ap.row(k);

            for (; row_first != row_last; ++row_first) {
                P[row_first->value] *= theta;
            }
        }
    }

    template<typename Xtype, typename ConstraintIterator>
    void compute_reduced_costs(Xtype& x,
                               ConstraintIterator first,
                               ConstraintIterator last) noexcept
    {
        for (; first != last; ++first) {
            const int k = constraint(first);
            auto [row_first, row_last] = ap.row(k);
            int i = 0;

            for (; row_first != row_last; ++row_first) {
                real sum_a_pi = Zero;
                real sum_a_p = Zero;

                auto [cbegin, cend] = ap.column(row_first->column);

                for (; cbegin != cend; ++cbegin) {
                    //auto a = static_cast<real>(A[cbegin->value]);
                    sum_a_pi += /*a **/ pi[cbegin->row];
                    sum_a_p += /*a **/ P[cbegin->value];
                }

                R[row_first->value].id = i++;
                R[row_first->value].value =
                  A[row_first->value] *
                  ((c(row_first->column, x) - sum_a_pi - sum_a_p));
                R[row_first->value].f = A[row_first->value];
            }
        }
    }

    template<typename ConstraintIterator>
    void sort_reduced_costs(ConstraintIterator first,
                            ConstraintIterator last) noexcept
    {
        for (; first != last; ++first) {
            const int k = constraint(first);
            auto [row_first, row_last] = ap.row(k);
            const auto length = row_last - row_first;

            calculator_sort<Mode>(
              &R[row_first->value], &R[row_first->value] + length, rng);
        }
    }

    template<typename Xtype, typename ConstraintIterator>
    void compute_objective_amplifier(Xtype& x,
                                     ConstraintIterator first,
                                     ConstraintIterator last,
                                     const real obj_amp) noexcept
    {
        for (; first != last; ++first) {
            const int k = constraint(first);
            auto [row_first, row_last] = ap.row(k);

            for (; row_first != row_last; ++row_first)
                R[row_first->value].value += obj_amp * c(row_first->column, x);
        }
    }

    template<typename ConstraintIterator>
    void select_variables(ConstraintIterator first,
                          ConstraintIterator last,
                          const real kappa,
                          const real delta) noexcept
    {
        for (; first != last; ++first) {
            const int k = constraint(first);
            auto [row_first, row_last] = ap.row(k);
            auto first_copy = row_first;

            bx_expects(row_last - row_first > 0);
            bx_expects(row_last - row_first < INT_MAX);

            const int row_length = static_cast<int>(row_last - row_first);
            int sum = 0;
            int selected = -1;
            int classic_selected = -1;
            int i = 0;
            bool found_a_solution = b[k].min <= sum && sum <= b[k].max;

            for (; first_copy != row_last; ++first_copy, ++i) {
                sum += R[first_copy->value].f;

                if (b[k].min <= sum && sum <= b[k].max) {
                    selected = i;
                    found_a_solution = true;
                }

                if (found_a_solution &&
                    stop_iterating<Mode>(R[first_copy->value].value, rng))
                    break;
            }

            // Classic algorithm use in Bastert et al. A Generalized Wedelin
            // Heuristic for Integer Programming in INFORMS Journal on
            // Computing.

            if (b[k].min == b[k].max) {
                classic_selected =
                  std::min(b[k].min + b[k].negative_factor, row_length) - 1;
            } else {
                int bkmin = b[k].min + b[k].negative_factor;
                int bkmax =
                  std::min(b[k].max + b[k].negative_factor, row_length);

                classic_selected = bkmax - 1;
                first_copy = row_first + bkmin;
                for (int i = bkmin; i <= bkmax; ++i, ++first_copy) {
                    if (stop_iterating<Mode>(R[first_copy->value].value,
                                             rng)) {
                        classic_selected = i - 1;
                        break;
                    }
                }
            }

            real pi = 0;
            real d = delta;

            if (selected > classic_selected) {
                found_a_solution = true;
            } else {
                found_a_solution = false;
                selected = classic_selected;
            }

            if (selected < 0) {
                d +=
                  (kappa / (One - kappa)) * (R[row_first->value].value * Two);
            } else if (selected + 1 >= row_length) {
                d += (kappa / (One - kappa)) *
                     (R[row_first->value].value * Middle);
            } else {
                pi = (R[row_first->value + selected].value +
                      R[row_first->value + selected + 1].value) /
                     Two;

                d += (kappa / (One - kappa)) *
                     (R[row_first->value + selected + 1].value -
                      R[row_first->value + selected].value);
            }

            buffer.emplace_back(pi, d, k, selected, found_a_solution);
        }
    }

    void print_PI() const
    {
        fmt::print("PI");

        for (int i = 0; i < m; ++i)
            fmt::print(" {} ({})", pi[i], i);

        fmt::print("\n");
    }

    template<typename Iterator>
    void print_R(Iterator first, Iterator last) const
    {
        fmt::print("R");

        for (; first != last; ++first)
            fmt::print(" [{}, {}, {}]",
                       R[first->value].value,
                       R[first->value].id,
                       R[first->value].f);

        fmt::print("\n");
    }

    template<typename Xtype>
    void affect(Xtype& x) noexcept
    {
        // print_PI();

        for (auto& elem : buffer) {
            // fmt::print("\nAffect: pi:{} d:{} k:{} select:{} use_loop: {}\n",
            //           elem.pi,
            //           elem.d,
            //           elem.k,
            //           elem.selected,
            //           elem.use_loop);

            pi[elem.k] += elem.pi;

            auto [row_first, row_last] = ap.row(elem.k);
            // print_R(row_first, row_last);

            auto first = row_first;
            int i = 0;

            for (; first != row_last && i <= elem.selected; ++i, ++first) {
                const auto var = row_first + R[first->value].id;

                if (elem.use_loop) {
                    x.set(var->column);
                    P[var->value] +=
                      R[first->value].is_negative_factor() ? -elem.d : elem.d;
                } else {
                    if (R[first->value].is_negative_factor()) {
                        x.unset(var->column);
                        P[var->value] -= elem.d;
                    } else {
                        x.set(var->column);
                        P[var->value] += elem.d;
                    }
                }
            }

            for (; first != row_last; ++first) {
                const auto var = row_first + R[first->value].id;

                if (elem.use_loop) {
                    x.unset(var->column);
                    P[var->value] +=
                      R[first->value].is_negative_factor() ? -elem.d : +elem.d;
                } else {
                    if (R[first->value].is_negative_factor()) {
                        x.set(var->column);
                        P[var->value] += elem.d;
                    } else {
                        x.unset(var->column);
                        P[var->value] -= elem.d;
                    }
                }
            }

            bx_expects(is_valid_constraint(*this, elem.k, x));
        }
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

        buffer.clear();

        if (theta != One)
            decrease_preference(first, last, theta);

        compute_reduced_costs(x, first, last);

        if (obj_amp != One)
            compute_objective_amplifier(x, first, last, obj_amp);

        sort_reduced_costs(first, last);
        select_variables(first, last, kappa, delta);

        //@TODO sort buffer structure according to the pi values for examples?
        std::sort(buffer.begin(),
                  buffer.end(),
                  [](const auto& left, const auto& right) {
                      return left.pi > right.pi;
                  });

        // std::shuffle(buffer.begin(), buffer.end(), rng);
        affect(x);

        return false;
    }

    template<typename Xtype, typename Iterator>
    bool compute_update_row(Xtype& x,
                            Iterator first,
                            Iterator last,
                            real kappa,
                            real delta,
                            real theta) noexcept
    {
        logger::log("push-update-row {} {} {}\n", kappa, delta, theta);

        buffer.clear();

        if (theta != One)
            decrease_preference(first, last, theta);

        compute_reduced_costs(x, first, last);
        sort_reduced_costs(first, last);
        select_variables(first, last, kappa, delta);

        //@TODO sort buffer structure according to the pi values for examples?
        std::sort(buffer.begin(),
                  buffer.end(),
                  [](const auto& left, const auto& right) {
                      return left.pi > right.pi;
                  });

        // std::shuffle(buffer.begin(), buffer.end(), rng);
        affect(x);

        return false;
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
