/* Copyright (C) 2017 INRA
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

#ifndef ORG_VLEPROJECT_LP_INEQUALITIES_1COEFF_HPP
#define ORG_VLEPROJECT_LP_INEQUALITIES_1COEFF_HPP

#include "lpcore-out"
#include "mitm.hpp"
#include "utils.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iterator>
#include <random>

namespace lp {
namespace inequalities_1coeff {

using A_type = Eigen::MatrixXi;
using b_type = Eigen::Matrix<double, 2, Eigen::Dynamic>;
using c_type = Eigen::VectorXf;
using x_type = Eigen::VectorXi;
using P_type = Eigen::MatrixXf;
using pi_type = Eigen::VectorXf;
using u_type = Eigen::VectorXi;

enum class constraint_order
{
    none,
    reversing,
    random_sorting,
    infeasibility_decr,
    infeasibility_incr
};

struct maximize_tag
{
};
struct minimize_tag
{
};

template <typename iteratorT>
void
calculator_sort(iteratorT begin, iteratorT end, minimize_tag)
{
    std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
        return lhs.value < rhs.value;
    });
}

template <typename iteratorT>
void
calculator_sort(iteratorT begin, iteratorT end, maximize_tag)
{
    std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
        return rhs.value < lhs.value;
    });
}

bool
stop_iterating(double value, minimize_tag)
{
    return value > 0;
}

bool
stop_iterating(double value, maximize_tag)
{
    return value < 0;
}

template <typename modeT>
struct constraint_calculator
{
    using mode_type = modeT;

    struct r_data
    {
        r_data(double value_, index index_)
          : value(value_)
          , id(index_)
        {
        }

        double value;
        index id;
    };

    A_type& A;
    b_type& b;
    c_type& c;
    x_type& x;
    P_type& P;
    pi_type& pi;
    std::vector<index> I; // Stores variables with non null coefficient.
    std::vector<r_data> r;
    std::vector<index> C; // Stores variables with negative coefficient.
    index m;
    index n;

    constraint_calculator(index k,
                          index m_,
                          index n_,
                          A_type& A_,
                          b_type& b_,
                          c_type& c_,
                          x_type& x_,
                          P_type& P_,
                          pi_type& pi_)
      : A(A_)
      , b(b_)
      , c(c_)
      , x(x_)
      , P(P_)
      , pi(pi_)
      , m(m_)
      , n(n_)
    {
        for (index i = 0; i != n; ++i) {
            if (A(k, i)) {
                I.emplace_back(i);
                r.emplace_back(0.0, i);
            }

            if (A(k, i) < 0)
                C.emplace_back(i);
        }
    }

    void update_row(index k, double kappa, double delta, double theta)
    {
        for (auto i : I)
            P(k, i) *= theta;

        for (index i = 0, endi = numeric_cast<index>(I.size()); i != endi;
             ++i) {
            double sum_a_pi{ 0 };
            double sum_a_p{ 0 };

            for (index h = 0; h != m; ++h) {
                if (A(h, I[i])) {
                    sum_a_pi += A(h, I[i]) * pi(h);
                    sum_a_p += A(h, I[i]) * P(h, I[i]);
                }
            }

            r[i].value = c(I[i]) - sum_a_pi - sum_a_p;
            r[i].id = I[i];
        }

        /*
         * Negate reduced costs and coefficients of these variables.
         */

        for (std::size_t i{ 0 }, e{ C.size() }; i != e; ++i) {
            int variable = C[i];

            auto it = std::find_if(r.begin(), r.end(), [variable](auto elem) {
                return elem.id == variable;
            });

            assert(it != r.end());

            it->value *= -1;
            A(k, variable) *= -1;
            P(k, variable) *= -1;
        }

        calculator_sort(r.begin(), r.end(), mode_type());

        b(0, k) += C.size();
        b(1, k) += C.size();

        /*
         * The bkmin and bkmax constraint bounds are not equal and can be
         * assigned to -infinity or +infinity. We have to scan the r vector
         * and search a value j such as b(0, k) <= Sum A(k, r[j]) < b(1,
         * k).
         */
        index i{ 0 };
        const index endi{ numeric_cast<index>(r.size()) };
        index sum{ 0 };

        if (b(0, k) == 0 and b(1, k) == 0) {
            for (index j{ 0 }; j != endi; ++j)
                x(r[j].id) = 0;

            goto undo;
        }

        {
            index low = -1, up = -1;

            for (; i != endi; ++i) {
                sum += A(k, r[i].id);

                if (b(0, k) <= sum) {
                    low = i;
                    up = i;
                    break;
                }
            }

            if (low == -1) {
                throw solver_error(solver_error::tag::unrealisable_constraint);
            }

            ++i;
            ++up;
            for (; i != endi; ++i) {
                sum += A(k, r[i].id);

                if (stop_iterating(r[i].value, mode_type()) or sum > b(1, k)) {
                    up -= 1;
                    break;
                }
            }

            i = up;
        }

        if (endi == 1) {
            pi(k) += ((r[0].value + r[0].value) / 2.0);

            double d =
              delta + ((kappa / (1.0 - kappa)) * (r[0].value - r[0].value));

            x(r[0].id) = 1;
            P(k, r[0].id) += d;

            for (index j{ 1 }; j != endi; ++j) {
                x(r[0].id) = 0;
                P(k, r[0].id) -= d;
            }

        } else {
            if (i >= endi)    // If the upper bound was not reached, we use
                i = endi - 1; // the last element in r.

            index first = i, second = i + 1;
            if (first < 0) {
                first = 0;
                second = 1;
            }

            if (second > endi) {
                second -= 2;
            }

            pi(k) += ((r[first].value + r[second].value) / 2.0);

            double d = delta + ((kappa / (1.0 - kappa)) *
                                (r[second].value - r[first].value));

            index j{ 0 };
            for (; j <= i; ++j) {
                x(r[j].id) = 1;
                P(k, r[j].id) += d;
            }

            for (; j != endi; ++j) {
                x(r[j].id) = 0;
                P(k, r[j].id) -= d;
            }
        }

    undo:
        b(0, k) -= C.size();
        b(1, k) -= C.size();

        /*
         * Clean up: correct negated costs and adjust value of negated
         * variables.
         */
        for (std::size_t i{ 0 }, e{ C.size() }; i != e; ++i) {
            A(k, C[i]) *= -1;
            P(k, C[i]) *= -1;

            x[C[i]] = 1 - x[C[i]];
        }
    }
};

template <typename modeT>
class solver
{
    std::deque<constraint_calculator<modeT>> row_updaters;
    std::vector<index> R;
    index m;
    index n;
    A_type A;
    b_type b;
    c_type c;
    x_type x;
    P_type P;
    pi_type pi;
    u_type u;
    const lp::problem& pb;
    bool solution_found;
    std::default_random_engine rng;

public:
    using mode_type = modeT;

    solver(const lp::problem& pb)
      : m(pb.equal_constraints.size() + pb.less_equal_constraints.size() +
          pb.greater_equal_constraints.size())
      , n(pb.vars.values.size())
      , A(A_type::Zero(m, n))
      , b(b_type::Zero(2, m))
      , c(c_type::Zero(n))
      , x(x_type::Zero(n))
      , P(P_type::Zero(m, n))
      , pi(pi_type::Zero(m))
      , u(u_type::Zero(n))
      , pb(pb)
      , solution_found(false)
      , rng(std::chrono::system_clock::now().time_since_epoch().count())
    {
        for (std::size_t i{ 0 }, e = pb.vars.values.size(); i != e; ++i)
            u(i) = pb.vars.values[i].max;

        index i = 0;
        for (const auto& elem : pb.equal_constraints) {
            for (const auto& cst : elem.elements)
                A(i, cst.variable_index) = cst.factor;

            b(0, i) = elem.min;
            b(1, i) = elem.max;
            ++i;
        }

        for (const auto& elem : pb.less_equal_constraints) {
            for (const auto& cst : elem.elements)
                A(i, cst.variable_index) = cst.factor;

            if (elem.min != elem.max) {
                b(0, i) = elem.min;
                b(1, i) = elem.max;
            } else {
                b(0, i) = -std::numeric_limits<double>::infinity();
                b(1, i) = elem.max;
            }

            ++i;
        }

        for (const auto& elem : pb.greater_equal_constraints) {
            for (const auto& cst : elem.elements)
                A(i, cst.variable_index) = cst.factor;

            if (elem.min != elem.max) {
                b(0, i) = elem.min;
                b(1, i) = elem.max;
            } else {
                b(0, i) = elem.min;
                b(1, i) = std::numeric_limits<double>::infinity();
            }

            ++i;
        }

        for (const auto& elem : pb.objective.elements) {
            c(elem.variable_index) += elem.factor;
            x(elem.variable_index) = c(elem.variable_index) <= 0 ? 1 : 0;
        }

        for (index k = 0; k != m; ++k)
            row_updaters.emplace_back(k, m, n, A, b, c, x, P, pi);
    }

    std::size_t compute_none(double kappa, double delta, double theta)
    {
        R.clear();

        for (index k{ 0 }; k != m; ++k) {
            int v = 0;

            for (index i{ 0 }; i != n; ++i)
                v += A(k, i) * x(i);

            if (not(b(0, k) <= v and v <= b(1, k)))
                R.push_back(k);
        }

        for (auto it{ R.cbegin() }, et{ R.cend() }; it != et; ++it)
            row_updaters[*it].update_row(*it, kappa, delta, theta);

        if (R.empty())
            solution_found = true;

        return R.size();
    }

    std::size_t compute_reversing(double kappa, double delta, double theta)
    {
        R.clear();

        for (index k{ 0 }; k != m; ++k) {
            int v = 0;

            for (index i{ 0 }; i != n; ++i)
                v += A(k, i) * x(i);

            if (not(b(0, k) <= v and v <= b(1, k)))
                R.push_back(k);
        }

        for (auto it{ R.crbegin() }, et{ R.crend() }; it != et; ++it)
            row_updaters[*it].update_row(*it, kappa, delta, theta);

        if (R.empty())
            solution_found = true;

        return R.size();
    }

    std::size_t compute_random_sorting(double kappa,
                                       double delta,
                                       double theta)
    {
        R.clear();

        for (index k{ 0 }; k != m; ++k) {
            int v = 0;

            for (index i{ 0 }; i != n; ++i)
                v += A(k, i) * x(i);

            if (not(b(0, k) <= v and v <= b(1, k)))
                R.push_back(k);
        }

        std::shuffle(R.begin(), R.end(), rng);

        for (index k : R)
            row_updaters[k].update_row(k, kappa, delta, theta);

        if (R.empty())
            solution_found = true;

        return R.size();
    }

    std::size_t compute_infeasibility_decr(double kappa,
                                           double delta,
                                           double theta)
    {
        std::vector<std::pair<index, index>> currentR;

        for (index k{ 0 }; k != m; ++k) {
            int v = 0;

            for (index i{ 0 }; i != n; ++i)
                v += A(k, i) * x(i);

            if (not(b(0, k) <= v and v <= b(1, k))) {
                currentR.emplace_back(k, 0);

                if (v < b(0, k))
                    currentR.back().second = b(0, k) - v;
                else
                    currentR.back().second = v - b(1, k);
            }
        }

        std::sort(currentR.begin(),
                  currentR.end(),
                  [](const auto& lhs, const auto& rhs) {

                      // [this](const auto& lhs, const auto& rhs) {
                      //     if (rhs.second == lhs.second) {
                      //         std::bernoulli_distribution d(.5);
                      //         return d(rng);
                      //     }

                      return rhs.second < lhs.second;
                  });

        std::bernoulli_distribution d(.5);
        std::size_t i{ 0 };
        const std::size_t e{ currentR.size() };

        for (; i != e; ++i) {
            if (i + 1 == e)
                break;

            if (currentR[i].second == currentR[i + 1].second)
                if (d(rng))
                    std::swap(currentR[i].first, currentR[i + 1].second);
        }

        for (i = 0; i != e; ++i)
            row_updaters[currentR[i].first].update_row(
              currentR[i].first, kappa, delta, theta);

        if (currentR.empty())
            solution_found = true;

        return currentR.size();
    }

    std::size_t compute_infeasibility_incr(double kappa,
                                           double delta,
                                           double theta)
    {
        std::vector<std::pair<index, index>> currentR;

        for (index k{ 0 }; k != m; ++k) {
            int v = 0;

            for (index i{ 0 }; i != n; ++i)
                v += A(k, i) * x(i);

            if (not(b(0, k) <= v and v <= b(1, k))) {
                currentR.emplace_back(k, 0);

                if (v < b(0, k))
                    currentR.back().second = b(0, k) - v;
                else
                    currentR.back().second = v - b(1, k);
            }
        }

        std::sort(currentR.begin(),
                  currentR.end(),
                  [](const auto& lhs, const auto& rhs) {

                      // [this](const auto& lhs, const auto& rhs) {
                      //     if (rhs.second == lhs.second) {
                      //         std::bernoulli_distribution d(.5);
                      //         return d(rng);
                      //     }

                      return lhs.second < rhs.second;
                  });

        std::bernoulli_distribution d(.5);
        std::size_t i{ 0 };
        const std::size_t e{ currentR.size() };

        for (; i != e; ++i) {
            if (i + 1 == e)
                break;

            if (currentR[i].second == currentR[i + 1].second)
                if (d(rng))
                    std::swap(currentR[i].first, currentR[i + 1].second);
        }

        for (i = 0; i != e; ++i)
            row_updaters[currentR[i].first].update_row(
              currentR[i].first, kappa, delta, theta);

        if (currentR.empty())
            solution_found = true;

        return currentR.size();
    }

    double compute_value() const
    {
        double ret = pb.objective.constant;

        for (auto& elem : pb.objective.elements)
            ret += elem.factor * x(elem.variable_index);

        return ret;
    }

    result results() const
    {
        result ret;
        ret.method = "inequalities_101coeff_wedelin";
        ret.variables = n;
        ret.constraints = m;
        ret.value = compute_value();
        ret.solution_found = solution_found;
        ret.variable_name.resize(n);
        ret.variable_value.resize(n, 0);

        for (index i = 0; i != n; ++i) {
            ret.variable_name[i] = pb.vars.names[i];
            ret.variable_value[i] = x(i);
        }

        return ret;
    }
};

template <typename modeT>
inline result
run(const problem& pb,
    constraint_order order,
    double theta,
    double delta,
    long int limit,
    double kappa_min,
    double kappa_step,
    double kappa_max,
    double alpha,
    long int w)
{
    using mode_type = modeT;

    auto begin = std::chrono::steady_clock::now();
    long int i2{ 0 };
    double kappa_old{ 0 };
    double kappa = kappa_min;

    solver<mode_type> slv(pb);
    result best;
    best.remaining_constraints = std::numeric_limits<index>::max();

    for (long int i{ 0 }; i < limit; ++i) {
        index remaining;

        switch (order) {
            case inequalities_1coeff::constraint_order::none:
                remaining = slv.compute_none(kappa, delta, theta);
                break;
            case inequalities_1coeff::constraint_order::reversing:
                remaining = slv.compute_reversing(kappa, delta, theta);
                break;
            case inequalities_1coeff::constraint_order::random_sorting:
                remaining = slv.compute_random_sorting(kappa, delta, theta);
                break;
            case inequalities_1coeff::constraint_order::infeasibility_decr:
                remaining =
                  slv.compute_infeasibility_decr(kappa, delta, theta);
                break;
            case inequalities_1coeff::constraint_order::infeasibility_incr:
                remaining =
                  slv.compute_infeasibility_incr(kappa, delta, theta);
                break;
        }

        auto current = slv.results();
        current.loop = i;
        current.remaining_constraints = remaining;
        current.begin = begin;
        current.end = std::chrono::steady_clock::now();

        if (remaining < best.remaining_constraints) {
            auto t = std::chrono::duration_cast<std::chrono::duration<double>>(
              current.end - current.begin);

            printf(
              "Constraints remaining: %ld at %fs\n", remaining, t.count());
            best = current;
        }

        if (current.solution_found) {
            std::cout << current << '\n';
            return current;
        }

        if (i2 <= w) {
            kappa = kappa_min;
            i2++;
        } else {
            i2 = 0;
            kappa =
              kappa_old +
              kappa_step * std::pow(static_cast<double>(remaining) /
                                      static_cast<double>(current.constraints),
                                    alpha);
            kappa_old = kappa;
        }

        if (kappa > kappa_max) {
            std::cout << "\nFail: kappa-max reached\n";
            return best;
        }
    }

    std::cout << "\nFail: limit reached\n";
    return best;
}

} // inequalities_1coeff

double
get_real(const std::map<std::string, parameter>& params,
         std::string param,
         double def)
{
    auto it = params.find(param);
    if (it == params.cend())
        return def;

    if (it->second.type != parameter::tag::real)
        return def;

    return it->second.d;
}

long int
get_integer(const std::map<std::string, parameter>& params,
            std::string param,
            double def)
{
    auto it = params.find(param);
    if (it == params.cend())
        return def;

    if (it->second.type != parameter::tag::integer)
        return def;

    return it->second.l;
}

inequalities_1coeff::constraint_order
get_constraint_order(const std::map<std::string, parameter>& params,
                     std::string param,
                     inequalities_1coeff::constraint_order def)
{
    auto it = params.find(param);
    if (it == params.cend())
        return def;

    if (it->second.type != parameter::tag::string)
        return def;

    if (it->second.s == "none")
        return inequalities_1coeff::constraint_order::none;
    if (it->second.s == "reversing")
        return inequalities_1coeff::constraint_order::reversing;
    if (it->second.s == "random-sorting")
        return inequalities_1coeff::constraint_order::random_sorting;
    if (it->second.s == "infeasibility-decr")
        return inequalities_1coeff::constraint_order::infeasibility_decr;
    if (it->second.s == "infeasibility-incr")
        return inequalities_1coeff::constraint_order::infeasibility_incr;

    return def;
}

const char*
constraint_order_to_string(inequalities_1coeff::constraint_order type)
{
    static const char* ret[] = { "none",
                                 "reversing",
                                 "random-sorting",
                                 "infeasibility-decr",
                                 "infeasibility-incr" };

    switch (type) {
        case inequalities_1coeff::constraint_order::none:
            return ret[0];
        case inequalities_1coeff::constraint_order::reversing:
            return ret[1];
        case inequalities_1coeff::constraint_order::random_sorting:
            return ret[2];
        case inequalities_1coeff::constraint_order::infeasibility_decr:
            return ret[3];
        case inequalities_1coeff::constraint_order::infeasibility_incr:
            return ret[4];
    }

    return nullptr;
}

inline result
inequalities_1coeff_wedelin(const problem& pb,
                            const std::map<std::string, parameter>& params)
{
    namespace ine_1 = lp::inequalities_1coeff;

    auto order = get_constraint_order(
      params, "constraint-order", ine_1::constraint_order::infeasibility_decr);

    double theta = get_real(params, "theta", 0.5);
    double delta = get_real(params, "delta", 0.5);
    long int limit = get_integer(params, "limit", 100l);
    double kappa_min = get_real(params, "kappa-min", 0.0);
    double kappa_step = get_real(params, "kappa-step", 0.0001);
    double kappa_max = get_real(params, "kappa-max", 0.6);
    double alpha = get_real(params, "alpha", 2.0);
    long int w = get_integer(params, "w", 20l);

    printf("inequalities_1coeff_wedelin constraint-order: %s theta: %f delta "
           "%f limit=%ld "
           "kappa-min: %f kappa-step: %f kappa-max: %f alpha: %f w: %ld\n",
           constraint_order_to_string(order),
           theta,
           delta,
           limit,
           kappa_min,
           kappa_step,
           kappa_max,
           alpha,
           w);

    if (pb.type == lp::objective_function_type::maximize)
        return ine_1::run<ine_1::maximize_tag>(pb,
                                               order,
                                               theta,
                                               delta,
                                               limit,
                                               kappa_min,
                                               kappa_step,
                                               kappa_max,
                                               alpha,
                                               w);

    return ine_1::run<ine_1::minimize_tag>(pb,
                                           order,
                                           theta,
                                           delta,
                                           limit,
                                           kappa_min,
                                           kappa_step,
                                           kappa_max,
                                           alpha,
                                           w);
}

} // lp

#endif
