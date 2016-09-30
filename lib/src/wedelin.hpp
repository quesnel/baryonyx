/* Copyright (C) 2016 INRA
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

#ifndef ORG_VLEPROJECT_LP_WEDELIN_HPP
#define ORG_VLEPROJECT_LP_WEDELIN_HPP

#include "mitm.hpp"
#include "utils.hpp"
#include <Eigen/Core>
#include <iterator>
#include <iostream>

namespace lp {
namespace details {

Eigen::MatrixXi
make_a(index m, index n, const problem& pb)
{
    Eigen::MatrixXi A = Eigen::MatrixXi::Zero(m, n);

    for (lp::index i = 0; i != m; ++i) {
        for (const auto& elem : pb.equal_constraints[i].elements) {
            auto j = elem.variable_index;
            A(i, j) = 1;
        }
    }

    return A;
}

Eigen::VectorXi
make_b(index m, const problem& p)
{
    Eigen::VectorXi b = Eigen::VectorXi::Zero(m);

    for (lp::index i = 0; i != m; ++i)
        b(i) = p.equal_constraints[i].value;

    return b;
}

Eigen::VectorXf
make_c(index n, const problem& p)
{
    Eigen::VectorXf c = Eigen::VectorXf::Zero(n);

    for (const auto& elem : p.objective.elements)
        c(elem.variable_index) = elem.factor;

    return c;
}


struct maximize_tag {};
struct minimize_tag {};

template<typename T>
struct solver_tag {
    using type = T;
};

template<typename iteratorT>
void sort(iteratorT begin, iteratorT end, minimize_tag)
{
    std::sort(begin, end,
              [](const auto& lhs, const auto& rhs)
              {
                  return lhs.value < rhs.value;
              });
}

template<typename iteratorT>
void sort(iteratorT begin, iteratorT end, maximize_tag)
{
    std::sort(begin, end,
              [](const auto& lhs, const auto& rhs)
              {
                  return rhs.value < lhs.value;
              });
}

template<typename modeT>
class default_algorithm
{
    struct r_data {
        r_data(double value_, index index_)
            : value(value_)
            , id(index_)
        {}

        double value;
        index id;
    };

    using tag = modeT;
    index m, n;
    const problem& pb;
    Eigen::MatrixXi A;
    Eigen::VectorXi b;
    Eigen::VectorXf c;
    Eigen::VectorXi x;
    Eigen::MatrixXf P;
    Eigen::VectorXf pi;
    std::vector<std::vector<index>> I;
    std::vector<std::vector<r_data>> r;
    double kappa, delta, theta;
    index loop;
    bool optimal;

    void update_row(int k)
    {
        for (auto i : I[k])
            P(k, i) *= theta;

        for (index i = 0, endi = numeric_cast<index>(I[k].size());
             i != endi; ++i) {
            double sum_a_pi {0};
            double sum_a_p {0};

            assert(k < numeric_cast<index>(I.size()));
            assert(i < numeric_cast<index>(I[k].size()));

            for (index h = 0; h != m; ++h) {
                if (A(h, I[k][i])) {
                    sum_a_pi += A(h, I[k][i]) * pi(h);
                    sum_a_p += A(h, I[k][i]) * P(h, I[k][i]);
                }
            }

            r[k][i].value = c(I[k][i]) - sum_a_pi - sum_a_p;
            r[k][i].id = I[k][i];
        }

        details::sort(r[k].begin(), r[k].end(), tag());

        pi(k) += ((r[k][b(k)].value + r[k][b(k) - 1].value) / 2.0);

        double d = delta
            + ((kappa / (1.0 - kappa))
               * (r[k][b(k) - 1].value - r[k][b(k)].value));

        for (index j = 0; j < b(k); ++j) {
            x(r[k][j].id) = 1;
            P(k, r[k][j].id) -= d;
        }

        for (index j = b(k), endj = numeric_cast<index>(I[k].size());
             j != endj; ++j) {
            x(r[k][j].id) = 0;
            P(k, r[k][j].id) -= -d;
        }
    }

public:
    default_algorithm(double kappa_, double delta_, double theta_,
                      long limit, const problem& pb_)
        : m(std::distance(pb_.equal_constraints.begin(),
                          pb_.equal_constraints.end()))
        , n(std::distance(pb_.vars.values.begin(),
                          pb_.vars.values.end()))
        , pb(pb_)
        , A(make_a(m, n, pb))
        , b(make_b(m, pb))
        , c(make_c(n, pb))
        , x(Eigen::VectorXi::Zero(n))
        , P(Eigen::MatrixXf::Zero(m, n))
        , pi(Eigen::VectorXf::Zero(m))
        , I(m)
        , r(m)
        , kappa(kappa_)
        , delta(delta_)
        , theta(theta_)
        , loop(0)
        , optimal(false)
    {
        Ensures(kappa >= 0 and kappa < 1, "kappa [0, 1[");
        Ensures(delta >= 0, "delta [0, +oo[");
        Ensures(theta >= 0 and theta <= 1, "theta [0, 1]");
        Ensures(m > 0, "equal_constraints number must be > 0");
        Ensures(n > 0, "variable number must be > 0");

        for (const auto& elem : pb.objective.elements)
            x(elem.variable_index) = c(elem.variable_index) <= 0;

        for (index i = 0; i != m; ++i) {
            for (const auto& elem : pb.equal_constraints[i].elements) {
                I[i].emplace_back(elem.variable_index);
                r[i].emplace_back(0.0, elem.variable_index);
            }
        }

        // std::cout << "A:\n" << A << '\n'
        //           << "b:\n" << b.transpose() << '\n'
        //           << "c:\n" << c.transpose() << '\n'
        //           << "x:\n" << x.transpose() << '\n'
        //           << '\n';

        std::vector<index> R;
        while (loop != limit) {
            for (index k {0}; k != m; ++k) {
                int v = 0;
                for (index i {0}; i != n; ++i)
                    v += A(k, i) * x(i);

                if (v != b(k))
                    R.push_back(k);
            }

            if (R.empty()) {
                optimal = true;
                return;
            }

            // std::cout << "P:\n" << P << '\n';
            // std::cout << "pi:\n" << pi.transpose() << '\n';

            for (auto k : R)
                update_row(k);

            R.clear();
            ++loop;
        }
    }

    double compute_value() const
    {
        double ret {0};

        for (auto& elem : pb.objective.elements)
            ret += elem.factor * x(elem.variable_index);

        return ret;
    }

    result results() const
    {
        std::cout << "result: " << n << '\n';

        result ret;
        ret.loop = loop;
        ret.optimal = optimal;
        ret.value  = compute_value();
        ret.variable_name.resize(n);
        ret.variable_value.resize(n, 0);

        for (index i = 0; i != n; ++i) {
            ret.variable_name[i] = pb.vars.names[i];
            if (x(i))
                ret.variable_value[i] = true;
        }

        return ret;
    }
};

} // namespace details

inline
result
simple_wedelin(double kappa, double delta, double theta,
               long limit, const problem& pb)
{
    using maximize_solver = details::default_algorithm<details::maximize_tag>;
    using minimize_solver = details::default_algorithm<details::minimize_tag>;

    if (pb.type == lp::objective_function_type::maximize) {
        maximize_solver solver(kappa, delta, theta, limit, pb);
        return solver.results();
    }

    minimize_solver solver(kappa, delta, theta, limit, pb);
    return solver.results();
}

} // namespace lp

#endif
