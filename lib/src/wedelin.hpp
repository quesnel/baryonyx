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

#ifndef ORG_VLEPROJECT_LP_WEDELIN
#define ORG_VLEPROJECT_LP_WEDELIN

#include "mitm.hpp"
#include "utils.hpp"
#include <Eigen/Core>
#include <iterator>

namespace lp {
namespace wedelin {

class default_algorithm
{
    index m, n;
    const problem& pb;
    Eigen::MatrixXi A;
    Eigen::VectorXi b;
    Eigen::VectorXf c;
    Eigen::VectorXi x;
    Eigen::VectorXi u;
    Eigen::MatrixXf P;
    Eigen::VectorXf pi;
    std::vector<std::vector<index>> I;
    double kappa, delta, theta;
    index loop;
    bool optimal;

    void update_row(int k)
    {
        std::vector<std::tuple<double, index>> r(I[k].size());

        for (auto i : I[k])
            P(k, i) *= theta;

        double sum_a_pi {0};
        double sum_a_p {0};

        for (auto i : I[k]) {
            for (index h = 0; h != m; ++h) {
                if (A(h, i)) {
                    sum_a_pi += A(h, i) * pi(h);
                    sum_a_p += A(h, i) * P(h, i);
                }
            }

            std::get<double>(r[i]) = c[i] - sum_a_pi - sum_a_p;
            std::get<index>(r[i]) = i;
        }

        std::sort(r.begin(), r.end(),
                [](const auto& lhs, const auto& rhs)
                {
                    return std::get<0>(lhs) < std::get<0>(rhs);
                });

        pi(k) += (std::get<double>(r[b(k)]) + std::get<double>(r[b(k) - 1])
                / 2.0);

        double d = delta
            + ((kappa / (1.0 - kappa))
                    * (std::get<double>(r[b(k)])
                       - std::get<double>(r[b(k) - 1])));

        for (index j = 0; j < b(k); ++j) {
            x(std::get<index>(r[j])) = 1;
            P(k, std::get<index>(r[j])) += d;
        }

        for (index j = b(k), endj = numeric_cast<index>(I[k].size());
             j != endj; ++j) {
            x(std::get<index>(r[j])) = 0;
            P(k, std::get<index>(r[j])) -= d;
        }
    }

public:
    default_algorithm(double kappa_, double delta_, double theta_,
                      long limit, const problem& pb_)
        : m(std::distance(pb.equal_constraints.begin(),
                          pb.equal_constraints.end()))
        , n(std::distance(pb.vars.values.begin(), pb.vars.values.end()))
        , pb(pb_)
        , A(make_a<int>(m, n, pb))
        , b(make_b<int>(m, pb))
        , c(make_c<float>(n, pb))
        , x(Eigen::VectorXi::Zero(n))
        , P(Eigen::MatrixXf::Zero(m, n))
        , pi(Eigen::VectorXf::Zero(m))
        , I(m)
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

        for (index i = 0; i != n; ++i)
            u(i) = pb.vars.values[i].max;

        for (index i = 0; i != m; ++i)
            for (const auto& elem : pb.equal_constraints[i].elements)
                I[i].emplace_back(elem.variable_index);

        std::vector<index> R;
        while (loop != limit) {
            for (index k {0}; k != m; ++k) {
                int v = 0;
                for (index i {0}; i != n; ++i)
                    v += A(k, i) * x(i);

                if (v != b(k))
                    R.push_back(v);
            }

            if (R.empty()) {
                optimal = true;
                return;
            }

            for (auto k : R)
                update_row(k);

            R.clear();
            ++loop;
        }
    }

    result results() const
    {
        result ret;
        ret.loop = loop;
        ret.optimal = optimal;

        for (index i = 0; i != n; ++i) {
            if (x(i)) {
                ret.variable_name[i] = pb.vars.names[i];
                ret.variable_value[i] = true;
            }
        }

        return ret;
    }
};

class general_problem
{
    index m, n;
    problem pb;                         // a copy to change the problem
                                        // (remove lower bound < 0)
    Eigen::MatrixXi A;
    Eigen::VectorXi b;
    Eigen::VectorXf c;
    Eigen::VectorXi x;
    Eigen::VectorXi u;
    Eigen::MatrixXf P;
    Eigen::VectorXf pi;
    std::vector<std::vector<index>> I; // Variable with non-zero coefficient
    std::vector<std::vector<index>> C; // Variable with negative coefficient
    double kappa, delta, theta;
    index loop;
    bool optimal;

    void make_c_i()
    {
        I.clear();
        C.clear();
        I.resize(m);
        C.resize(m);

        lp::index i {0};
        for (const auto& cst : pb.equal_constraints) {
            for (const auto& elem : cst.elements) {
                auto j = elem.variable_index;
                if (elem.factor < 0)
                    C[i].push_back(j);
                if (elem.factor)
                    I[i].push_back(j);
            }

            ++i;
        }

        for (const auto& cst : pb.greater_equal_constraints) {
            for (const auto& elem : cst.elements) {
                auto j = elem.variable_index;
                if (elem.factor < 0)
                    C[i].push_back(j);
                if (elem.factor)
                    I[i].push_back(j);
            }

            ++i;
        }

        for (const auto& cst : pb.less_equal_constraints) {
            for (const auto& elem : cst.elements) {
                auto j = elem.variable_index;
                if (elem.factor < 0)
                    C[i].push_back(j);
                if (elem.factor)
                    I[i].push_back(j);
            }

            ++i;
        }

        Ensures(i == n, "make_c_i i != n");
    }

    void update_row(int k)
    {
        std::vector<std::tuple<double, index>> r(I[k].size());

        for (auto i : I[k])
            P(k, i) *= theta;

        double sum_a_pi {0};
        double sum_a_p {0};

        for (auto i : I[k]) {
            for (index h = 0; h != m; ++h) {
                if (A(h, i)) {
                    sum_a_pi += A(h, i) * pi(h);
                    sum_a_p += A(h, i) * P(h, i);
                }
            }

            std::get<double>(r[i]) = c[i] - sum_a_pi - sum_a_p;
            std::get<index>(r[i]) = i;
        }

        // negate reduced costs and coefficients of -1 coefficient.
        for (auto i : C[k]) {
            std::get<double>(r[i]) = - std::get<double>(r[i]);
            A(k, i) = - A(k, i);
            P(k, i) = - P(k, i);
        }

        int a_ki_ui {0};
        for (auto i : C[k])
            a_ki_ui += A(k, i) * u(i);

        b(k, 0) += a_ki_ui;
        b(k, 1) += a_ki_ui;

        std::sort(r.begin(), r.end(),
                  [](const auto& lhs, const auto& rhs)
                  {
                      return std::get<0>(lhs) < std::get<0>(rhs);
                  });




        pi(k) += (std::get<double>(r[b(k)]) + std::get<double>(r[b(k) - 1])
                / 2.0);

        double d = delta
            + ((kappa / (1.0 - kappa))
                    * (std::get<double>(r[b(k)])
                       - std::get<double>(r[b(k) - 1])));

        for (index j = 0; j < b(k); ++j) {
            x(std::get<index>(r[j])) = 1;
            P(k, std::get<index>(r[j])) += d;
        }

        for (index j = b(k), endj = numeric_cast<index>(I[k].size());
             j != endj; ++j) {
            x(std::get<index>(r[j])) = 0;
            P(k, std::get<index>(r[j])) -= d;
        }




        // undo adjustment of row lower and upper bound
        b(k, 0) += a_ki_ui;
        b(k, 1) += a_ki_ui;

        // clean up: correct negated costs and adjust value of negated
        // variables.
        for (auto i : C[k]) {
            A(k, i) = - A(k, i);
            P(k, i) = - P(k, i);
            x(i) = u(i) - x(i);
        }
    }

public:
    general_problem(double kappa_, double delta_, double theta_,
            long limit, const problem& pb_)
        : m(numeric_cast<index>(pb.equal_constraints.size())
            + numeric_cast<index>(pb.greater_equal_constraints.size())
            + numeric_cast<index>(pb.less_equal_constraints.size()))
        , n(numeric_cast<index>(pb.vars.values.size()))
        , pb(adapt_problem(pb_))
        , A(make_inequality_a<int>(m, n, pb))
        , b(make_inequality_b<int>(m, n, pb))
        , c(make_c<float>(n, pb))
        , x(Eigen::VectorXi::Zero(n))
        , u(Eigen::VectorXi::Zero(n))
        , P(Eigen::MatrixXf::Zero(m, n))
        , pi(Eigen::VectorXf::Zero(m))
        , I(m)
        , kappa(kappa_)
        , delta(delta_)
        , theta(theta_)
        , loop(0)
        , optimal(false)
    {
        Ensures(kappa >= 0 and kappa < 1, "kappa [0, 1[");
        Ensures(delta >= 0, "delta [0, +oo[");
        Ensures(theta >= 0 and theta <= 1, "theta [0, 1]");
        Ensures(m > 0, "constraints number must be > 0");
        Ensures(n > 0, "variable number must be > 0");

        for (const auto& elem : pb.objective.elements)
            x(elem.variable_index) = c(elem.variable_index) <= 0;

        // build the C_k, I_k vectors store -1/1 and -1 coefficient.
        make_c_i();

        std::vector<index> R;
        while (loop != limit) {
            for (index k {0}; k != m; ++k) {
                int v = 0;
                for (index i {0}; i != n; ++i)
                    v += A(k, i) * x(i);

                if (v != b(k))
                    R.push_back(v);
            }

            if (R.empty()) {
                optimal = true;
                return;
            }

            for (auto k : R)
                update_row(k);

            R.clear();
            ++loop;
        }
    }

    result results() const
    {
        result ret;
        ret.loop = loop;
        ret.optimal = optimal;

        for (index i = 0; i != n; ++i) {
            if (x(i)) {
                ret.variable_name[i] = pb.vars.names[i];
                ret.variable_value[i] = true;
            }
        }

        return ret;
    }
};

} // namespace details
} // namespace lp

#endif
