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

namespace lp {
namespace details {

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
                  return std::get<0>(lhs) < std::get<0>(rhs);
              });
}

template<typename iteratorT>
void sort(iteratorT begin, iteratorT end, maximize_tag)
{
    std::sort(begin, end,
              [](const auto& lhs, const auto& rhs)
              {
                  return std::get<0>(rhs) < std::get<0>(lhs);
              });
}

template<typename modeT>
class default_algorithm
{
    using tag = modeT;
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

        details::sort(r.begin(), r.end(), tag());

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
        : m(std::distance(pb_.equal_constraints.begin(),
                          pb_.equal_constraints.end()))
        , n(std::distance(pb_.vars.values.begin(),
                          pb_.vars.values.end()))
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

} // namespace details

inline
result simple_wedelin(double kappa, double delta, double theta,
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
