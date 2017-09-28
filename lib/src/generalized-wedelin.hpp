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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_GENERALIZED_WEDELIN_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_GENERALIZED_WEDELIN_HPP

#include "utils.hpp"

#include <iterator>

namespace baryonyx {

// namespace generalized {

// using A_type = Eigen::MatrixXi;
// using b_type = Eigen::Matrix<double, 2, Eigen::Dynamic>;
// using c_type = Eigen::VectorXf;
// using x_type = Eigen::VectorXi;
// using P_type = Eigen::MatrixXf;
// using pi_type = Eigen::VectorXf;
// using u_type = Eigen::VectorXi;

// struct maximize_tag
// {
// };
// struct minimize_tag
// {
// };

// template <typename iteratorT>
// void
// simple_constraint_calculator_sort(iteratorT begin, iteratorT end,
// minimize_tag)
// {
//     std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
//         return lhs.value < rhs.value;
//     });
// }

// template <typename iteratorT>
// void
// simple_constraint_calculator_sort(iteratorT begin, iteratorT end,
// maximize_tag)
// {
//     std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
//         return rhs.value < lhs.value;
//     });
// }

// template <typename iteratorT>
// void
// advanded_constraint_calculator_sort(iteratorT begin,
//                                     iteratorT end,
//                                     minimize_tag)
// {
//     std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
//         return lhs.value < rhs.value;
//     });
// }

// template <typename iteratorT>
// void
// advanded_constraint_calculator_sort(iteratorT begin,
//                                     iteratorT end,
//                                     maximize_tag)
// {
//     std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
//         return rhs.value < lhs.value;
//     });
// }

// template <typename modeT>
// struct constraint_calculator
// {
//     using mode_type = modeT;

//     virtual void update_row(index k,
//                             double kappa,
//                             double delta,
//                             double theta) = 0;

//     virtual ~constraint_calculator()
//     {
//     }
// };

// template <typename modeT>
// struct simple_constraint_calculator : public constraint_calculator<modeT>
// {
//     using mode_type = typename constraint_calculator<modeT>::mode_type;

//     struct r_data
//     {
//         r_data(double value_, index index_)
//           : value(value_)
//           , id(index_)
//         {
//         }

//         double value;
//         index id;
//     };

//     A_type& A;
//     b_type& b;
//     c_type& c;
//     x_type& x;
//     P_type& P;
//     pi_type& pi;
//     std::vector<index> I;
//     std::vector<r_data> r;
//     index m;
//     index n;

//     simple_constraint_calculator(const problem& pb,
//                                  index m_,
//                                  index n_,
//                                  A_type& A_,
//                                  b_type& b_,
//                                  c_type& c_,
//                                  x_type& x_,
//                                  P_type& P_,
//                                  pi_type& pi_)
//       : A(A_)
//       , b(b_)
//       , c(c_)
//       , x(x_)
//       , P(P_)
//       , pi(pi_)
//       , m(m_)
//       , n(n_)
//     {
//         for (const auto& elem : pb.objective.elements)
//             x(elem.variable_index) = c(elem.variable_index) <= 0;

//         for (index i = 0; i != m; ++i) {
//             for (const auto& elem : pb.equal_constraints[i].elements) {
//                 I.emplace_back(elem.variable_index);
//                 r.emplace_back(0.0, elem.variable_index);
//             }
//         }
//     }

//     void update_row(index k, double kappa, double delta, double theta)
//     override
//     {
//         for (auto i : I)
//             P(k, i) *= theta;

//         for (index i = 0, endi = numeric_cast<index>(I.size()); i != endi;
//              ++i) {
//             double sum_a_pi{ 0 };
//             double sum_a_p{ 0 };

//             for (index h = 0; h != m; ++h) {
//                 if (A(h, I[i])) {
//                     sum_a_pi += A(h, I[i]) * pi(h);
//                     sum_a_p += A(h, I[i]) * P(h, I[i]);
//                 }
//             }

//             r[i].value = c(I[i]) - sum_a_pi - sum_a_p;
//             r[i].id = I[i];
//         }

//         simple_constraint_calculator_sort(r.begin(), r.end(), mode_type());

//         pi(k) += ((r[b(k)].value + r[b(k) - 1].value) / 2.0);

//         double d = delta + ((kappa / (1.0 - kappa)) *
//                             (r[b(k) - 1].value - r[b(k)].value));

//         for (index j = 0; j < b(k); ++j) {
//             x(r[j].id) = 1;
//             P(k, r[j].id) -= d;
//         }

//         for (index j = b(k), endj = numeric_cast<index>(I.size()); j !=
//         endj;
//              ++j) {
//             x(r[j].id) = 0;
//             P(k, r[j].id) -= -d;
//         }
//     }
// };

// template <typename modeT>
// struct advanded_constraint_calculator : public constraint_calculator<modeT>
// {
//     using mode_type = typename constraint_calculator<modeT>::mode_type;

//     struct r_data
//     {
//         r_data(double value_, index index_)
//           : value(value_)
//           , id(index_)
//         {
//         }
//         double value;
//         index id;
//     };

//     A_type& A;
//     b_type& b;
//     c_type& c;
//     x_type& x;
//     P_type& P;
//     pi_type& pi;
//     u_type& u;
//     std::vector<index> I;
//     std::vector<r_data> r;
//     std::vector<index> C;
//     index m;
//     index n;

//     advanded_constraint_calculator(const problem& pb,
//                                    index k,
//                                    index m_,
//                                    index n_,
//                                    A_type& A_,
//                                    b_type& b_,
//                                    c_type& c_,
//                                    x_type& x_,
//                                    P_type& P_,
//                                    pi_type& pi_,
//                                    u_type& u_)
//       : A(A_)
//       , b(b_)
//       , c(c_)
//       , x(x_)
//       , P(P_)
//       , pi(pi_)
//       , u(u_)
//       , m(m_)
//       , n(n_)
//     {
//         for (const auto& elem : pb.objective.elements)
//             x(elem.variable_index) = c(elem.variable_index) <= 0;

//         for (index i = 0; i != m; ++i) {
//             for (const auto& elem : pb.equal_constraints[i].elements) {
//                 I.emplace_back(elem.variable_index);

//                 r.emplace_back(0.0, elem.variable_index);
//             }
//         }

//         for (index i{ 0 }; i != n; ++i)
//             if (A(k, i) < 0)
//                 C.emplace_back(i);
//     }

//     void update_row(index k,
//                     double /*kappa*/,
//                     double /*delta*/,
//                     double theta) override
//     {
//         for (auto i : I)
//             P(k, i) *= theta;

//         for (index i = 0, endi = numeric_cast<index>(I.size()); i != endi;
//              ++i) {
//             double sum_a_pi{ 0 };
//             double sum_a_p{ 0 };

//             for (index h = 0; h != m; ++h) {
//                 if (A(h, I[i])) {
//                     sum_a_pi += A(h, I[i]) * pi(h);
//                     sum_a_p += A(h, I[i]) * P(h, I[i]);
//                 }
//             }

//             r[i].value = c(I[i]) - sum_a_pi - sum_a_p;
//             r[i].id = I[i];
//         }

//         for (auto i : C) {
//             r[i].value *= -1;
//             A(k, i) *= -1;
//             P(k, i) *= -1;
//         }

//         auto sum_aki_ui{ 0 };
//         for (auto i : C)
//             sum_aki_ui += A(k, i) * u(i);

//         b(k, 0) += sum_aki_ui;
//         b(k, 1) += sum_aki_ui;

//         advanded_constraint_calculator_sort(r.begin(), r.end(),
//         mode_type());

//         // pi(k) += ((r[b(k)].value + r[b(k) - 1].value) / 2.0);

//         // double d = delta + ((kappa / (1.0 - kappa)) *
//         //                     (r[b(k) - 1].value - r[b(k)].value));

//         // for (index j = 0; j < b(k); ++j) {
//         //     x(r[j].id) = 1;
//         //     P(k, r[j].id) -= d;
//         // }

//         // for (index j = b(k), endj = numeric_cast<index>(I.size()); j !=
//         // endj;
//         //      ++j) {
//         //     x(r[j].id) = 0;
//         //     P(k, r[j].id) -= -d;
//         // }

//         // Undo adjustment of row lower and upper bound.
//         b(k, 0) -= sum_aki_ui;
//         b(k, 1) -= sum_aki_ui;

//         // Clean-up: correct negated costs and adjust value of negated
//         // variables.
//         for (auto i : C) {
//             A(k, i) *= -1;
//             P(k, i) *= -1;
//             x(i) = u(i) - x(i);
//         }
//     }
// };

// // struct generalized_constraint_calculator : public constraint_calculator
// // {
// // public:
// // };

// template <typename modeT>
// class solver
// {
//     std::vector<constraint_calculator<modeT>*> row_updaters;
//     std::vector<index> R;
//     index m;
//     index n;
//     A_type A;
//     b_type b;
//     c_type c;
//     x_type x;
//     P_type P;
//     pi_type pi;
//     u_type u;

// public:
//     using mode_type = modeT;

//     solver(const lp::problem& pb)
//       : m(pb.equal_constraints.size() + pb.less_equal_constraints.size() +
//           pb.greater_equal_constraints.size())
//       , n(pb.vars.values.size())
//       , A(A_type::Zero(m, n))
//       , b(b_type::Zero(2, m))
//       , c(c_type::Zero(n))
//       , x(x_type::Zero(n))
//       , P(P_type::Zero(m, n))
//       , pi(pi_type::Zero(m))
//       , u(u_type::Zero(n))
//     {
//         for (std::size_t i{ 0 }, e = pb.vars.values.size(); i != e; ++i)
//             u(i) = pb.vars.values[i].max;

//         index i = 0;
//         for (const auto& elem : pb.equal_constraints) {
//             for (const auto& cst : elem.elements)
//                 A(i, cst.variable_index) = cst.factor;

//             b(0, i) = elem.value;
//             b(1, i) = elem.value;
//             ++i;
//         }

//         for (const auto& elem : pb.less_equal_constraints) {
//             for (const auto& cst : elem.elements)
//                 A(i, cst.variable_index) = cst.factor;

//             b(0, i) = -std::numeric_limits<double>::infinity();
//             b(1, i) = elem.value;
//             ++i;
//         }

//         for (const auto& elem : pb.greater_equal_constraints) {
//             for (const auto& cst : elem.elements)
//                 A(i, cst.variable_index) = cst.factor;

//             b(0, i) = elem.value;
//             b(1, i) = std::numeric_limits<double>::infinity();
//             ++i;
//         }

//         for (const auto& elem : pb.objective.elements)
//             c(elem.variable_index) = elem.factor;

//         index k = 0;
//         for (; k != m; ++k) {
//             bool simple = b(0, k) == b(1, k);

//             if (simple) {
//                 for (index i{ 0 }; i != n; ++i) {
//                     if (A(k, i) != 0 and A(k, i) != 1) {
//                         simple = false;
//                         break;
//                     }

//                     if (A(k, i) == 1) {
//                         if (pb.vars.values[i].type != variable_type::binary)
//                         {
//                             simple = false;
//                             break;
//                         }
//                     }
//                 }
//             }

//             if (simple) {
//                 row_updaters.emplace_back(
//                   new simple_constraint_calculator<mode_type>(
//                     pb, m, n, A, b, c, x, P, pi));
//             } else {
//                 bool advanced = true;

//                 for (index i{ 0 }; i != n; ++i) {
//                     if (A(k, i) != 0 and A(k, i) != 1 and A(k, i) != -1) {
//                         advanced = false;
//                         break;
//                     }
//                 }

//                 if (advanced) {
//                     row_updaters.emplace_back(
//                       new advanded_constraint_calculator<mode_type>(
//                         pb, k, m, n, A, b, c, x, P, pi, u));
//                 } else
//                     throw solver_failure(
//                       solver_error_tag::no_solver_available);
//             }
//         }
//     }

//     virtual ~solver()
//     {
//         std::for_each(row_updaters.begin(),
//                       row_updaters.end(),
//                       [](auto* elem) { delete elem; });
//     }

//     bool compute(double kappa, double delta, double theta)
//     {
//         R.clear();

//         for (index k{ 0 }; k != m; ++k) {
//             int v = 0;

//             for (index i{ 0 }; i != n; ++i)
//                 v += A(k, i) * x(i);

//             if (v != b(k))
//                 R.push_back(k);
//         }

//         if (R.empty())
//             return true;

//         for (index k : R)
//             row_updaters[k]->update_row(k, kappa, delta, theta);

//         return false;
//     }

//     result results() const
//     {
//         result ret;
//         ret.status = result_status::success;

//         return ret;
//     }
// };

// /**
//  * Remove lower bound < 0 without changing the problem's structure other
//  than
//  * modifying the right-hand side and adding a constant to the objective
//  * function.
//  *
//  * \param pb Initial problem definition.
//  *
//  * \return New problem without lower bound variable less than 0.
//  */
// inline problem
// adapt_problem(const problem& pb)
// {
//     problem ret(pb);
//     int constant{ 0 };

//     for (std::size_t i{ 0 }, e{ ret.vars.values.size() }; i != e; ++i) {
//         if (ret.vars.values[i].min != 0) {
//             constant += ret.vars.values[i].min;
//             ret.vars.values[i].max -= ret.vars.values[i].min;
//             ret.vars.values[i].min = 0;
//         }
//     }

//     assert(pb.less_constraints.empty() and pb.greater_constraints.empty());

//     ret.objective.constant = constant;

//     return ret;
// }

// } // generalized

inline result
generalized_wedelin(double /* kappa */,
                    double /* delta */,
                    double /* theta */,
                    long /* limit */,
                    const problem& /* pb */)
{
    //     using maximize_solver =
    //     generalized::solver<generalized::maximize_tag>;
    //     using minimize_solver =
    //     generalized::solver<generalized::minimize_tag>;
    //     long int i{ 0 };

    //     if (pb.type == lp::objective_function_type::maximize) {
    //         maximize_solver slv(generalized::adapt_problem(pb));
    //         while (not slv.compute(kappa, delta, theta and i++ != limit))
    //             ;

    //         return slv.results();
    //     }

    //     minimize_solver slv(generalized::adapt_problem(pb));

    //     while (not slv.compute(kappa, delta, theta) and i++ != limit)
    //         ;

    //     return slv.results();

    return baryonyx::result{};
}
}

#endif
