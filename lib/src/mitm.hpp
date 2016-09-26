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

#ifndef ORG_VLEPROJECT_LP_MITM_HPP
#define ORG_VLEPROJECT_LP_MITM_HPP

#include <lpcore>
#include <Eigen/Core>
#include <type_traits>
#include "utils.hpp"

namespace lp {

template<typename T>
Eigen::Matrix<T, 2, Eigen::Dynamic>
make_a(index m, index n, const problem& pb)
{
    using matrix = Eigen::Matrix<T, 2, Eigen::Dynamic>;

    matrix A(matrix::Zero(m, n));

    for (lp::index i = 0; i != m; ++i) {
        for (const auto& elem : pb.equal_constraints[i].elements) {
            auto j = elem.variable_index;
            A(i, j) = true;
        }
    }

    return A;
}

template<typename T>
Eigen::Matrix<T, 1, Eigen::Dynamic>
make_b(index m, const problem& p)
{
    using vector = Eigen::Matrix<T, 1, Eigen::Dynamic>;

    vector b(vector::Zero(m));

    for (lp::index i = 0; i != m; ++i)
        b(i) = p.equal_constraints[i].value;

    return b;
}


template<typename T>
Eigen::Matrix<T, 1, Eigen::Dynamic>
make_c(index n, const problem& p)
{
    using vector = Eigen::Matrix<T, 1, Eigen::Dynamic>;

    vector c(vector::Zero(n));

    for (const auto& elem : p.objective.elements)
        c(elem.variable_index) = elem.factor;

    return c;
}

template<typename T>
Eigen::Matrix<T, 2, Eigen::Dynamic>
make_inequality_a(index m, index n, const problem& pb)
{
    using matrix = Eigen::Matrix<T, 2, Eigen::Dynamic>;

    matrix A(matrix::Zero(m, n));

    lp::index i {0};
    for (const auto& cst : pb.equal_constraints) {
        for (const auto& elem : cst.elements) {
            auto j = elem.variable_index;
            A(i, j) = elem.factor;
        }
        ++i;
    }

    for (const auto& cst : pb.greater_equal_constraints) {
        for (const auto& elem : cst.elements) {
            auto j = elem.variable_index;
            A(i, j) = elem.factor;
        }
        ++i;
    }

    for (const auto& cst : pb.less_equal_constraints) {
        for (const auto& elem : cst.elements) {
            auto j = elem.variable_index;
            A(i, j) = elem.factor;
        }
        ++i;
    }

    Ensures(i == n, "make_inequality i != n");

    return A;
}

template<typename T>
Eigen::Matrix<T, 2, Eigen::Dynamic>
make_inequality_b(index m, index n, const problem& p)
{
    using matrix = Eigen::Matrix<T, 2, Eigen::Dynamic>;

    matrix b(matrix::Zero(2, m));

    lp::index i = 0;
    for (const auto& elem : p.equal_constraints) {
        b(0, i) = elem.value;
        b(1, i) = elem.value;
        ++i;
    }

    for (const auto& elem : p.greater_equal_constraints) {
        b(0, i) = elem.value;
        b(1, i) = std::numeric_limits<int>::max();
        ++i;
    }

    for (const auto& elem : p.greater_equal_constraints) {
        b(0, i) = std::numeric_limits<int>::min();
        b(1, i) = elem.value;
        ++i;
    }

    Ensures(i == n, "make_inequality_b != n");

    return b;
}

/**
 * Remove lower bound < 0 without changing the problem's structure other
 * than modifying the right-hand side and adding a constant to the
 * objective function.
 *
 * \param pb Initial problem definition.
 *
 * \return New problem without lower bound variable less than 0.
 */
inline
problem adapt_problem(const problem& pb)
{
    problem ret {pb};
    int constant {0};

    for (std::size_t i {0}, e{ret.vars.values.size()}; i != e; ++i) {
        if (ret.vars.values[i].min < 0) {
            constant += ret.vars.values[i].min;
            ret.vars.values[i].max -= ret.vars.values[i].min;
            ret.vars.values[i].min = 0;
        }
    }

    ret.objective.constant = constant;

    return ret;
}

result mitm(const problem& pb, const std::vector<parameter>& params);

}

#endif
