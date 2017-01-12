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

#include "mitm.hpp"
#include "generalized-wedelin.hpp"
#include "inequalities-1coeff.hpp"
#include "utils.hpp"
#include "wedelin.hpp"
#include <Eigen/Core>
#include <iterator>

namespace lp {

std::tuple<double, double, double, long>
get_parameters(const std::map<std::string, parameter>& params)
{
    double kappa{ 0.001 }, delta{ 0.001 }, theta{ 0.0001 };
    long limit{ 1000 };

    {
        auto it = params.find("kappa");
        if (it->second.type == parameter::tag::real)
            kappa = it->second.d;
    }

    {
        auto it = params.find("theta");
        if (it->second.type == parameter::tag::real)
            theta = it->second.d;
    }

    {
        auto it = params.find("delta");
        if (it->second.type == parameter::tag::real)
            delta = it->second.d;
    }

    {
        auto it = params.find("limit");
        if (it->second.type == parameter::tag::integer)
            limit = it->second.l;
    }

    std::printf("Solve: kappa(%f) theta(%f) delta(%f) - limit(%ld)\n",
                kappa,
                theta,
                delta,
                limit);

    return std::make_tuple(kappa, delta, theta, limit);
}

template <typename variableT>
bool
is_boolean_variable(const variableT& vars)
{
    for (const auto& elem : vars)
        if (elem.type != variable_type::binary)
            return false;

    return true;
}

template <typename constraintsT>
bool
is_boolean_coefficient(const constraintsT csts)
{
    for (const auto& cst : csts)
        for (const auto& elem : cst.elements)
            if (elem.factor < 0 or elem.factor > 1)
                return false;

    return true;
}

template <typename variableT>
bool
is_integer_variable(const variableT& vars)
{
    for (const auto& elem : vars)
        if (elem.type != variable_type::general)
            return false;

    return true;
}

template <typename constraintsT>
bool
is_101_coefficient(const constraintsT csts)
{
    for (const auto& cst : csts)
        for (const auto& elem : cst.elements)
            if (elem.factor < -1 or elem.factor > 1)
                return false;

    return true;
}

result
mitm(const problem& pb, const std::map<std::string, parameter>& params)
{
    double kappa, delta, theta;
    long limit;

    std::tie(kappa, delta, theta, limit) = get_parameters(params);

    if (pb.greater_constraints.empty() and
        pb.greater_equal_constraints.empty() and
        pb.less_constraints.empty() and pb.less_equal_constraints.empty() and
        is_boolean_coefficient(pb.equal_constraints) and
        is_boolean_variable(pb.vars.values)) {
        printf("simple_wedelin\n");
        return simple_wedelin(kappa, delta, theta, limit, pb);
    }

    // printf("test: %d\n", pb.equal_constraints.empty());
    // printf("test: %d\n", pb.greater_equal_constraints.empty());
    // printf("test: %d\n", pb.less_equal_constraints.empty());

    if ((not pb.equal_constraints.empty() or
         not pb.greater_equal_constraints.empty() or
         not pb.less_equal_constraints.empty())) {
        printf("test: %d\n", pb.greater_constraints.empty());
        printf("test: %d\n", pb.less_constraints.empty());
        printf("test: %d\n", is_101_coefficient(pb.equal_constraints));
        printf("test: %d\n", is_101_coefficient(pb.greater_equal_constraints));
        printf("test: %d\n", is_101_coefficient(pb.less_equal_constraints));
        printf("test: %d\n", is_boolean_variable(pb.vars.values));

        for (int i{ 0 },
             e(numeric_cast<int>(pb.less_equal_constraints.size()));
             i != e;
             ++i)
            for (const auto& elem : pb.less_equal_constraints[i].elements)
                if (elem.factor < -1 or elem.factor > 1)
                    printf("Error %d (factor=%d\n", i, elem.factor);

        if (pb.greater_constraints.empty() and pb.less_constraints.empty() and
            is_101_coefficient(pb.equal_constraints) and
            is_101_coefficient(pb.greater_equal_constraints) and
            is_101_coefficient(pb.less_equal_constraints) and
            is_boolean_variable(pb.vars.values)) {
            printf("inequalities-1coeff\n");
            return lp::inequalities_1coeff_wedelin(
              kappa, delta, theta, limit, pb);
        }
    }

    if ((is_101_coefficient(pb.equal_constraints) or
         is_101_coefficient(pb.greater_constraints) or
         is_101_coefficient(pb.greater_equal_constraints) or
         is_101_coefficient(pb.less_constraints) or
         is_101_coefficient(pb.less_equal_constraints)) and
        is_integer_variable(pb.vars.values)) {
        printf("generalized_wedelin\n");
        return lp::generalized_wedelin(kappa, delta, theta, limit, pb);
    }

    throw lp::solver_error(solver_error::tag::no_solver_available);
}
} // namespace lp
