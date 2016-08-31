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
#include "wedelin.hpp"
#include "utils.hpp"
#include <Eigen/Core>
#include <iterator>

namespace lp {
std::tuple<double, double, double, long>
get_parameters(const std::vector<parameter>& params)
{
    double kappa {0.001}, delta{0.001}, theta{0.0001};
    long limit {1000};

    for (std::size_t i{0}, e{params.size()}; i != e; ++i) {
        if (params[i].name == "kappa")
            kappa = (params[i].type == parameter::tag::real) ?
                params[i].d : kappa;
        else if (params[i].name == "theta")
            theta = (params[i].type == parameter::tag::real) ?
                params[i].d : theta;
        else if (params[i].name == "delta")
            delta = (params[i].type == parameter::tag::real) ?
                params[i].d : delta;
        else if (params[i].name == "limit")
            limit = (params[i].type == parameter::tag::integer) ?
                params[i].l : limit;
    }

    std::printf("Solve: kappa(%f) theta(%f) delta(%f) - limit(%ld)\n",
                kappa, theta, delta, limit);

    return std::make_tuple(kappa, delta, theta, limit);
}

template<typename variableT>
bool is_boolean_variable(const variableT& vars)
{
    for (const auto& elem : vars)
        if (elem.type != variable_type::binary)
            return false;

    return true;
}

template<typename constraintsT>
bool is_boolean_coefficient(const constraintsT csts)
{
    for (const auto& cst : csts)
        for (const auto& elem : cst.elements)
            if (elem.factor < 0 or elem.factor > 1)
                return false;

    return true;
}

template<typename variableT>
bool is_integer_variable(const variableT& vars)
{
    for (const auto& elem : vars)
        if (elem.type != variable_type::general)
            return false;

    return true;
}

template<typename constraintsT>
bool is_101_coefficient(const constraintsT csts)
{
    for (const auto& cst : csts)
        for (const auto& elem : cst.elements)
            if (elem.factor < 1 or elem.factor > 1)
                return false;

    return true;
}

result mitm(const problem& pb, const std::vector<parameter>& params)
{
    double kappa, delta, theta;
    long limit;

    std::tie(kappa, delta, theta, limit) = get_parameters(params);

    if (pb.greater_constraints.empty()
        and pb.greater_equal_constraints.empty()
        and pb.less_constraints.empty()
        and pb.less_equal_constraints.empty()) {

        if (is_boolean_coefficient(pb.equal_constraints)
            and is_boolean_variable(pb.vars.values)) {
            wedelin::default_algorithm solver(kappa, delta, theta, limit, pb);
            return solver.results();
        }
    }

    if (is_101_coefficient(pb.equal_constraints)
        and is_101_coefficient(pb.greater_constraints)
        and is_101_coefficient(pb.greater_equal_constraints)
        and is_101_coefficient(pb.less_constraints)
        and is_101_coefficient(pb.less_equal_constraints)
        and is_integer_variable(pb.vars.values)) {
        wedelin::general_problem solver(kappa, delta, theta, limit, pb);
        return solver.results();
    }

    throw lp::solver_error(solver_error::tag::no_solver_available);
}

}
