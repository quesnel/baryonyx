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

#ifndef ORG_VLEPROJECT_LP_LPFORMAT_CONSISTENCY_HPP
#define ORG_VLEPROJECT_LP_LPFORMAT_CONSISTENCY_HPP

#include <algorithm>
#include <lpcore>

namespace lp {
namespace details {

inline void
are_variables_used(const problem& pb)
{
    std::vector<bool> vars(pb.vars.names.size(), false);

    for (const auto& elem : pb.objective.elements)
        vars[elem.variable_index] = true;

    for (const auto& cst : pb.equal_constraints)
        for (const auto& elem : cst.elements)
            vars[elem.variable_index] = true;
    for (const auto& cst : pb.greater_constraints)
        for (const auto& elem : cst.elements)
            vars[elem.variable_index] = true;
    for (const auto& cst : pb.greater_equal_constraints)
        for (const auto& elem : cst.elements)
            vars[elem.variable_index] = true;
    for (const auto& cst : pb.less_constraints)
        for (const auto& elem : cst.elements)
            vars[elem.variable_index] = true;
    for (const auto& cst : pb.less_equal_constraints)
        for (const auto& elem : cst.elements)
            vars[elem.variable_index] = true;

    auto it = std::find(std::begin(vars), std::end(vars), false);
    if (it != std::end(vars))
        throw problem_definition_error(
          pb.vars.names[std::distance(std::begin(vars), it)],
          problem_definition_error::tag::variable_not_used);
}

inline void
are_bounds_correct(const problem& pb)
{
    for (std::size_t i{ 0 }, e{ pb.vars.values.size() }; i != e; ++i)
        if (pb.vars.values[i].min > pb.vars.values[i].max)
            throw problem_definition_error(
              pb.vars.names[i], problem_definition_error::tag::bad_bound);
}

} // namespace details

void
check(const problem& pb)
{
    details::are_variables_used(pb);
    details::are_bounds_correct(pb);
}

} // namespace lp

#endif
