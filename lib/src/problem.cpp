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

#include "problem.hpp"

namespace baryonyx {

void
clear(raw_problem& pb)
{
    std::vector<objective_function_element>().swap(pb.objective.elements);
    std::vector<objective_quadratic_element>().swap(pb.objective.qelements);

    std::vector<constraint>().swap(pb.equal_constraints);
    std::vector<constraint>().swap(pb.greater_constraints);
    std::vector<constraint>().swap(pb.less_constraints);

    std::vector<std::string_view>().swap(pb.vars.names);
    std::vector<variable_value>().swap(pb.vars.values);

    pb.type = objective_function_type::maximize;
}

} // baryonyx namespace
