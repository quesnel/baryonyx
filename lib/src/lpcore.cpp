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

#include "lpformat-consistency.hpp"
#include "lpformat-io.hpp"
#include "mitm.hpp"
#include <fstream>
#include <lpcore>

namespace lp {

problem
make_problem(const std::string& filename)
{
    std::ifstream ifs;
    ifs.exceptions(std::ifstream::badbit);
    ifs.open(filename);

    return details::read_problem(ifs);
}

problem
make_problem(std::istream& is)
{
    is.exceptions(std::ifstream::badbit);

    return details::read_problem(is);
}

std::ostream&
operator<<(std::ostream& os, const problem& p)
{
    details::problem_writer pw(p, os);

    return os;
}

result
solve(problem& pb)
{
    check(pb);

    std::map<std::string, parameter> params;

    return mitm(pb, params);
}

result
solve(problem& pb, const std::map<std::string, parameter>& params)
{
    check(pb);

    return mitm(pb, params);
}

template <typename functionT, typename variablesT>
int
compute_function(const functionT& fct, const variablesT& vars) noexcept
{
    int v{ 0 };

    for (auto& f : fct)
        v += f.factor * vars[f.variable_index];

    return v;
}

bool
is_valid_solution(const problem& pb,
                  const std::deque<int>& variable_value) noexcept
{
    for (auto& cst : pb.equal_constraints) {
        if (compute_function(cst.elements, variable_value) != cst.value) {
            printf("constraint %s (=) fails\n", cst.label.c_str());
            return false;
        }
    }

    for (auto& cst : pb.greater_constraints) {
        if (compute_function(cst.elements, variable_value) <= cst.value) {
            printf("constraint %s (>) fails\n", cst.label.c_str());
            return false;
        }
    }

    for (auto& cst : pb.greater_equal_constraints) {
        if (compute_function(cst.elements, variable_value) < cst.value) {
            printf("constraint %s (>=) fails\n", cst.label.c_str());
            return false;
        }
    }

    for (auto& cst : pb.less_constraints) {
        if (compute_function(cst.elements, variable_value) >= cst.value) {
            printf("constraint %s (<) fails\n", cst.label.c_str());
            return false;
        }
    }

    for (auto& cst : pb.greater_constraints) {
        if (compute_function(cst.elements, variable_value) > cst.value) {
            printf("constraint %s (<=) fails\n", cst.label.c_str());
            return false;
        }
    }

    return true;
}
}
