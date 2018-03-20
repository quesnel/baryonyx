/* Copyright (C) 2016-2018 INRA
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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_COMMON_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_COMMON_HPP

#include <baryonyx/core>

#include <tuple>

#include <cassert>

namespace baryonyx {
namespace itm {

struct maximize_tag
{};

struct minimize_tag
{};

/**
 * x_type is a std::vector<bool> instead of a baryonyx::fixed_array<bool> to
 * use the specialized version of vector, which is used for elements of type
 * bool and optimizes for space.
 */
using x_type = std::vector<bool>;

inline double
best_solution_value(const baryonyx::result& res) noexcept
{
    assert(res.status == baryonyx::result_status::success);
    assert(not res.solutions.empty());

    return res.solutions.back().value;
}

inline bool
is_better_solution(const baryonyx::result& lhs,
                   const baryonyx::result& rhs,
                   maximize_tag) noexcept
{
    return best_solution_value(lhs) > best_solution_value(rhs);
}

inline bool
is_better_solution(const baryonyx::result& lhs,
                   const baryonyx::result& rhs,
                   minimize_tag) noexcept
{
    return best_solution_value(lhs) < best_solution_value(rhs);
}

inline std::ostream&
operator<<(std::ostream& os, const baryonyx::affected_variables& var)
{
    std::size_t i = 0, e = var.names.size();

    assert(var.names.size() == var.values.size());

    for (; i != e; ++i)
        os << var.names[i] << ": " << (var.values[i] ? 1 : 0) << '\n';

    return os;
}

struct best_solution_writer
{
    const baryonyx::result& res;

    best_solution_writer(const baryonyx::result& res_)
      : res(res_)
    {}
};

inline std::ostream&
operator<<(std::ostream& os, const best_solution_writer& writer)
{
    assert(writer.res.status == baryonyx::result_status::success);
    assert(not writer.res.solutions.empty());
    assert(writer.res.variable_name.size() ==
           writer.res.solutions.back().variables.size());

    std::size_t i = 0, e = writer.res.variable_name.size();

    for (; i != e; ++i)
        os << writer.res.variable_name[i] << ": "
           << (writer.res.solutions.back().variables[i] ? 1 : 0) << '\n';

    return os;
}
}
}

#endif
