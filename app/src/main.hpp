/* Copyright (C) 2016-2019 INRA
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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_MAIN_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_MAIN_HPP

#include <baryonyx/core>

#include <string>

#include <cerrno>
#include <climits>
#include <cmath>

/** @brief Perform a benchmark according to the benchmark description in the
 * json file @e filepath.
 *
 * @param ctx Context with all parameters use to perform optimization.
 * @param filepath Description file.
 * @param name The name of the solver (e.g: cplex-10.0.3, baryonyx-0.2)
 * @return @c true if the processing to the benchmark success, @c false
 *     otherwise.
 *
 * @note Implementation details are in the benchmark.cpp file.
 */
bool
benchmark(const baryonyx::context_ptr& ctx,
          std::string filepath,
          std::string name);

inline double
to_double(std::string s, double bad_value) noexcept
{
    char* c;
    errno = 0;
    double value = std::strtod(s.c_str(), &c);

    if ((errno == ERANGE && (value == HUGE_VAL || value == -HUGE_VAL)) ||
        (value == 0.0 && c == s.c_str()))
        return bad_value;

    return value;
}

inline int
to_int(std::string s, int bad_value) noexcept
{
    char* c;
    errno = 0;
    long value = std::strtol(s.c_str(), &c, 10);

    if ((errno == ERANGE && (value == LONG_MIN || value == LONG_MAX)) ||
        (value == 0 && c == s.c_str()))
        return bad_value;

    if (value < INT_MIN)
        return INT_MIN;

    if (value > INT_MAX)
        return INT_MAX;

    return static_cast<int>(value);
}

constexpr const char*
file_format_error_format(baryonyx::file_format_error_tag failure) noexcept
{
    constexpr const char* const tag[] = {
        "end of file",     "unknown",
        "already defined", "incomplete",
        "bad name",        "bad operator",
        "bad integer",     "bad objective function type",
        "bad bound",       "bad function element",
        "bad constraint"
    };

    return tag[static_cast<int>(failure)];
}

constexpr const char*
problem_definition_error_format(
  baryonyx::problem_definition_error_tag failure) noexcept
{
    constexpr const char* const tag[] = {
        "empty variables",
        "empty objective function",
        "variable not used",
        "bad bound",
        "multiple constraints with different value"
    };

    return tag[static_cast<int>(failure)];
}

constexpr const char*
solver_error_format(baryonyx::solver_error_tag failure) noexcept
{
    constexpr const char* tag[] = { "no solver available",
                                    "unrealisable constraint",
                                    "not enough memory" };

    return tag[static_cast<int>(failure)];
}

#endif
