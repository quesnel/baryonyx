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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_MAIN_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_MAIN_HPP

#include <baryonyx/core>

#include <algorithm>
#include <optional>
#include <string>

#include <fmt/format.h>

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

/**
 * @brief Convert a string_view into a double.
 *
 * @note waiting for std::fron_chars or boost::qi dependencies
 */
constexpr inline std::optional<double>
to_double(std::string_view s) noexcept
{
    // if (auto [p, ec] =
    //       std::from_chars(s.data(), s.data() + s.size(), value, 10);
    //     ec != std::errc())
    //     return bad_value;
    // return value;

    constexpr std::size_t size_limit = 512;

    if (s.size() > size_limit)
        return std::nullopt;

    char buffer[size_limit + 1] = { '\0' };
    std::size_t i = 0;
    std::size_t e = std::min(s.size(), size_limit);

    for (i = 0; i != e; ++i)
        buffer[i] = s[i];

    double result = 0;
    if (auto read = std::sscanf(buffer, "%lf", &result); read)
        return result;
    else
        return std::nullopt;
}

/**
 * @brief Convert a string_view into a integer.
 *
 * @note waiting for std::fron_chars or boost::qi dependencies
 */
constexpr inline std::optional<int>
to_int(std::string_view s) noexcept
{
    constexpr std::size_t size_limit = 512;

    if (s.size() > size_limit)
        return std::nullopt;

    char buffer[size_limit + 1] = { '\0' };
    std::size_t i = 0;
    std::size_t e = std::min(s.size(), size_limit);

    for (i = 0; i != e; ++i)
        buffer[i] = s[i];

    int result = 0;
    if (auto read = std::sscanf(buffer, "%d", &result); read > 0)
        return result;
    else
        return std::nullopt;
}

/**
 * @brief Convert a string_view into a integer.
 *
 * @note waiting for std::fron_chars or boost::qi dependencies
 */
constexpr inline std::optional<long int>
to_long(const std::string_view s) noexcept
{
    constexpr std::size_t size_limit = 512;

    if (s.size() > size_limit)
        return std::nullopt;

    char buffer[size_limit + 1] = { '\0' };
    std::size_t i = 0;
    std::size_t e = std::min(s.size(), size_limit);

    for (i = 0; i != e; ++i)
        buffer[i] = s[i];

    long int result = 0;
    if (auto read = std::sscanf(buffer, "%ld", &result); read > 0)
        return result;
    else
        return std::nullopt;
}

#endif
