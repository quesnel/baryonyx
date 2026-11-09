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
 * OR IMPLIED, INCLUDING BUT !LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_UTILS_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_UTILS_HPP

#include <baryonyx/core>

#include "private.hpp"

#include <chrono>
#include <concepts>
#include <functional>
#include <iterator>
#include <type_traits>
#include <utility>

#include <climits>
#include <cmath>

namespace baryonyx {

/** Casts nonnegative integer to unsigned.
 */
template<std::integral Integer>
constexpr std::make_unsigned_t<Integer>
to_unsigned(Integer value)
{
    bx_expects(value >= 0);

    return static_cast<std::make_unsigned_t<Integer>>(value);
}

/**
 * @brief Get a sub-string without any @c std::isspace characters at left.
 *
 * @param s The string-view to clean up.
 *
 * @return An update sub-string or the same string-view.
 */
inline std::string_view
left_trim(std::string_view s) noexcept
{
    auto found = s.find_first_not_of(" \t\n\v\f\r");

    if (found == std::string_view::npos)
        return {};

    return s.substr(found, std::string_view::npos);
}

/**
 * @brief Get a sub-string without any @c std::isspace characters at right.
 *
 * @param s The string-view to clean up.
 *
 * @return An update sub-string or the same string-view.
 */
inline std::string_view
right_trim(std::string_view s) noexcept
{
    auto found = s.find_last_not_of(" \t\n\v\f\r");

    if (found == std::string_view::npos)
        return {};

    return s.substr(0, found + 1);
}

/**
 * @brief Get a sub-string without any @c std::isspace characters at left and
 * right
 *
 * @param s The string-view to clean up.
 *
 * @return An update sub-string or the same string-view.
 */
inline std::string_view
trim(std::string_view s) noexcept
{
    return left_trim(right_trim(s));
}

/**
 * @brief Size of a container or a C array as an @c int.
 * @details Baryonyx guarantees that the number of variables and constraints
 *     is lower than @c INT_MAX, so the narrowing is safe.
 */
template<class C>
constexpr int
length(const C& c) noexcept
{
    return static_cast<int>(std::size(c));
}

template<typename Xtype, typename Ytype>
constexpr Ytype
linearize(Xtype p1x, Ytype p1y, Xtype p2x, Ytype p2y, Xtype value) noexcept
{
    bx_expects(p1x - p2x != 0);

    auto m = static_cast<double>(p1y - p2y) / static_cast<double>(p1x - p2x);
    auto p = static_cast<double>(p2y) - m * static_cast<double>(p2x);

    return static_cast<Ytype>(std::ceil(m * static_cast<double>(value) + p));
}

template<std::floating_point T>
inline bool
is_essentially_equal(const T v1, const T v2, const T epsilon)
{
    return std::fabs(v1 - v2) <=
           std::min(std::fabs(v1), std::fabs(v2)) * epsilon;
}

/**
 * @brief Check if the duration between @c begin and @c end is greater than
 * @c limit in second.
 *
 * @details This function checks if @c end - @c begin is greater than @c
 * limit.
 *
 * @code
 * auto begin = std::chrono::steady_clock::now();
 * // computation
 * auto end = std::chrono::steady_clock::now();
 *
 * if (is_time_limit(10.0, begin, end)) {
 *    std::cout << "computation takes more than 10s.\n";
 * }
 * @endcode
 *
 * @param limit Minimal duration in second to test.
 * @param begin Time point from the begining of a compuation.
 * @param end Time point of the current compuation.
 */
inline bool
is_time_limit(double limit,
              std::chrono::steady_clock::time_point begin,
              std::chrono::steady_clock::time_point end) noexcept
{
    using std::chrono::duration;
    using std::chrono::duration_cast;

    if (limit <= 0)
        return false;

    return duration_cast<duration<double>>(end - begin).count() > limit;
}

/**
 * @brief Converts the integer @c s into @c Target.
 *
 * @exception numeric_cast_failure if @c s is out of the @c Target range.
 */
template<std::integral Target, std::integral Source>
constexpr Target
numeric_cast(Source s)
{
    if (!std::in_range<Target>(s))
        throw numeric_cast_failure();

    return static_cast<Target>(s);
}

/**
 * @brief A class to compute time spent during object life.
 *
 * @code
 * {
 *    ...
 *    timer t([](double t){std::cout << t << "s."; });
 *    ...
 * } // Show time spend since timer object instantiation.
 * @endcode
 */
class timer
{
private:
    std::chrono::time_point<std::chrono::steady_clock> m_start;
    std::function<void(double)> m_fct;

public:
    /**
     * @brief Build timer with output function.
     *
     * @param fct Output function called in destructor.
     */
    template<typename Function>
    timer(Function fct)
      : m_start(std::chrono::steady_clock::now())
      , m_fct(fct)
    {}

    double time_elapsed() const
    {
        namespace sc = std::chrono;

        auto diff = std::chrono::steady_clock::now() - m_start;
        auto dc = sc::duration_cast<sc::duration<double, std::ratio<1>>>(diff);

        return dc.count();
    }

    ~timer() noexcept
    {
        try {
            if (m_fct)
                m_fct(time_elapsed());
        } catch (...) {
        }
    }
};

} // namespace baryonyx

#endif
