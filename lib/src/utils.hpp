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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_UTILS_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_UTILS_HPP

#include <baryonyx/core>

#include <chrono>

#include <cmath>
#include <cstdarg>
#include <cstdio>

namespace baryonyx {

inline void
Expects(bool condition)
{
    if (not condition)
        throw precondition_failure("precondition failure");
}

inline void
Ensures(bool condition)
{
    if (not condition)
        throw postcondition_failure("postcondition failure");
}

inline void
Expects(bool condition, const char* s)
{
    if (not condition)
        throw precondition_failure(s);
}

inline void
Ensures(bool condition, const char* s)
{
    if (not condition)
        throw postcondition_failure(s);
}

template<typename T>
inline bool
is_essentially_equal(const T v1, const T v2, const T epsilon)
{
    static_assert(std::is_floating_point<T>::value,
                  "is_essentially_equal required a float/double "
                  "as template arguement");

    return fabs((v1) - (v2)) <=
           ((fabs(v1) > fabs(v2) ? fabs(v2) : fabs(v1)) * (epsilon));
}

std::string
stringf(const char* format, ...) noexcept BARYONYX_FORMAT(1, 2);

std::string
vstringf(const char* format, va_list ap) noexcept;

class timer_profiler
{
public:
    timer_profiler(std::shared_ptr<baryonyx::context> ctx)
      : m_ctx(ctx)
      , m_s(std::chrono::steady_clock::now())
    {
    }

    timer_profiler()
      : m_s(std::chrono::steady_clock::now())
    {
    }

    ~timer_profiler()
    {
        auto end = std::chrono::steady_clock::now();

        if (m_ctx.get())
            m_ctx->info(
              "%fs\n",
              std::chrono::duration_cast<std::chrono::duration<double>>(end -
                                                                        m_s)
                .count());
        else
            fprintf(stderr,
                    "%fs\n",
                    std::chrono::duration_cast<std::chrono::duration<double>>(
                      end - m_s)
                      .count());
    }

private:
    std::shared_ptr<baryonyx::context> m_ctx;
    std::chrono::time_point<std::chrono::steady_clock> m_s;
};

/**
 * @brief Check if the @c limit is reached.
 * @details This function checks if @c end - @c begin is greater than @c limit.
 *
 * @code
 * auto begin = std::chrono::steady_clock::now();
 * // [..] Jbos.
 * if (is_time_limit(10.0, begin, std::chrono::steady_clock::now())
 *    return;
 * @endcode
 *
 * @param limit Duration in second.
 * @param begin Start time.
 * @param end End time or current time.
 */
inline bool
is_time_limit(double limit,
              std::chrono::steady_clock::time_point begin,
              std::chrono::steady_clock::time_point end) noexcept
{
    if (limit <= 0.0)
        return false;

    return std::chrono::duration_cast<std::chrono::duration<double>>(end -
                                                                     begin)
             .count() > limit;
}

} // namespace baryonyx

#endif
