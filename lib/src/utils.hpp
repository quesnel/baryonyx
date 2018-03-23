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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_UTILS_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_UTILS_HPP

#include <baryonyx/core>

#include "private.hpp"

#include <fmt/format.h>
#include <fmt/printf.h>

#include <chrono>

#include <cassert>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <utility>

namespace baryonyx {

namespace details {

enum class colors
{
    Default = 0,
    Black,
    Red,
    Green,
    Yellow,
    Blue,
    Magenta,
    Cyan,
    Light_gray,
    Dark_gray,
    Light_red,
    Light_green,
    Light_yellow,
    Light_blue,
    Light_magenta,
    Light_cyan,
    White,
    No_color_change
};

enum class setters
{
    Reset = 0,
    Bold,
    Dim,
    Underlined,
    No_setter_change
};

static const char colors_str[][10] = {
    "\033[39m", "\033[30m", "\033[31m", "\033[32m", "\033[33m", "\033[34m",
    "\033[35m", "\033[36m", "\033[37m", "\033[90m", "\033[91m", "\033[92m",
    "\033[93m", "\033[94m", "\033[95m", "\033[96m", "\033[97m"
};

static const char setters_str[][10] = { "\033[0m",
                                        "\033[1m",
                                        "\033[2m",
                                        "\033[4m" };

constexpr const char*
color_to_str(baryonyx::details::colors c) noexcept
{
    assert(static_cast<int>(c) < static_cast<int>(array_lenght(colors_str)));

    return colors_str[static_cast<int>(c)];
}

constexpr inline const char*
setter_to_str(baryonyx::details::setters s) noexcept
{

    assert(static_cast<int>(s) < static_cast<int>(array_lenght(setters_str)));

    return setters_str[static_cast<int>(s)];
}

struct log_color
{
    FILE* f;

    constexpr log_color(FILE* f_, context::message_type type)
      : f(f_)
    {
        switch (type) {
        case context::message_type::emerg:
            fmt::print(f, color_to_str(colors::Red));
            fmt::print(f, setter_to_str(setters::Bold));
            break;
        case context::message_type::alert:
            fmt::print(f, color_to_str(colors::Red));
            fmt::print(f, setter_to_str(setters::Reset));
            break;
        case context::message_type::crit:
            fmt::print(f, color_to_str(colors::Magenta));
            fmt::print(f, setter_to_str(setters::Bold));
            break;
        case context::message_type::err:
            fmt::print(f, color_to_str(colors::Magenta));
            fmt::print(f, setter_to_str(setters::Reset));
            break;
        case context::message_type::warning:
            fmt::print(f, color_to_str(colors::Yellow));
            fmt::print(f, setter_to_str(setters::Bold));
            break;
        case context::message_type::notice:
            fmt::print(f, color_to_str(colors::Yellow));
            fmt::print(f, setter_to_str(setters::Underlined));
            break;
        case context::message_type::info:
            fmt::print(f, color_to_str(colors::Default));
            fmt::print(f, setter_to_str(setters::Reset));
            break;
        case context::message_type::debug:
            fmt::print(f, color_to_str(colors::Dark_gray));
            fmt::print(f, setter_to_str(setters::Reset));
            break;
        }
    }

    ~log_color()
    {
        try {
            fmt::print(f, color_to_str(colors::Default));
            fmt::print(f, setter_to_str(setters::Reset));
        } catch (...) {
        }
    }
};

} // namespace details

inline bool
is_loggable(context::message_type current_level,
            context::message_type level) noexcept
{
    return static_cast<int>(current_level) >= static_cast<int>(level);
}

template<typename... Args>
void
log(const context_ptr& ctx,
    context::message_type level,
    const char* fmt,
    const Args&... args)
{
    if (not is_loggable(ctx->log_priority, level))
        return;

    if (ctx->logger == context::logger_type::c_file) {
        if (ctx->color_cfile_logger) {
            details::log_color lc(ctx->cfile_logger, level);
            fmt::print(ctx->cfile_logger, fmt, args...);
        } else {
            fmt::print(ctx->cfile_logger, fmt, args...);
        }
    } else {
        ctx->string_logger(static_cast<int>(level), fmt::format(fmt, args...));
    }
}

template<typename... Args>
void
log(const context_ptr& ctx, context::message_type level, const char* msg)
{
    if (not is_loggable(ctx->log_priority, level))
        return;

    if (ctx->logger == context::logger_type::c_file) {
        if (ctx->color_cfile_logger) {
            details::log_color lc(ctx->cfile_logger, level);
            fmt::print(ctx->cfile_logger, msg);
        } else {
            fmt::print(ctx->cfile_logger, msg);
        }
    } else {
        ctx->string_logger(static_cast<int>(level), fmt::format(msg));
    }
}

template<typename... Args>
void
log(baryonyx::context* ctx,
    context::message_type level,
    const char* fmt,
    const Args&... args)
{
    if (not is_loggable(ctx->log_priority, level))
        return;

    if (ctx->logger == context::logger_type::c_file) {
        if (ctx->color_cfile_logger) {
            details::log_color lc(ctx->cfile_logger, level);
            fmt::print(ctx->cfile_logger, fmt, args...);
        } else {
            fmt::print(ctx->cfile_logger, fmt, args...);
        }
    } else {
        ctx->string_logger(static_cast<int>(level), fmt::format(fmt, args...));
    }
}

template<typename... Args>
void
log(baryonyx::context* ctx, context::message_type level, const char* msg)
{
    if (not is_loggable(ctx->log_priority, level))
        return;

    if (ctx->logger == context::logger_type::c_file) {
        if (ctx->color_cfile_logger) {
            details::log_color lc(ctx->cfile_logger, level);
            fmt::print(ctx->cfile_logger, msg);
        } else {
            fmt::print(ctx->cfile_logger, msg);
        }
    } else {
        ctx->string_logger(static_cast<int>(level), fmt::format(msg));
    }
}

template<typename... Args>
void
emerg(const context_ptr& ctx, const char* fmt, const Args&... args)
{
    log(ctx, context::message_type::emerg, fmt, args...);
}

template<typename... Args>
void
alert(const context_ptr& ctx, const char* fmt, const Args&... args)
{
    log(ctx, context::message_type::alert, fmt, args...);
}

template<typename... Args>
void
crit(const context_ptr& ctx, const char* fmt, const Args&... args)
{
    log(ctx, context::message_type::crit, fmt, args...);
}

template<typename... Args>
void
error(const context_ptr& ctx, const char* fmt, const Args&... args)
{
    log(ctx, context::message_type::err, fmt, args...);
}

template<typename... Args>
void
warning(const context_ptr& ctx, const char* fmt, const Args&... args)
{
    log(ctx, context::message_type::warning, fmt, args...);
}

template<typename... Args>
void
notice(const context_ptr& ctx, const char* fmt, const Args&... args)
{
    log(ctx, context::message_type::notice, fmt, args...);
}

template<typename... Args>
void
info(const context_ptr& ctx, const char* fmt, const Args&... args)
{
    log(ctx, context::message_type::info, fmt, args...);
}

template<typename... Args>
void
debug(const context_ptr& ctx, const char* fmt, const Args&... args)
{
#ifndef BARYONYX_DISABLE_LOGGING
    //
    // Default, the logging system is active and the call to the @c log
    // function are send to the logger functor. Define BARYONYX_DISABLE_LOGGING
    // as preprocessor value to hide all logging message..
    //
    log(ctx, context::message_type::debug, fmt, args...);
#else
    (void)ctx;
    (void)fmt;
    (void)arg1;
#endif
}

template<typename Arg1, typename... Args>
void
emerg(const context_ptr& ctx,
      const char* fmt,
      const Arg1& arg1,
      const Args&... args)
{
    log(ctx, context::message_type::emerg, fmt, arg1, args...);
}

template<typename Arg1, typename... Args>
void
alert(const context_ptr& ctx,
      const char* fmt,
      const Arg1& arg1,
      const Args&... args)
{
    log(ctx, context::message_type::alert, fmt, arg1, args...);
}

template<typename Arg1, typename... Args>
void
crit(const context_ptr& ctx,
     const char* fmt,
     const Arg1& arg1,
     const Args&... args)
{
    log(ctx, context::message_type::crit, fmt, arg1, args...);
}

template<typename Arg1, typename... Args>
void
error(const context_ptr& ctx,
      const char* fmt,
      const Arg1& arg1,
      const Args&... args)
{
    log(ctx, context::message_type::err, fmt, arg1, args...);
}

template<typename Arg1, typename... Args>
void
warning(const context_ptr& ctx,
        const char* fmt,
        const Arg1& arg1,
        const Args&... args)
{
    log(ctx, context::message_type::warning, fmt, arg1, args...);
}

template<typename Arg1, typename... Args>
void
notice(const context_ptr& ctx,
       const char* fmt,
       const Arg1& arg1,
       const Args&... args)
{
    log(ctx, context::message_type::notice, fmt, arg1, args...);
}

template<typename Arg1, typename... Args>
void
info(const context_ptr& ctx,
     const char* fmt,
     const Arg1& arg1,
     const Args&... args)
{
    log(ctx, context::message_type::info, fmt, arg1, args...);
}

template<typename Arg1, typename... Args>
void
debug(const context_ptr& ctx,
      const char* fmt,
      const Arg1& arg1,
      const Args&... args)
{
#ifndef BARYONYX_DISABLE_LOGGING
    //
    // Default, the logging system is active and the call to the @c log
    // function are send to the logger functor. Define BARYONYX_DISABLE_LOGGING
    // as preprocessor value to hide all logging message..
    //
    log(ctx, context::message_type::debug, fmt, arg1, args...);
#else
    (void)ctx;
    (void)fmt;
    (void)arg1;
    (void)args;
#endif
}

template<typename T>
void
log(const context_ptr& ctx, context::message_type level, const T& msg)
{
    if (not is_loggable(ctx->log_priority, level))
        return;

    if (ctx->logger == context::logger_type::c_file) {
        if (ctx->color_cfile_logger) {
            details::log_color lc(ctx->cfile_logger, level);
            fmt::print(ctx->cfile_logger, "{}", msg);
        } else {
            fmt::print(ctx->cfile_logger, "{}", msg);
        }
    } else {
        ctx->string_logger(static_cast<int>(level), fmt::format("{}", msg));
    }
}

template<typename T>
void
log(baryonyx::context* ctx, context::message_type level, const T& msg)
{
    if (not is_loggable(ctx->log_priority, level))
        return;

    if (ctx->logger == context::logger_type::c_file) {
        if (ctx->color_cfile_logger) {
            details::log_color lc(ctx->cfile_logger, level);
            fmt::print(ctx->cfile_logger, "{}", msg);
        } else {
            fmt::print(ctx->cfile_logger, "{}", msg);
        }
    } else {
        ctx->string_logger(static_cast<int>(level), fmt::format("{}", msg));
    }
}

////////////////////////////////////////////////

template<typename T>
void
emerg(const context_ptr& ctx, const T& msg)
{
    log(ctx, context::message_type::emerg, msg);
}

template<typename T>
void
alert(const context_ptr& ctx, const T& msg)
{
    log(ctx, context::message_type::alert, msg);
}

template<typename T>
void
crit(const context_ptr& ctx, const T& msg)
{
    log(ctx, context::message_type::crit, msg);
}

template<typename T>
void
error(const context_ptr& ctx, const T& msg)
{
    log(ctx, context::message_type::err, msg);
}

template<typename T>
void
warning(const context_ptr& ctx, const T& msg)
{
    log(ctx, context::message_type::warning, msg);
}

template<typename T>
void
notice(const context_ptr& ctx, const T& msg)
{
    log(ctx, context::message_type::notice, msg);
}

template<typename T>
void
info(const context_ptr& ctx, const T& msg)
{
    log(ctx, context::message_type::info, msg);
}

template<typename T>
void
debug(const context_ptr& ctx, const T& msg)
{
#ifndef BARYONYX_DISABLE_LOGGING
    //
    // Default, the logging system is active and the call to the @c log
    // function are send to the logger functor. Define BARYONYX_DISABLE_LOGGING
    // as preprocessor value to hide all logging message..
    //
    log(ctx, context::message_type::debug, msg);
#else
    (void)ctx;
    (void)msg;
#endif
}

////////////////////////////////////////////////

/**
 * @brief Compute the length of the @c container.
 * @details Return the @c size provided by the @c C::size() but cast it into a
 *     @c int. This is a specific baryonyx function, we know that number of
 *         variables and constraints are lower than the  @c int max value
 *         (INT_MAX).
 *
 * @code
 * std::vector<int> v(z);
 *
 * for (int i = 0, e = length(v); i != e; ++i)
 *     ...
 * @endcode
 *
 * @param c The container to request to size.
 * @tparam T The of container value (must provide a @c size() member).
 */
template<class C>
constexpr int
length(const C& c) noexcept
{
    assert(c.size() <= static_cast<std::size_t>(INT_MAX));

    return static_cast<int>(c.size());
}

/**
 * @details Reference to @c lo if @c v is less than @c lo, reference to @c hi
 *     if @c hi is less than @c v, otherwise reference to @c v.
 *
 * @code
 * assert(baryonyx::clamp(0.0, 0.0, 1.0) == 0.0);
 * assert(baryonyx::clamp(1.0, 0.0, 1.0) == 1.0);
 * assert(baryonyx::clamp(-0.5, 0.0, 1.0) == 0.0);
 * assert(baryonyx::clamp(1.5, 0.0, 1.0) == 1.0);
 * assert(baryonyx::clamp(168, -128, +127) == 127);
 * assert(baryonyx::clamp(168, 0, +255) == 168);
 * assert(baryonyx::clamp(128, -128, +127) == 127);
 * assert(baryonyx::clamp(128, 0, +255) == 128);
 * @endcode
 *
 * @param v The value to clamp.
 * @param lo The low value to compare.
 * @param hi The high value to compare.
 * @return value, low or high.
 */
template<class T>
constexpr const T&
clamp(const T& v, const T& lo, const T& hi)
{
    assert(lo < hi);

    return v < lo ? lo : v > hi ? hi : v;
}

/**
 * @brief Compute the length of the C array.
 * @details Return the size of the C array but cast it into a @c int. This is a
 *     specific baryonyx function, we know that number of variables and
 *     constraints are lower than the  @c int max value (INT_MAX).
 *
 * @code
 * int v[150];
 * for (int i = 0, e = length(v); i != e; ++i)
 *     ...
 * @endcode
 *
 * @param v The container to return size.
 * @tparam T The type of the C array.
 * @tparam N The size of the C array.
 */
template<class T, size_t N>
constexpr int
length(const T (&array)[N]) noexcept
{
    (void)array;

    assert(N <= INT_MAX);

    return static_cast<int>(N);
}

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

/**
 * @brief Check if the duration between @c begin and @c end is greater than @c
 *     limit in second.
 *
 * @details This function checks if @c end - @c begin is greater than @c limit.
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

} // namespace baryonyx

#endif
