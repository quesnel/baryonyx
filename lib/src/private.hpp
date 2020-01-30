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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_PRIVATE_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_PRIVATE_HPP

#include <baryonyx/core>

#include <utility>

#include <fmt/color.h>
#include <fmt/format.h>

#include <cstdio>

#define bx_stringify_detail(x) #x
#define bx_stringify(x) bx_stringify_detail(x)

#if defined(__clang__) || defined(__GNUC__)
#define bx_always_inline __attribute__((always_inline))
#define bx_likely(x) __builtin_expect(!!(x), 1)
#define bx_unlikely(x) __builtin_expect(!!(x), 0)
#else
#define bx_always_inline
#define bx_likely(x) (!!(x))
#define bx_unlikely(x) (!!(x))
#endif

namespace baryonyx {

struct context
{
    enum class message_type
    {
        emerg = 0, ///< system is unusable
        alert,     ///< action must be taken immediately
        crit,      ///< critical conditions
        err,       ///< error conditions
        warning,   ///< warning conditions
        notice,    ///< normal, but significant, condition
        info,      ///< informational message
        debug,     ///< debug-level message
    };

    static const fmt::text_style message_style[];

    static message_type convert_verbose_level(int level) noexcept
    {
        return static_cast<context::message_type>(std::clamp(level, 0, 7));
    }

    context(int verbose_level = 6)
      : log_priority(convert_verbose_level(verbose_level))
    {}

    solver_parameters parameters;
    std::string method;

    solver_started_cb start;
    solver_updated_cb update;
    solver_finished_cb finish;

    message_type log_priority = context::message_type::info;
};

template<typename... Args>
void
to_log([[maybe_unused]] std::FILE* os,
       [[maybe_unused]] const std::string_view fmt,
       [[maybe_unused]] const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
#ifdef BARYONYX_ENABLE_DEBUG
    fmt::print(os, fmt, args...);
#endif
#endif
}

template<typename... Args>
void
to_log([[maybe_unused]] std::FILE* os,
       [[maybe_unused]] unsigned indent,
       [[maybe_unused]] const std::string_view fmt,
       [[maybe_unused]] const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
#ifdef BARYONYX_ENABLE_DEBUG
    fmt::print(os, "{:{}}", "", indent);
    fmt::print(os, fmt, args...);
#endif
#endif
}

inline bool
is_loggable(context::message_type current_level, int level) noexcept
{
    return static_cast<int>(current_level) >= level;
}

template<typename... Args>
void
log(FILE* out, int style, const char* fmt, const Args&... args)
{
    fmt::print(out, context::message_style[style], fmt, args...);
}

template<typename... Args>
void
log(FILE* out, int style, const char* msg)
{
    fmt::print(out, context::message_style[style], msg);
}

template<typename T>
void
log(FILE* out, int style, const T& msg)
{
    fmt::print(out, context::message_style[style], "{}", msg);
}

namespace detail {

struct sink
{
    template<typename... Args>
    sink(const Args&...)
    {}
};

}

template<typename... Args>
void
emerg(const context& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 0))
        return;

    log(stderr, 0, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
}

template<typename... Args>
void
alert(const context& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 1))
        return;

    log(stderr, 1, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
}

template<typename... Args>
void
crit(const context& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 2))
        return;

    log(stderr, 2, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
}

template<typename... Args>
void
error(const context& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 3))
        return;

    log(stderr, 3, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
}

template<typename... Args>
void
warning(const context& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 4))
        return;

    log(stderr, 4, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
}

template<typename... Args>
void
notice(const context& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 5))
        return;

    log(stdout, 5, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
}

template<typename... Args>
void
info(const context& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 6))
        return;

    log(stdout, 6, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
}

template<typename... Args>
void
debug(const context& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
#ifdef BARYONYX_ENABLE_DEBUG
    if (!is_loggable(ctx.log_priority, 7))
        return;

    log(stdout, 7, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
#else
    detail::sink(ctx, fmt, args...);
#endif
}

template<typename Arg1, typename... Args>
void
emerg(const context& ctx,
      const char* fmt,
      const Arg1& arg1,
      const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 0))
        return;

    log(stderr, 0, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
}

template<typename Arg1, typename... Args>
void
alert(const context& ctx,
      const char* fmt,
      const Arg1& arg1,
      const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 1))
        return;

    log(stderr, 1, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
}

template<typename Arg1, typename... Args>
void
crit(const context& ctx,
     const char* fmt,
     const Arg1& arg1,
     const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 2))
        return;

    log(stderr, 2, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
}

template<typename Arg1, typename... Args>
void
error(const context& ctx,
      const char* fmt,
      const Arg1& arg1,
      const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 3))
        return;

    log(stderr, 3, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
}

template<typename Arg1, typename... Args>
void
warning(const context& ctx,
        const char* fmt,
        const Arg1& arg1,
        const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 4))
        return;

    log(stderr, 4, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
}

template<typename Arg1, typename... Args>
void
notice(const context& ctx,
       const char* fmt,
       const Arg1& arg1,
       const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 5))
        return;

    log(stdout, 5, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
}

template<typename Arg1, typename... Args>
void
info(const context& ctx,
     const char* fmt,
     const Arg1& arg1,
     const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 6))
        return;

    log(stdout, 6, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
}

template<typename Arg1, typename... Args>
void
debug(const context& ctx,
      const char* fmt,
      const Arg1& arg1,
      const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
#ifdef BARYONYX_ENABLE_DEBUG
    if (!is_loggable(ctx.log_priority, 7))
        return;

    log(stdout, 7, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
}

////////////////////////////////////////////////

template<typename T>
void
emerg(const context& ctx, const T& msg)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 0))
        return;

    log(stderr, 0, msg);
#else
    detail::sink(ctx, msg);
#endif
}

template<typename T>
void
alert(const context& ctx, const T& msg)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 1))
        return;

    log(stderr, 1, msg);
#else
    detail::sink(ctx, msg);
#endif
}

template<typename T>
void
crit(const context& ctx, const T& msg)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 2))
        return;

    log(stderr, 2, msg);
#else
    detail::sink(ctx, msg);
#endif
}

template<typename T>
void
error(const context& ctx, const T& msg)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 3))
        return;

    log(stderr, 3, msg);
#else
    detail::sink(ctx, msg);
#endif
}

template<typename T>
void
warning(const context& ctx, const T& msg)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 4))
        return;

    log(stderr, 4, msg);
#else
    detail::sink(ctx, msg);
#endif
}

template<typename T>
void
notice(const context& ctx, const T& msg)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 5))
        return;

    log(stdout, 5, msg);
#else
    detail::sink(ctx, msg);
#endif
}

template<typename T>
void
info(const context& ctx, const T& msg)
{
#ifdef BARYONYX_ENABLE_LOG
    if (!is_loggable(ctx.log_priority, 6))
        return;

    log(stdout, 6, msg);
#else
    detail::sink(ctx, msg);
#endif
}

template<typename T>
void
debug(const context& ctx, const T& msg)
{
#ifdef BARYONYX_ENABLE_LOG
#ifndef BARYONYX_ENABLE_DEBUG
    if (!is_loggable(ctx.log_priority, 7))
        return;

    log(stdout, 7, msg);
#else
    detail::sink(ctx, msg);
#endif
#else
    detail::sink(ctx, msg);
#endif
}

} // namespace baryonyx

#endif
