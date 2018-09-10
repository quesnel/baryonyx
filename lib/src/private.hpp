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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_PRIVATE_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_PRIVATE_HPP

#include <baryonyx/core>

#include <unordered_map>
#include <utility>

#include <fmt/format.h>

#include <cstdio>

#ifndef _WIN32
#include <unistd.h>
#endif

#define bx_stringify_detail(x) #x
#define bx_stringify(x) bx_stringify_detail(x)

#if defined(__clang__) || defined(__GNUC__)
#define bx_always_inline __attribute__((always_inline))
#define bx_likely(x) __builtin_expect(!!(x), 1)
#define bx_unlikely(x) __builtin_expect(!!(x), 0)
#else
#define bx_always_inline
#define likely(x) (!!(x))
#define unlikely(x) (!!(x))
#endif

namespace baryonyx {

using string_logger_functor = std::function<void(int, std::string)>;

inline const char*
to_string(solver_parameters::pre_constraint_order type) noexcept

{
    static const char* ret[] = { "none",
                                 "reversing",
                                 "random-sorting",
                                 "infeasibility-decr",
                                 "infeasibility-incr" };

    return ret[static_cast<int>(type)];
}

// Context

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

    enum class logger_type
    {
        c_file, ///< log are send to a C FILE* structure.
        string  ///< log are store to the string_logger_functor.
    };

    context(int verbose_level = 6)
      : cfile_logger(stdout)
    {
        if (verbose_level != 6)
            log_priority = static_cast<context::message_type>(
              verbose_level < 0 ? 0 : verbose_level > 7 ? 7 : verbose_level);

#ifndef _WIN32
        if (::isatty(::fileno(cfile_logger)))
            color_cfile_logger = true;
#endif
    }

    context(FILE* f, int verbose_level = 6)
      : cfile_logger(f ? f : stdout)
    {
        if (verbose_level != 6)
            log_priority = static_cast<context::message_type>(
              verbose_level < 0 ? 0 : verbose_level > 7 ? 7 : verbose_level);

#ifndef _WIN32
        if (::isatty(::fileno(cfile_logger)))
            color_cfile_logger = true;
#endif
    }

    context(string_logger_functor logger_, int verbose_level_ = 6)
      : string_logger(std::move(logger_))
      , cfile_logger(nullptr)
      , logger(context::logger_type::string)
    {
        if (verbose_level_ != 6)
            log_priority = static_cast<context::message_type>(
              verbose_level_ < 0 ? 0
                                 : verbose_level_ > 7 ? 7 : verbose_level_);
    }

    solver_parameters parameters;
    std::string method;

    string_logger_functor string_logger;
    FILE* cfile_logger = stdout;
    message_type log_priority = context::message_type::info;
    logger_type logger = context::logger_type::c_file;
    bool color_cfile_logger = false;
};

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

#ifdef __unix__
#define error_symbol "\u26d4 "
#define warning_symbol "\u26a0 "
#endif

constexpr const char*
color_to_str(baryonyx::details::colors c) noexcept
{
    return colors_str[static_cast<int>(c)];
}

constexpr inline const char*
setter_to_str(baryonyx::details::setters s) noexcept
{
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
#ifdef __unix__
            fmt::print(f, error_symbol);
#endif
            break;
        case context::message_type::alert:
            fmt::print(f, color_to_str(colors::Red));
#ifdef __unix__
            fmt::print(f, error_symbol);
#endif
            break;
        case context::message_type::crit:
            fmt::print(f, color_to_str(colors::Magenta));
            fmt::print(f, setter_to_str(setters::Bold));
#ifdef __unix__
            fmt::print(f, error_symbol);
#endif
            break;
        case context::message_type::err:
            fmt::print(f, color_to_str(colors::Magenta));
#ifdef __unix__
            fmt::print(f, error_symbol);
#endif
            break;
        case context::message_type::warning:
            fmt::print(f, color_to_str(colors::Yellow));
            fmt::print(f, setter_to_str(setters::Bold));
#ifdef __unix__
            fmt::print(f, warning_symbol);
#endif
            break;
        case context::message_type::notice:
            fmt::print(f, color_to_str(colors::Yellow));
            break;
        case context::message_type::info:
            fmt::print(f, color_to_str(colors::Default));
            break;
        case context::message_type::debug:
            fmt::print(f, color_to_str(colors::Dark_gray));
            break;
        }
    }

    ~log_color() noexcept
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
    if (!is_loggable(ctx->log_priority, level))
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
    if (!is_loggable(ctx->log_priority, level))
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
emerg(const context_ptr& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::emerg, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
}

template<typename... Args>
void
alert(const context_ptr& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::alert, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
}

template<typename... Args>
void
crit(const context_ptr& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::crit, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
}

template<typename... Args>
void
error(const context_ptr& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::err, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
}

template<typename... Args>
void
warning(const context_ptr& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::warning, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
}

template<typename... Args>
void
notice(const context_ptr& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::notice, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
}

template<typename... Args>
void
info(const context_ptr& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::info, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
}

template<typename... Args>
void
debug(const context_ptr& ctx, const char* fmt, const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
#ifdef BARYONYX_ENABLE_DEBUG
    log(ctx, context::message_type::debug, fmt, args...);
#else
    detail::sink(ctx, fmt, args...);
#endif
#endif
}

template<typename Arg1, typename... Args>
void
emerg(const context_ptr& ctx,
      const char* fmt,
      const Arg1& arg1,
      const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::emerg, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
}

template<typename Arg1, typename... Args>
void
alert(const context_ptr& ctx,
      const char* fmt,
      const Arg1& arg1,
      const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::alert, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
}

template<typename Arg1, typename... Args>
void
crit(const context_ptr& ctx,
     const char* fmt,
     const Arg1& arg1,
     const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::crit, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
}

template<typename Arg1, typename... Args>
void
error(const context_ptr& ctx,
      const char* fmt,
      const Arg1& arg1,
      const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::err, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
}

template<typename Arg1, typename... Args>
void
warning(const context_ptr& ctx,
        const char* fmt,
        const Arg1& arg1,
        const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::warning, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
}

template<typename Arg1, typename... Args>
void
notice(const context_ptr& ctx,
       const char* fmt,
       const Arg1& arg1,
       const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::notice, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
}

template<typename Arg1, typename... Args>
void
info(const context_ptr& ctx,
     const char* fmt,
     const Arg1& arg1,
     const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::info, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
}

template<typename Arg1, typename... Args>
void
debug(const context_ptr& ctx,
      const char* fmt,
      const Arg1& arg1,
      const Args&... args)
{
#ifdef BARYONYX_ENABLE_LOG
#ifdef BARYONYX_ENABLE_DEBUG
    log(ctx, context::message_type::debug, fmt, arg1, args...);
#else
    detail::sink(ctx, fmt, arg1, args...);
#endif
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
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::emerg, msg);
#else
    detail::sink(ctx, msg);
#endif
}

template<typename T>
void
alert(const context_ptr& ctx, const T& msg)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::alert, msg);
#else
    detail::sink(ctx, msg);
#endif
}

template<typename T>
void
crit(const context_ptr& ctx, const T& msg)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::crit, msg);
#else
    detail::sink(ctx, msg);
#endif
}

template<typename T>
void
error(const context_ptr& ctx, const T& msg)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::err, msg);
#else
    detail::sink(ctx, msg);
#endif
}

template<typename T>
void
warning(const context_ptr& ctx, const T& msg)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::warning, msg);
#else
    detail::sink(ctx, msg);
#endif
}

template<typename T>
void
notice(const context_ptr& ctx, const T& msg)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::notice, msg);
#else
    detail::sink(ctx, msg);
#endif
}

template<typename T>
void
info(const context_ptr& ctx, const T& msg)
{
#ifdef BARYONYX_ENABLE_LOG
    log(ctx, context::message_type::info, msg);
#else
    detail::sink(ctx, msg);
#endif
}

template<typename T>
void
debug(const context_ptr& ctx, const T& msg)
{
#ifdef BARYONYX_ENABLE_LOG
#ifndef BARYONYX_ENABLE_DEBUG
    log(ctx, context::message_type::debug, msg);
#else
    detail::sink(ctx, msg);
#endif
#endif
}

/**
 * @brief Shows value of all baryonyx parameters.
 *
 */
void
print(const context_ptr& ctx);

} // namespace baryonyx

#endif
