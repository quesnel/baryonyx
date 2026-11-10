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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_PRIVATE_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_PRIVATE_HPP

#include <baryonyx/core>

#include <algorithm>
#include <utility>

#include <fmt/color.h>
#include <fmt/format.h>

#include <cstdio>

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
    {
    }

    solver_parameters parameters;

    solver_started_cb start;
    solver_updated_cb update;
    solver_finished_cb finish;

    message_type log_priority = context::message_type::info;

    // Set in the runs of the tuned optimizer, whose candidates must not start
    // a tuned optimization again.
    bool tuning = false;
};

namespace detail {

template<context::message_type Level, typename... Args>
void
log([[maybe_unused]] const context& ctx,
     [[maybe_unused]] fmt::format_string<Args...> fmt,
     [[maybe_unused]] Args&&... args)
{
#ifdef BARYONYX_ENABLE_LOG
#ifndef BARYONYX_ENABLE_DEBUG
    if constexpr (Level == context::message_type::debug)
        return;
#endif

    if (ctx.log_priority < Level)
        return;

    auto* out = Level <= context::message_type::warning ? stderr : stdout;
    fmt::print(out,
               context::message_style[static_cast<int>(Level)],
               fmt,
               std::forward<Args>(args)...);
#endif
}

} // namespace detail

template<typename... Args>
void
emerg(const context& ctx, fmt::format_string<Args...> fmt, Args&&... args)
{
    detail::log<context::message_type::emerg>(
      ctx, fmt, std::forward<Args>(args)...);
}

template<typename... Args>
void
alert(const context& ctx, fmt::format_string<Args...> fmt, Args&&... args)
{
    detail::log<context::message_type::alert>(
      ctx, fmt, std::forward<Args>(args)...);
}

template<typename... Args>
void
crit(const context& ctx, fmt::format_string<Args...> fmt, Args&&... args)
{
    detail::log<context::message_type::crit>(
      ctx, fmt, std::forward<Args>(args)...);
}

template<typename... Args>
void
error(const context& ctx, fmt::format_string<Args...> fmt, Args&&... args)
{
    detail::log<context::message_type::err>(
      ctx, fmt, std::forward<Args>(args)...);
}

template<typename... Args>
void
warning(const context& ctx, fmt::format_string<Args...> fmt, Args&&... args)
{
    detail::log<context::message_type::warning>(
      ctx, fmt, std::forward<Args>(args)...);
}

template<typename... Args>
void
notice(const context& ctx, fmt::format_string<Args...> fmt, Args&&... args)
{
    detail::log<context::message_type::notice>(
      ctx, fmt, std::forward<Args>(args)...);
}

template<typename... Args>
void
info(const context& ctx, fmt::format_string<Args...> fmt, Args&&... args)
{
    detail::log<context::message_type::info>(
      ctx, fmt, std::forward<Args>(args)...);
}

template<typename... Args>
void
debug(const context& ctx, fmt::format_string<Args...> fmt, Args&&... args)
{
    detail::log<context::message_type::debug>(
      ctx, fmt, std::forward<Args>(args)...);
}

} // namespace baryonyx

#endif
