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

#ifndef _WIN32
#include <unistd.h>
#endif

namespace baryonyx {

using string_logger_functor = std::function<void(int, std::string)>;

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

    context(string_logger_functor logger, int verbose_level = 6)
      : string_logger(std::move(logger))
      , cfile_logger(nullptr)
      , logger(context::logger_type::string)
    {
        if (verbose_level != 6)
            log_priority = static_cast<context::message_type>(
              verbose_level < 0 ? 0 : verbose_level > 7 ? 7 : verbose_level);
    }

    solver_parameters parameters;
    std::string method;

    string_logger_functor string_logger;
    FILE* cfile_logger = stdout;
    message_type log_priority = context::message_type::info;
    logger_type logger = context::logger_type::c_file;
    bool color_cfile_logger = false;
};

/**
 * @brief Shows value of all baryonyx parameters.
 *
 */
void
print(const context_ptr& ctx);

template<typename... Args>
void
log(const context_ptr& ctx,
    context::message_type level,
    const char* fmt,
    const Args&... args);

template<typename... Args>
void
log(const context_ptr& ctx, context::message_type level, const char* msg);

template<typename... Args>
void
log(baryonyx::context* ctx,
    context::message_type level,
    const char* fmt,
    const Args&... args);

template<typename... Args>
void
log(baryonyx::context* ctx, context::message_type level, const char* msg);

template<typename T>
void
log(const context_ptr& ctx, context::message_type level, const T& msg);

template<typename T>
void
log(baryonyx::context* ctx, context::message_type level, const T& msg);

baryonyx::problem
read_problem(std::istream& is);

baryonyx::result
read_result(std::istream& is);

bool
write_problem(std::ostream& os, const baryonyx::problem& pb);

bool
check_consistency(const baryonyx::problem& pb);

baryonyx::result
solver_select(const baryonyx::context_ptr& ctx, baryonyx::problem& pb);

baryonyx::result
optimizer_select(const baryonyx::context_ptr& ctx, baryonyx::problem& pb);

} // namespace baryonyx

#endif
