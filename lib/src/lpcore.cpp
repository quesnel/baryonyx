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

#include "lpformat-consistency.hpp"
#include "lpformat-io.hpp"
#include "mitm.hpp"

#include <algorithm>
#include <fstream>
#include <lpcore>

#ifdef __unix__
#include <unistd.h>
#endif

namespace lp {

class standard_stream_logger : public context::logger
{
public:
    void write(int priority,
               const char* file,
               int line,
               const char* fn,
               const char* format,
               va_list args) noexcept override
    {
        if (priority > 5)
            vfprintf(stdout, format, args);
        else {
            fprintf(stderr,
                    "lp: %d at %d in function '%s' from file %s: ",
                    priority,
                    line,
                    fn,
                    file);
            vfprintf(stderr, format, args);
        }
    }

    void write(lp::context::message_type m,
               const char* format,
               va_list args) noexcept override
    {
#ifdef __unix__
        if (::isatty(STDOUT_FILENO)) {
            switch (m) {
            case context::message_type::emerg:
            case context::message_type::alert:
            case context::message_type::crit:
            case context::message_type::err:
                ::puts("\033[30m\033[2m");
                break;
            case context::message_type::warning:
                ::puts("\033[32m\033[1m");
                break;
            case context::message_type::notice:
                break;
            case context::message_type::info:
                break;
            case context::message_type::debug:
                ::puts("\033[33m\033[1m");
                break;
            }

            vfprintf(stdout, format, args);

            switch (m) {
            case context::message_type::emerg:
            case context::message_type::alert:
            case context::message_type::crit:
            case context::message_type::err:
                ::puts("\033[30m\033[0m");
                break;
            case context::message_type::warning:
                ::puts("\033[30m\033[0m");
                break;
            case context::message_type::notice:
                break;
            case context::message_type::info:
                break;
            case context::message_type::debug:
                ::puts("\033[30m\033[0m");
                break;
            }

        } else {
            vfprintf(stdout, format, args);
        }
#else
        vfprintf(stdout, format, args);
#endif
    }
};

void
context::set_log_priority(int priority) noexcept
{
    m_log_priority = priority < 0 ? 0 : 7 > priority ? 7 : priority;
}

int
context::get_log_priority() const noexcept
{
    return m_log_priority;
}

void
context::set_standard_stream_logger() noexcept
{
    set_logger(std::make_unique<standard_stream_logger>());
}

void
context::set_logger(std::unique_ptr<logger> function) noexcept
{
    m_logger = std::move(function);
}

#ifndef LP_DISABLE_LOGGING
void
context::log(message_type type, const char* format, ...) noexcept
{
    if (not m_logger)
        return;

    va_list args;

    va_start(args, format);
    m_logger->write(type, format, args);
    va_end(args);
}

void
context::log(int priority,
             const char* file,
             int line,
             const char* fn,
             const char* format,
             ...) noexcept
{
    if (not m_logger)
        return;

    va_list args;

    va_start(args, format);
    m_logger->write(priority, file, line, fn, format, args);
    va_end(args);
}

void
context::info(const char* format, ...) noexcept
{
    if (not m_logger)
        return;

    va_list args;

    va_start(args, format);
    m_logger->write(context::message_type::info, format, args);
    va_end(args);
}

void
context::debug(const char* format, ...) noexcept
{
    if (not m_logger)
        return;

    va_list args;

    va_start(args, format);
    m_logger->write(context::message_type::debug, format, args);
    va_end(args);
}

void
context::warning(const char* format, ...) noexcept
{
    if (not m_logger)
        return;

    va_list args;

    va_start(args, format);
    m_logger->write(context::message_type::warning, format, args);
    va_end(args);
}

void
context::error(const char* format, ...) noexcept
{
    if (not m_logger)
        return;

    va_list args;

    va_start(args, format);
    m_logger->write(context::message_type::err, format, args);
    va_end(args);
}

#endif

problem
make_problem(std::shared_ptr<lp::context> ctx, const std::string& filename)
{
    ctx->info("problem read from file `%s'\n", filename.c_str());

    std::ifstream ifs;
    ifs.exceptions(std::ifstream::badbit);
    ifs.open(filename);

    return details::read_problem(ifs);
}

problem
make_problem(std::shared_ptr<lp::context> ctx, std::istream& is)
{
    ctx->info("problem read from stream\n");

    is.exceptions(std::ifstream::badbit);

    return details::read_problem(is);
}

std::ostream&
operator<<(std::ostream& os, const problem& p)
{
    details::problem_writer pw(p, os);

    return os;
}

result
solve(std::shared_ptr<lp::context> ctx, problem& pb)
{
    check(pb);

    std::map<std::string, parameter> params;

    return mitm_solve(ctx, pb, params);
}

result
solve(std::shared_ptr<lp::context> ctx,
      problem& pb,
      const std::map<std::string, parameter>& params)
{
    check(pb);

    return mitm_solve(ctx, pb, params);
}

result
optimize(std::shared_ptr<lp::context> ctx,
         problem& pb,
         const std::map<std::string, parameter>& params)
{
    check(pb);

    return mitm_optimize(ctx, pb, params);
}

result
optimize(std::shared_ptr<lp::context> ctx, problem& pb)
{
    check(pb);

    std::map<std::string, parameter> params;

    return mitm_optimize(ctx, pb, params);
}

template <typename functionT, typename variablesT>
int
compute_function(const functionT& fct, const variablesT& vars) noexcept
{
    int v{ 0 };

    for (auto& f : fct)
        v += f.factor * vars[f.variable_index];

    return v;
}

bool
is_valid_solution(const problem& pb, const std::vector<int>& variable_value)
{
    Expects(not variable_value.empty(), "variables vector empty");

    for (auto& cst : pb.equal_constraints) {
        if (compute_function(cst.elements, variable_value) != cst.value) {
            printf("constraint %s (=) fails\n", cst.label.c_str());
            return false;
        }
    }

    for (auto& cst : pb.greater_constraints) {
        if (compute_function(cst.elements, variable_value) <= cst.value) {
            printf("constraint %s (>) fails\n", cst.label.c_str());
            return false;
        }
    }

    for (auto& cst : pb.greater_equal_constraints) {
        if (compute_function(cst.elements, variable_value) < cst.value) {
            printf("constraint %s (>=) fails\n", cst.label.c_str());
            return false;
        }
    }

    for (auto& cst : pb.less_constraints) {
        if (compute_function(cst.elements, variable_value) >= cst.value) {
            printf("constraint %s (<) fails\n", cst.label.c_str());
            return false;
        }
    }

    for (auto& cst : pb.greater_constraints) {
        if (compute_function(cst.elements, variable_value) > cst.value) {
            printf("constraint %s (<=) fails\n", cst.label.c_str());
            return false;
        }
    }

    return true;
}

double
compute_solution(const problem& pb, const std::vector<int>& variable_value)
{
    Expects(not variable_value.empty(), "variables vector empty");

    double ret = pb.objective.constant;

    for (auto& elem : pb.objective.elements)
        ret += elem.factor * variable_value[elem.variable_index];

    return ret;
}
}
