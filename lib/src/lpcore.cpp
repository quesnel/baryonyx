/* Copyright (C) 2017 INRA
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

#include <baryonyx/core>

#include "lpformat-consistency.hpp"
#include "lpformat-io.hpp"
#include "mitm.hpp"
#include "utils.hpp"

#include <algorithm>
#include <fstream>

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>

#ifdef __unix__
#include <getopt.h>
#include <sys/types.h>
#include <unistd.h>
#endif

namespace {

double
to_double(const char* s, double bad_value) noexcept
{
    char* c;
    errno = 0;
    double value = std::strtod(s, &c);

    if ((errno == ERANGE and (value == HUGE_VAL or value == -HUGE_VAL)) or
        (value == 0.0 and c == s))
        return bad_value;

    return value;
}

long
to_long(const char* s, long bad_value) noexcept
{
    char* c;
    errno = 0;
    long value = std::strtol(s, &c, 10);

    if ((errno == ERANGE and (value == LONG_MIN or value == LONG_MAX)) or
        (value == 0 and c == s))
        return bad_value;

    return value;
}

std::tuple<std::string, baryonyx::parameter>
split_param(const char* param) noexcept
{
    std::string name, value;

    while (*param) {
        if (isalpha(*param) or *param == '_' or *param == '-')
            name += *param;
        else
            break;

        param++;
    }

    if (*param and (*param == ':' or *param == '=')) {
        param++;

        while (*param)
            value += *param++;
    }

    auto valuel = to_long(value.c_str(), LONG_MIN);
    auto valued = to_double(value.c_str(), -HUGE_VAL);

    double tmp;
    if (valued != -HUGE_VAL and std::modf(valued, &tmp))
        return std::make_tuple(name, baryonyx::parameter(valued));

    if (valuel != LONG_MIN)
        return std::make_tuple(name, baryonyx::parameter(valuel));

    return std::make_tuple(name, baryonyx::parameter(value));
}

void
help(std::shared_ptr<baryonyx::context> ctx) noexcept
{
    ctx->info("--help|-h                   This help message\n"
              "--param|-p [name]:[value]   Add a new parameter (name is"
              " [a-z][A-Z]_ value can be a double, an integer otherwise a"
              " string.\n"
              "--optimize|-O               Optimize model (default "
              "feasibility search only)\n"
              "--limit int                 Set limit\n"
              "--quiet                     Remove any verbose message\n"
              "--verbose|-v int            Set verbose level\n");
}
}

namespace baryonyx {

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

    void write(baryonyx::context::message_type m,
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

int
context::parse(int argc, char* argv[]) noexcept
{
    const char* const short_opts = "Ohp:l:qv:";
    const struct option long_opts[] = { { "optimize", 0, nullptr, 'O' },
                                        { "help", 0, nullptr, 'h' },
                                        { "param", 1, nullptr, 'p' },
                                        { "limit", 1, nullptr, 'l' },
                                        { "quiet", 0, nullptr, 'q' },
                                        { "verbose", 1, nullptr, 'v' },
                                        { 0, 0, nullptr, 0 } };

    int opt_index;
    int verbose = 1;
    int quiet = 0;

    auto ctx = std::make_shared<baryonyx::context>();
    ctx->set_standard_stream_logger();

    for (;;) {
        const auto opt =
          getopt_long(argc, argv, short_opts, long_opts, &opt_index);
        if (opt == -1)
            break;

        switch (opt) {
        case 0:
            break;
        case 'O':
            m_optimize = true;
            break;
        case 'l':
            m_parameters["limit"] = ::to_long(::optarg, 1000l);
            break;
        case 'h':
            ::help(ctx);
            break;
        case 'p': {
            std::string name;
            baryonyx::parameter value;
            std::tie(name, value) = ::split_param(::optarg);
            m_parameters[name] = value;
        } break;
        case 'q':
            quiet = 1;
            break;
        case '?':
        default:
            ctx->error("Unknown command line option\n");
            return -1;
        };
    }

    if (quiet)                    // priority to quiet over the
        ctx->set_log_priority(3); // verbose mode.
    else if (verbose >= 0 and verbose <= 7)
        ctx->set_log_priority(verbose);

    return ::optind;
}

void
context::set_parameter(const std::string& name, double p) noexcept
{
    if (name.empty())
        return;

    if (not std::isalnum(name[0]))
        return;

    m_parameters[name] = p;
}

void
context::set_parameter(const std::string& name, long int p) noexcept
{
    if (name.empty())
        return;

    if (not std::isalnum(name[0]))
        return;

    m_parameters[name] = p;
}

void
context::set_parameter(const std::string& name, std::string p) noexcept
{
    if (name.empty())
        return;

    if (not std::isalnum(name[0]))
        return;

    m_parameters[name] = std::move(p);
}

double
context::get_real_parameter(const std::string& name, double def) const noexcept
{
    auto it = m_parameters.find(name);
    if (it == m_parameters.cend())
        return def;

    if (it->second.type == parameter::tag::real)
        return it->second.d;

    if (it->second.type == parameter::tag::integer)
        return static_cast<double>(it->second.l);

    warning("fail to convert parameter %s\n", name.c_str());

    return def;
}

long int
context::get_integer_parameter(const std::string& name, long int def) const
  noexcept
{
    auto it = m_parameters.find(name);
    if (it == m_parameters.cend())
        return def;

    if (it->second.type == parameter::tag::integer)
        return it->second.l;

    if (it->second.type == parameter::tag::real)
        return static_cast<long>(it->second.d);

    warning("fail to convert parameter %s\n", name.c_str());

    return def;
}

std::string
context::get_string_parameter(const std::string& name, std::string def) const
  noexcept
{
    auto it = m_parameters.find(name);
    if (it == m_parameters.cend())
        return def;

    if (it->second.type == parameter::tag::string)
        return it->second.s;

    if (it->second.type == parameter::tag::real)
        return std::to_string(it->second.d);

    if (it->second.type == parameter::tag::integer)
        return std::to_string(it->second.l);

    warning("fail to convert parameter %s\n", name.c_str());

    return def;
}

void
context::set_log_priority(int priority) noexcept
{
    m_log_priority = priority < 0 ? 0 : 7 < priority ? 7 : priority;
}

int
context::get_log_priority() const noexcept
{
    return m_log_priority;
}

bool
context::optimize() const noexcept
{
    return m_optimize;
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
//
// Default, the logging system is active and the call to the @c log function
// are send to the logger functor. Define LP_DISABLE_LOGGING as preprocessor
// value to hide all logging message..
//
void
context::log(message_type type, const char* format, ...) const noexcept
{
    if (not m_logger)
        return;

    switch (type) {
    case baryonyx::context::message_type::emerg:
        break;
    case baryonyx::context::message_type::alert:
        if (m_log_priority < 1)
            return;
        break;
    case baryonyx::context::message_type::crit:
        if (m_log_priority < 2)
            return;
        break;
    case baryonyx::context::message_type::err:
        if (m_log_priority < 3)
            return;
        break;
    case baryonyx::context::message_type::warning:
        if (m_log_priority < 4)
            return;
        break;
    case baryonyx::context::message_type::notice:
        if (m_log_priority < 5)
            return;
        break;
    case baryonyx::context::message_type::info:
        if (m_log_priority < 6)
            return;
        break;
    case baryonyx::context::message_type::debug:
        if (m_log_priority < 7)
            return;
        break;
    }

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
             ...) const noexcept
{
    if (not m_logger and m_log_priority < priority)
        return;

    va_list args;

    va_start(args, format);
    m_logger->write(priority, file, line, fn, format, args);
    va_end(args);
}

void
context::info(const char* format, ...) const noexcept
{
    if (not m_logger or m_log_priority < 6)
        return;

    va_list args;

    va_start(args, format);
    m_logger->write(context::message_type::info, format, args);
    va_end(args);
}

void
context::debug(const char* format, ...) const noexcept
{
    if (not m_logger or m_log_priority < 7)
        return;

    va_list args;

    va_start(args, format);
    m_logger->write(context::message_type::debug, format, args);
    va_end(args);
}

void
context::warning(const char* format, ...) const noexcept
{
    if (not m_logger and m_log_priority < 4)
        return;

    va_list args;

    va_start(args, format);
    m_logger->write(context::message_type::warning, format, args);
    va_end(args);
}

void
context::error(const char* format, ...) const noexcept
{
    if (not m_logger and m_log_priority < 3)
        return;

    va_list args;

    va_start(args, format);
    m_logger->write(context::message_type::err, format, args);
    va_end(args);
}
#else
void
context::log(message_type, const char*, va_list) const noexcept
{
}

void
context::log(int, const char*, int, const char*, const char*, va_list) const
  noexcept
{
}

void
context::info(const char*, va_list) const noexcept
{
}

void
context::warning(const char*, va_list) const noexcept
{
}

void
context::debug(const char*, va_list) const noexcept
{
}

void
context::error(const char*, va_list) const noexcept
{
}
#endif

problem
make_problem(std::shared_ptr<baryonyx::context> ctx, const std::string& filename)
{
    ctx->info("problem read from file `%s'\n", filename.c_str());

    std::ifstream ifs;
    ifs.exceptions(std::ifstream::badbit);
    ifs.open(filename);

    return details::read_problem(ifs);
}

problem
make_problem(std::shared_ptr<baryonyx::context> ctx, std::istream& is)
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
solve(std::shared_ptr<baryonyx::context> ctx, problem& pb)
{
    check(pb);

    return mitm_solve(ctx, pb);
}

result
optimize(std::shared_ptr<baryonyx::context> ctx, problem& pb)
{
    check(pb);

    return mitm_optimize(ctx, pb);
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

} // namespace baryonyx
