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

#include "private.hpp"
#include "utils.hpp"

#include <algorithm>
#include <fstream>

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <getopt.h>

#ifndef _WIN32
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

int
to_int(const char* s, int bad_value) noexcept
{
    char* c;
    errno = 0;
    long value = std::strtol(s, &c, 10);

    if ((errno == ERANGE and (value == LONG_MIN or value == LONG_MAX)) or
        (value == 0 and c == s))
        return bad_value;

    if (value < INT_MIN)
        return INT_MIN;

    if (value > INT_MAX)
        return INT_MAX;

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

    auto valuel = to_int(value.c_str(), INT_MIN);
    auto valued = to_double(value.c_str(), -HUGE_VAL);

    double tmp;
    if (valued != -HUGE_VAL and std::modf(valued, &tmp))
        return std::make_tuple(name, baryonyx::parameter(valued));

    if (valuel != INT_MIN)
        return std::make_tuple(name, baryonyx::parameter(valuel));

    return std::make_tuple(name, baryonyx::parameter(value));
}

void
help(baryonyx::context* ctx) noexcept
{
    ctx->info(
      "Baryonyx v%d.%d.%d", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);

    if (VERSION_TWEAK)
        ctx->info("-%d", VERSION_TWEAK);

    ctx->info("\nGeneral options:\n"
              "  --help|-h                   This help message\n"
              "  --param|-p [name][:|=][value]   Add a new parameter (name is"
              " [a-z][A-Z]_) value can be a double, an integer otherwise a"
              " string.\n"
              "  --optimize|-O               Optimize model (default "
              "feasibility search only)\n"
              "  --check filename.sol        Check if the solution is correct."
              "\n"
              "  --quiet                     Remove any verbose message\n"
              "  --verbose|-v int            Set verbose level\n\n"
              "Parameter list for in the middle heuristic\n"
              "solver parameters:\n"
              "  - preprocessing: none variables-number variables-weight "
              "constraints-weight implied\n"
              "  - constraint-order: none reversing random-sorting "
              "infeasibility-decr ingeasibility-incr\n"
              "  - time-limit: real [0, +oo[ in second\n"
              "  - theta: real [0, 1]\n"
              "  - delta: real [0, +oo[\n"
              "  - limit: integer [0, +oo[\n"
              "  - kappa-min: real [0, 1[\n"
              "  - kappa-step: real [0, 1[\n"
              "  - kappa-max: real [0, 1[\n"
              "  - alpha: integer [0, 2]\n"
              "  - w: integer [0, +oo[\n"
              "  - norm: l1 l2 inf none rng\n"
              "  - print-level: [0, 2]\n"
              "  - floating-point-type: float, double, longdouble\n"
              "optimizer parameters:\n"
              "  - pushes-limit: integer [0, +oo[\n"
              "  - pushing-objective-amplifier: real [0, +oo[\n"
              "  - pushing-iteration-limit: integer [0, +oo[\n"
              "  - pushing-k-factor: real [0, +oo[\n");
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
#ifndef __WIN32
        if (::isatty(STDOUT_FILENO)) {
            switch (m) {
            case context::message_type::emerg:
            case context::message_type::alert:
            case context::message_type::crit:
            case context::message_type::err:
                ::printf("\033[30m\033[2m");
                break;
            case context::message_type::warning:
                ::printf("\033[32m\033[1m");
                break;
            case context::message_type::notice:
                break;
            case context::message_type::info:
                break;
            case context::message_type::debug:
                ::printf("\033[33m\033[1m");
                break;
            }

            vfprintf(stdout, format, args);

            switch (m) {
            case context::message_type::emerg:
            case context::message_type::alert:
            case context::message_type::crit:
            case context::message_type::err:
                ::printf("\033[30m\033[0m");
                break;
            case context::message_type::warning:
                ::printf("\033[30m\033[0m");
                break;
            case context::message_type::notice:
                break;
            case context::message_type::info:
                break;
            case context::message_type::debug:
                ::printf("\033[30m\033[0m");
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
    const char* const short_opts = "OC:hp:l:qv:";
    const struct option long_opts[] = {
        { "optimize", 0, nullptr, 'O' }, { "check", 1, nullptr, 'C' },
        { "help", 0, nullptr, 'h' },     { "param", 1, nullptr, 'p' },
        { "limit", 1, nullptr, 'l' },    { "quiet", 0, nullptr, 'q' },
        { "verbose", 1, nullptr, 'v' },  { 0, 0, nullptr, 0 }
    };

    int opt_index;
    int verbose = 6;
    int quiet = 0;

    for (;;) {
        const auto opt =
          getopt_long(argc, argv, short_opts, long_opts, &opt_index);
        if (opt == -1)
            break;

        switch (opt) {
        case 0:
            break;
        case 'C':
            m_check = true;
            m_parameters["check-filename"] = std::string(::optarg);
            break;
        case 'O':
            m_optimize = true;
            break;
        case 'h':
            ::help(this);
            break;
        case 'p': {
            std::string name;
            baryonyx::parameter value;
            std::tie(name, value) = ::split_param(::optarg);
            m_parameters[name] = value;
        } break;
        case 'l':
            m_parameters["limit"] = ::to_int(::optarg, 1000);
            break;
        case 'q':
            quiet = 1;
            break;
        case 'v':
            verbose = ::to_int(::optarg, 3);
            break;
        case '?':
        default:
            error("Unknown command line option\n");
            return -1;
        };
    }

    if (quiet)
        set_log_priority(4); // verbose mode.
    else if (verbose >= 0 and verbose <= 7)
        set_log_priority(verbose);

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
context::set_parameter(const std::string& name, int p) noexcept
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

int
context::get_integer_parameter(const std::string& name, int def) const noexcept
{
    auto it = m_parameters.find(name);
    if (it == m_parameters.cend())
        return def;

    if (it->second.type == parameter::tag::integer)
        return it->second.l;

    if (it->second.type == parameter::tag::real)
        return static_cast<int>(it->second.d);

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

bool
context::check() const noexcept
{
    return m_check;
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

#ifndef BARYONYX_DISABLE_LOGGING
//
// Default, the logging system is active and the call to the @c log function
// are send to the logger functor. Define BARYONYX_DISABLE_LOGGING as
// preprocessor value to hide all logging message..
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
make_problem(std::shared_ptr<baryonyx::context> ctx,
             const std::string& filename)
{
    ctx->info("problem reads from file `%s'\n", filename.c_str());

    std::ifstream ifs;
    ifs.exceptions(std::ifstream::badbit);
    ifs.open(filename);

    return baryonyx_private::read_problem(ifs);
}

problem
make_problem(std::shared_ptr<baryonyx::context> ctx, std::istream& is)
{
    ctx->info("problem reads from stream\n");

    is.exceptions(std::ifstream::badbit);

    return baryonyx_private::read_problem(is);
}

result
make_result(std::shared_ptr<baryonyx::context> ctx,
            const std::string& filename)
{
    ctx->info("solution reads from file `%s'\n", filename.c_str());

    std::ifstream ifs;
    ifs.exceptions(std::ifstream::badbit);
    ifs.open(filename);

    return baryonyx_private::read_result(ifs);
}

result
make_result(std::shared_ptr<baryonyx::context> ctx, std::istream& is)
{
    ctx->info("solution reads from stream\n");

    is.exceptions(std::ifstream::badbit);

    return baryonyx_private::read_result(is);
}

std::ostream&
operator<<(std::ostream& os, const problem& p)
{
    if (not baryonyx_private::write_problem(os, p))
        os.setstate(std::ios_base::failbit);

    return os;
}

result
solve(std::shared_ptr<baryonyx::context> ctx, problem& pb)
{
    baryonyx_private::check_consistency(pb);

    return baryonyx_private::solve(ctx, pb);
}

result
optimize(std::shared_ptr<baryonyx::context> ctx, problem& pb)
{
    baryonyx_private::check_consistency(pb);

    return baryonyx_private::optimize(ctx, pb);
}

template<typename functionT, typename variablesT>
static int
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
    std::size_t i, e;

    for (i = 0, e = pb.equal_constraints.size(); i != e; ++i) {
        if (compute_function(pb.equal_constraints[i].elements,
                             variable_value) != pb.equal_constraints[i].value)
            return false;
    }

    for (i = 0, e = pb.less_constraints.size(); i != e; ++i) {
        if (compute_function(pb.less_constraints[i].elements, variable_value) >
            pb.less_constraints[i].value)
            return false;
    }

    for (i = 0, e = pb.greater_constraints.size(); i != e; ++i) {
        if (compute_function(pb.greater_constraints[i].elements,
                             variable_value) < pb.greater_constraints[i].value)
            return false;
    }

    return true;
}

double
compute_solution(const problem& pb, const std::vector<int>& variable_value)
{
    Expects(not variable_value.empty(), "variables vector empty");

    double ret = pb.objective.value;

    for (auto& elem : pb.objective.elements)
        ret += elem.factor * variable_value[elem.variable_index];

    return ret;
}

static std::vector<int>
make_variable_value(const problem& pb, const result& r)
{
    std::unordered_map<std::string, int> cache;

    std::transform(
      r.variable_name.cbegin(),
      r.variable_name.cend(),
      r.variable_value.cbegin(),
      std::inserter(cache, cache.begin()),
      [](const auto& name, int value) { return std::make_pair(name, value); });

    std::vector<int> ret(pb.vars.names.size(), INT_MAX);

    for (std::size_t i = 0, e = pb.vars.names.size(); i != e; ++i) {
        auto it = cache.find(pb.vars.names[i]);
        Expects(it != cache.end());
        ret[i] = it->second;
    }

    return ret;
}

bool
is_valid_solution(const problem& pb, const result& r)
{
    Expects(pb.vars.names.size() == pb.vars.values.size());
    Expects(pb.vars.names.size() == r.variable_name.size());
    Expects(r.variable_value.size() == r.variable_name.size());
    Expects(pb.affected_vars.names.empty());
    Expects(pb.affected_vars.values.empty());

    auto variable_value = make_variable_value(pb, r);

    return is_valid_solution(pb, variable_value);
}

double
compute_solution(const problem& pb, const result& r)
{
    Expects(pb.vars.names.size() == pb.vars.values.size());
    Expects(pb.vars.names.size() == r.variable_name.size());
    Expects(r.variable_value.size() == r.variable_name.size());
    Expects(pb.affected_vars.names.empty());
    Expects(pb.affected_vars.values.empty());

    auto variable_value = make_variable_value(pb, r);

    return compute_solution(pb, variable_value);
}

} // namespace baryonyx
