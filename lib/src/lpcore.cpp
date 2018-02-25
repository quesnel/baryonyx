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

#include <baryonyx/core>

#include "private.hpp"
#include "utils.hpp"

#include <algorithm>
#include <fstream>
#include <utility>

#include <fmt/ostream.h>
#include <fmt/printf.h>

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstring>

#ifndef _WIN32
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

    return static_cast<int>(value);
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
    baryonyx::log(ctx,
                  baryonyx::context::message_type::info,
                  "Baryonyx v{}.{}.{}",
                  VERSION_MAJOR,
                  VERSION_MINOR,
                  VERSION_PATCH);

    if (VERSION_TWEAK)
        baryonyx::log(
          ctx, baryonyx::context::message_type::info, "-{}", VERSION_TWEAK);

    baryonyx::log(
      ctx,
      baryonyx::context::message_type::info,
      "\nGeneral options:\n"
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
      " * Global parameters"
      "  - limit: integer ]-oo, +oo[ in loop number\n"
      "  - time-limit: real [0, +oo[ in seconds\n"
      "  - floating-point-type: float double longdouble\n"
      "  - print-level: [0, 2]\n"
      " * In The Middle parameters\n"
      "  - preprocessing: none variables-number variables-weight "
      "constraints-weight implied\n"
      "  - constraint-order: none reversing random-sorting "
      "infeasibility-decr infeasibility-incr\n"
      "  - theta: real [0, 1]\n"
      "  - delta: real [0, +oo[\n"
      "  - kappa-min: real [0, 1[\n"
      "  - kappa-step: real [0, 1[\n"
      "  - kappa-max: real [0, 1[\n"
      "  - alpha: integer [0, 2]\n"
      "  - w: integer [0, +oo[\n"
      "  - norm: l1 l2 inf none rng\n"
      " * Pushes system parameters\n"
      "  - pushes-limit: integer [0, +oo[\n"
      "  - pushing-objective-amplifier: real [0, +oo[\n"
      "  - pushing-iteration-limit: integer [0, +oo[\n"
      "  - pushing-k-factor: real [0, +oo[\n"
      " * Initialization parameters\n"
      "  - init-policy: bastert random best\n"
      "  - init-random: real [0, 1]\n");
}
}

namespace baryonyx {

context::context()
  : m_cfile_logger(stdout)
  , m_log_priority(context::message_type::info)
  , m_logger(context::logger_type::c_file)
{}

context::context(FILE* f)
  : m_cfile_logger(f ? f : stdout)
  , m_log_priority(context::message_type::info)
  , m_logger(context::logger_type::c_file)
{}

context::context(string_logger_functor logger)
  : m_string_logger(logger)
  , m_cfile_logger(nullptr)
  , m_log_priority(context::message_type::info)
  , m_logger(context::logger_type::string)
{}

int
context::parse(int argc, char* argv[]) noexcept
{
    const char* const short_opts = "OC:hp:l:qv:";
    const struct option long_opts[] = {
        { "optimize", 0, nullptr, 'O' }, { "check", 1, nullptr, 'C' },
        { "help", 0, nullptr, 'h' },     { "param", 1, nullptr, 'p' },
        { "limit", 1, nullptr, 'l' },    { "quiet", 0, nullptr, 'q' },
        { "verbose", 1, nullptr, 'v' },  { nullptr, 0, nullptr, 0 }
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
            baryonyx::log(this,
                          baryonyx::context::message_type::err,
                          "Unknown command line option\n");
            return -1;
        };
    }

    if (quiet)
        set_log_priority(message_type::notice);
    else if (verbose >= 0 and verbose <= 7)
        set_log_priority(static_cast<message_type>(verbose));

    return ::optind;
}

void
context::set_parameters(std::unordered_map<std::string, parameter>&& params)
{
    m_parameters = params;
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

    m_parameters[name] = p;
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

    baryonyx::log(const_cast<baryonyx::context*>(this),
                  context::message_type::warning,
                  "fail to convert parameter {}\n",
                  name);

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

    baryonyx::log(const_cast<baryonyx::context*>(this),
                  context::message_type::warning,
                  "fail to convert parameter {}\n",
                  name);

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

    baryonyx::log(const_cast<baryonyx::context*>(this),
                  context::message_type::warning,
                  "fail to convert parameter {}\n",
                  name);

    return def;
}

problem
make_problem(const std::shared_ptr<baryonyx::context>& ctx,
             const std::string& filename)
{
    info(ctx, "problem reads from file `{}'\n", filename);

    std::ifstream ifs;
    ifs.exceptions(std::ifstream::badbit);
    ifs.open(filename);

    return baryonyx_private::read_problem(ifs);
}

problem
make_problem(const std::shared_ptr<baryonyx::context>& ctx, std::istream& is)
{
    info(ctx, "problem reads from stream\n");

    is.exceptions(std::ifstream::badbit);

    return baryonyx_private::read_problem(is);
}

result
make_result(const std::shared_ptr<baryonyx::context>& ctx,
            const std::string& filename)
{
    info(ctx, "solution reads from file {}\n", filename);

    std::ifstream ifs;
    ifs.exceptions(std::ifstream::badbit);
    ifs.open(filename);

    return baryonyx_private::read_result(ifs);
}

result
make_result(const std::shared_ptr<baryonyx::context>& ctx, std::istream& is)
{
    info(ctx, "solution reads from stream\n");

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

    return baryonyx_private::optimize(std::move(ctx), pb);
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
