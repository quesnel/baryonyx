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

namespace baryonyx {

context::context()
  : m_cfile_logger(stdout)
  , m_log_priority(context::message_type::info)
  , m_logger(context::logger_type::c_file)
  , m_optimize(false)
  , m_check(false)
{
}

context::context(FILE* f)
  : m_cfile_logger(f ? f : stdout)
  , m_log_priority(context::message_type::info)
  , m_logger(context::logger_type::c_file)
  , m_optimize(false)
  , m_check(false)
{
}

context::context(string_logger_functor logger)
  : m_string_logger(logger)
  , m_cfile_logger(nullptr)
  , m_log_priority(context::message_type::info)
  , m_logger(context::logger_type::string)
  , m_optimize(false)
  , m_check(false)
{
}

void
context::set_parameters(std::unordered_map<std::string, parameter>&& params)
{
    m_parameters = params;

    std::unordered_map<std::string, baryonyx::parameter>::const_iterator it;

    it = m_parameters.find("check-filename");
    if (it != m_parameters.cend())
        m_check = true;

    it = m_parameters.find("optimize");
    if (it != m_parameters.cend())
        m_optimize = true;

    int quiet = 0;
    it = m_parameters.find("quiet");
    if (it != m_parameters.cend())
        quiet = 1;

    int verbose = 6;
    it = m_parameters.find("verbose");
    if (it != m_parameters.cend() and
        it->second.type == parameter::tag::integer)
        verbose = baryonyx::clamp(it->second.l, 0, 7);

    if (quiet)
        set_log_priority(message_type::notice);
    else if (verbose >= 0 and verbose <= 7)
        set_log_priority(static_cast<message_type>(verbose));
}

bool
context::try_set_important_parameter(std::string name,
                                     baryonyx::parameter value)
{
    if (name == "check-filename" and value.type == parameter::tag::string) {
        m_check = true;
        m_parameters["check-filename"] = value;
        return true;
    }

    if (name == "optimize") {
        m_optimize = true;
        return true;
    }

    if (name == "quiet") {
        set_log_priority(message_type::notice);
        return true;
    }

    if (name == "verbose" and value.type == parameter::tag::integer) {
        set_log_priority(
          static_cast<message_type>(baryonyx::clamp(value.l, 0, 7)));
        return true;
    }

    return false;
}

void
context::set_parameter(const std::string& name, double p) noexcept
{
    if (name.empty())
        return;

    if (not std::isalnum(name[0]))
        return;

    if (not try_set_important_parameter(name, baryonyx::parameter(p)))
        m_parameters[name] = p;
}

void
context::set_parameter(const std::string& name, int p) noexcept
{
    if (name.empty())
        return;

    if (not std::isalnum(name[0]))
        return;

    if (not try_set_important_parameter(name, baryonyx::parameter(p)))
        m_parameters[name] = p;
}

void
context::set_parameter(const std::string& name, std::string p) noexcept
{
    if (name.empty())
        return;

    if (not std::isalnum(name[0]))
        return;

    if (not try_set_important_parameter(name, baryonyx::parameter(p)))
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
