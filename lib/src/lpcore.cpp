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
#include <iterator>
#include <string>
#include <utility>

#include <fmt/ostream.h>
#include <fmt/printf.h>

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstring>

namespace baryonyx {

void
context_deleter(context* ctx)
{
    delete ctx;
}

context_ptr
make_context(FILE* f, int verbose_level)
{
    auto pointer = new context(f, verbose_level);

    return context_ptr(pointer, &context_deleter);
}

context_ptr
make_context(string_logger_functor logger, int verbose_level)
{
    auto pointer = new context(logger, verbose_level);

    return context_ptr(pointer, &context_deleter);
}

context_ptr
copy_context(const context_ptr& ctx, FILE* f)
{
    auto pointer = new context(*ctx);

    pointer->string_logger = nullptr;
    pointer->cfile_logger = f ? f : stdout;
    pointer->logger = context::logger_type::c_file;

    return context_ptr(pointer, &context_deleter);
}

context_ptr
copy_context(const context_ptr& ctx, string_logger_functor logger)
{
    auto pointer = new context(*ctx);

    pointer->string_logger = logger;
    pointer->cfile_logger = nullptr;
    pointer->logger = context::logger_type::string;

    return context_ptr(pointer, &context_deleter);
}

context_ptr
copy_context(const context_ptr& ctx, FILE* f, int verbose_level)
{
    auto pointer = new context(*ctx);

    pointer->string_logger = nullptr;
    pointer->cfile_logger = f ? f : stdout;
    pointer->logger = context::logger_type::c_file;
    pointer->log_priority = static_cast<context::message_type>(
      verbose_level < 0 ? 0 : verbose_level > 7 ? 7 : verbose_level);

    return context_ptr(pointer, &context_deleter);
}

context_ptr
copy_context(const context_ptr& ctx,
             string_logger_functor logger,
             int verbose_level)
{
    auto pointer = new context(*ctx);

    pointer->string_logger = logger;
    pointer->cfile_logger = nullptr;
    pointer->logger = context::logger_type::string;
    pointer->log_priority = static_cast<context::message_type>(
      verbose_level < 0 ? 0 : verbose_level > 7 ? 7 : verbose_level);

    return context_ptr(pointer, &context_deleter);
}

problem
make_problem(const baryonyx::context_ptr& ctx, const std::string& filename)
{
    info(ctx, "problem reads from file `{}'\n", filename);

    std::ifstream ifs;
    ifs.exceptions(std::ifstream::badbit);
    ifs.open(filename);

    return read_problem(ifs);
}

problem
make_problem(const baryonyx::context_ptr& ctx, std::istream& is)
{
    info(ctx, "problem reads from stream\n");

    is.exceptions(std::ifstream::badbit);

    return read_problem(is);
}

result
make_result(const baryonyx::context_ptr& ctx, const std::string& filename)
{
    info(ctx, "solution reads from file {}\n", filename);

    std::ifstream ifs;
    ifs.exceptions(std::ifstream::badbit);
    ifs.open(filename);

    return read_result(ifs);
}

result
make_result(const baryonyx::context_ptr& ctx, std::istream& is)
{
    info(ctx, "solution reads from stream\n");

    is.exceptions(std::ifstream::badbit);

    return read_result(is);
}

std::ostream&
operator<<(std::ostream& os, const problem& p)
{
    if (not write_problem(os, p))
        os.setstate(std::ios_base::failbit);

    return os;
}

result
solve(const baryonyx::context_ptr& ctx, const problem& pb)
{
    check_consistency(pb);

    return solver_select(ctx, pb);
}

result
optimize(const baryonyx::context_ptr& ctx, const problem& pb)
{
    check_consistency(pb);

    return optimizer_select(ctx, pb);
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
is_valid_solution(const problem& pb, const std::vector<bool>& variable_value)
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
compute_solution(const problem& pb, const std::vector<bool>& variable_value)
{
    Expects(not variable_value.empty(), "variables vector empty");

    double ret = pb.objective.value;

    for (auto& elem : pb.objective.elements)
        ret += elem.factor * variable_value[elem.variable_index];

    return ret;
}

static std::vector<bool>
make_variable_value(const problem& pb, const result& r)
{
    if (not r or r.solutions.empty())
        return {};

    std::unordered_map<std::string, int> cache;

    std::transform(
      r.variable_name.cbegin(),
      r.variable_name.cend(),
      r.solutions.back().variables.cbegin(),
      std::inserter(cache, cache.begin()),
      [](const auto& name, int value) { return std::make_pair(name, value); });

    std::vector<bool> ret(pb.vars.names.size(), false);

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
    if (not r or r.solutions.empty())
        return false;

    Expects(pb.vars.names.size() == pb.vars.values.size());
    Expects(pb.vars.names.size() == r.variable_name.size());
    Expects(r.solutions.back().variables.size() == r.variable_name.size());
    Expects(pb.affected_vars.names.empty());
    Expects(pb.affected_vars.values.empty());

    auto variable_value = make_variable_value(pb, r);

    return is_valid_solution(pb, variable_value);
}

double
compute_solution(const problem& pb, const result& r)
{
    Expects(r and not r.solutions.empty());
    Expects(pb.vars.names.size() == pb.vars.values.size());
    Expects(pb.vars.names.size() == r.variable_name.size());
    Expects(r.solutions.back().variables.size() == r.variable_name.size());
    Expects(pb.affected_vars.names.empty());
    Expects(pb.affected_vars.values.empty());

    auto variable_value = make_variable_value(pb, r);

    return compute_solution(pb, variable_value);
}

} // namespace baryonyx
