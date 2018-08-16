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

#include "debug.hpp"
#include "private.hpp"
#include "problem.hpp"
#include "utils.hpp"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <string>
#include <utility>

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

raw_problem
make_problem(const baryonyx::context_ptr& ctx, const std::string& filename)
{
    info(ctx, "problem reads from file `{}'\n", filename);

    raw_problem pb;

    std::ifstream ifs;
    ifs.open(filename);

    ifs >> pb;

    return pb;
}

result
solve(const baryonyx::context_ptr& ctx,
      const raw_problem& rawpb,
      preprocessor_options pp_opt)
{
    return (pp_opt == preprocessor_options::all)
             ? solver_select(ctx, preprocess(ctx, rawpb))
             : solver_select(ctx, unpreprocess(ctx, rawpb));
}

result
optimize(const baryonyx::context_ptr& ctx,
         const raw_problem& rawpb,
         preprocessor_options pp_opt)
{
    return (pp_opt == preprocessor_options::all)
             ? optimizer_select(ctx, preprocess(ctx, rawpb))
             : optimizer_select(ctx, unpreprocess(ctx, rawpb));
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
is_valid_solution_impl(const raw_problem& pb,
                       const std::vector<bool>& variable_value)
{
    bx_expects(!variable_value.empty());
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
compute_solution_impl(const raw_problem& pb,
                      const std::vector<bool>& variable_value)
{
    bx_expects(!variable_value.empty());

    double ret = pb.objective.value;

    for (auto& elem : pb.objective.elements)
        ret += elem.factor * variable_value[elem.variable_index];

    return ret;
}

template<typename Problem>
static std::vector<bool>
make_variable_value(const Problem& pb, const result& r)
{
    if (!r || r.solutions.empty())
        return {};

    std::unordered_map<std::string, bool> cache;

    bx_ensures(r.affected_vars.names.size() == r.affected_vars.values.size());

    std::transform(r.affected_vars.names.cbegin(),
                   r.affected_vars.names.cend(),
                   r.affected_vars.values.cbegin(),
                   std::inserter(cache, cache.begin()),
                   [](const auto& name, bool value) {
                       return std::make_pair(name, value);
                   });

    bx_ensures(r.variable_name.size() == r.solutions.back().variables.size());

    std::transform(r.variable_name.cbegin(),
                   r.variable_name.cend(),
                   r.solutions.back().variables.cbegin(),
                   std::inserter(cache, cache.begin()),
                   [](const auto& name, bool value) {
                       return std::make_pair(name, value);
                   });

    std::vector<bool> ret(pb.vars.names.size(), false);

    for (std::size_t i = 0, e = pb.vars.names.size(); i != e; ++i) {
        auto it = cache.find(pb.vars.names[i]);
        bx_expects(it != cache.end());
        ret[i] = it->second;
    }

    return ret;
}

bool
is_valid_solution(const raw_problem& pb, const result& r)
{
    if (!r || r.solutions.empty())
        return false;

    bx_expects(pb.vars.names.size() == pb.vars.values.size());
    bx_expects(pb.vars.names.size() ==
               r.variable_name.size() + r.affected_vars.names.size());
    bx_expects(r.solutions.back().variables.size() == r.variable_name.size());

    return is_valid_solution_impl(pb, make_variable_value(pb, r));
}

double
compute_solution(const raw_problem& pb, const result& r)
{
    bx_expects(r && !r.solutions.empty());
    bx_expects(pb.vars.names.size() == pb.vars.values.size());
    bx_expects(pb.vars.names.size() == r.variable_name.size());
    bx_expects(r.solutions.back().variables.size() == r.variable_name.size());

    return compute_solution_impl(pb, make_variable_value(pb, r));
}

} // namespace baryonyx
