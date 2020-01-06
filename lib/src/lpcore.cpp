/* Copyright (C) 2016-2019 INRA
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

#include <baryonyx/core-utils>
#include <baryonyx/core>

#include "debug.hpp"
#include "itm.hpp"
#include "private.hpp"
#include "problem.hpp"
#include "utils.hpp"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <string>
#include <unordered_map>
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
make_context(int verbose_level)
{
    return context_ptr(new context(verbose_level), &context_deleter);
}

void
context_register(const context_ptr& ctx,
                 solver_started_cb start,
                 solver_updated_cb update,
                 solver_finished_cb finish)
{
    ctx->start = std::move(start);
    ctx->update = std::move(update);
    ctx->finish = std::move(finish);
}

raw_problem
make_problem(const baryonyx::context_ptr& ctx, const std::string& filename)
{
    info(ctx, "problem reads from file `{}'\n", filename);

    std::ifstream ifs(filename);
    if (!ifs.is_open())
        return raw_problem(baryonyx::file_format_error_tag::file_not_found);

    return make_problem(ctx, ifs);
}

result
solve(const baryonyx::context_ptr& ctx, const raw_problem& rawpb)
{
    return (ctx->parameters.preprocessor ==
            solver_parameters::preprocessor_options::all)
             ? itm::solve(ctx, preprocess(ctx, rawpb))
             : itm::solve(ctx, unpreprocess(ctx, rawpb));
}

result
optimize(const baryonyx::context_ptr& ctx, const raw_problem& rawpb)
{
    if ((ctx->parameters.mode & solver_parameters::mode_type::branch) ==
        solver_parameters::mode_type::branch)
        return (ctx->parameters.preprocessor ==
                solver_parameters::preprocessor_options::all)
                 ? itm::branch_optimize(ctx, preprocess(ctx, rawpb))
                 : itm::branch_optimize(ctx, unpreprocess(ctx, rawpb));

    if ((ctx->parameters.mode & solver_parameters::mode_type::nlopt) ==
        solver_parameters::mode_type::nlopt)
        return (ctx->parameters.preprocessor ==
                solver_parameters::preprocessor_options::all)
                 ? itm::nlopt_optimize(ctx, preprocess(ctx, rawpb))
                 : itm::nlopt_optimize(ctx, unpreprocess(ctx, rawpb));

    if ((ctx->parameters.mode & solver_parameters::mode_type::manual) ==
        solver_parameters::mode_type::manual)
        return (ctx->parameters.preprocessor ==
                solver_parameters::preprocessor_options::all)
                 ? itm::manual_optimize(ctx, preprocess(ctx, rawpb))
                 : itm::manual_optimize(ctx, unpreprocess(ctx, rawpb));

    return (ctx->parameters.preprocessor ==
            solver_parameters::preprocessor_options::all)
             ? itm::optimize(ctx, preprocess(ctx, rawpb))
             : itm::optimize(ctx, unpreprocess(ctx, rawpb));
}

template<typename functionT, typename variablesT>
static int
compute_function(const functionT& fct, const variablesT& vars) noexcept
{
    int v{ 0 };

    for (auto& f : fct)
        v += vars[f.variable_index] * f.factor;

    /* Check if this is more efficient to remove the condition branch and
       replace with the 0-1 mult f.factor.

       if (vars[f.variable_index])
        v += f.factor;

     */

    return v;
}

static bool
is_valid_solution_impl(const raw_problem& pb,
                       const std::vector<var_value>& variable_value)
{
    bx_expects(!variable_value.empty());
    bx_expects(variable_value.size() == pb.vars.names.size());
    bx_expects(variable_value.size() == pb.vars.values.size());

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

static double
compute_solution_impl(const raw_problem& pb,
                      const std::vector<var_value>& variable_value)
{
    bx_expects(!variable_value.empty());

    double ret = pb.objective.value;

    for (auto& elem : pb.objective.elements)
        ret += elem.factor * variable_value[elem.variable_index];

    for (auto& elem : pb.objective.qelements)
        ret += elem.factor * variable_value[elem.variable_index_a] *
               variable_value[elem.variable_index_b];

    return ret;
}

template<typename Problem>
static std::vector<var_value>
make_variable_value(const Problem& pb, const result& r)
{
    if (!r || r.solutions.empty())
        return {};

    std::unordered_map<std::string_view, bool> cache;

    bx_ensures(r.affected_vars.names.size() == r.affected_vars.values.size());

    for (size_t i = 0, e = r.affected_vars.names.size(); i != e; ++i)
        cache[r.affected_vars.names[i]] = r.affected_vars.values[i];

    bx_ensures(r.variable_name.size() == r.solutions.back().variables.size());

    for (size_t i = 0, e = r.variable_name.size(); i != e; ++i)
        cache[r.variable_name[i]] = r.solutions.back().variables[i];

    std::vector<var_value> ret(pb.vars.names.size(), false);

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
    bx_expects(pb.vars.names.size() ==
               r.affected_vars.names.size() + r.variable_name.size());
    bx_expects(r.solutions.back().variables.size() == r.variable_name.size());

    return compute_solution_impl(pb, make_variable_value(pb, r));
}

} // namespace baryonyx
