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

#include <baryonyx/core-compare>
#include <baryonyx/core>

#include <fstream>
#include <iterator>
#include <thread>

#include "generalized-wedelin.hpp"
#include "itm.hpp"
#include "private.hpp"
#include "utils.hpp"
#include "wedelin.hpp"

namespace bx = baryonyx;

//
// Get number of thread to use in optimizer from parameters list or from the
// standard thread API. If an error occurred, this function returns 1.
//
static inline int
get_thread_number(std::shared_ptr<bx::context>& ctx) noexcept
{
    auto t = ctx->get_integer_parameter("thread",
                                        std::thread::hardware_concurrency());

    return t <= 0 ? 1 : bx::numeric_cast<int>(t);
}

static std::tuple<double, double, double, int>
get_parameters(std::shared_ptr<bx::context>& ctx) noexcept
{
    auto kappa = ctx->get_real_parameter("kappa", 0.001);
    auto theta = ctx->get_real_parameter("theta", 0.001);
    auto delta = ctx->get_real_parameter("delta", 0.001);
    auto limit = ctx->get_real_parameter("limit", 1000l);

    return std::make_tuple(kappa, delta, theta, limit);
}

template<typename variableT>
static bool
is_boolean_variable(const variableT& vars)
{
    for (const auto& elem : vars)
        if (elem.type != bx::variable_type::binary)
            return false;

    return true;
}

template<typename constraintsT>
static bool
is_boolean_coefficient(const constraintsT csts)
{
    for (const auto& cst : csts)
        for (const auto& elem : cst.elements)
            if (elem.factor < 0 or elem.factor > 1)
                return false;

    return true;
}

template<typename variableT>
static bool
is_integer_variable(const variableT& vars)
{
    for (const auto& elem : vars)
        if (elem.type != bx::variable_type::general)
            return false;

    return true;
}

template<typename constraintsT>
static bool
is_101_coefficient(const constraintsT csts)
{
    for (const auto& cst : csts)
        for (const auto& elem : cst.elements)
            if (elem.factor < -1 or elem.factor > 1)
                return false;

    return true;
}

namespace baryonyx_private {

bx::result
solve(std::shared_ptr<bx::context> ctx, bx::problem& pb)
{
    baryonyx_private::preprocess(ctx, pb);

    if (pb.greater_constraints.empty() and pb.less_constraints.empty() and
        is_boolean_coefficient(pb.equal_constraints) and
        is_boolean_variable(pb.vars.values))
        return baryonyx::itm::inequalities_1coeff_wedelin_solve(ctx, pb);

    if (is_101_coefficient(pb.equal_constraints) and
        is_101_coefficient(pb.less_constraints) and
        is_101_coefficient(pb.greater_constraints) and
        is_boolean_variable(pb.vars.values))
        return baryonyx::itm::inequalities_1coeff_wedelin_solve(ctx, pb);

    if ((is_101_coefficient(pb.equal_constraints) or
         is_101_coefficient(pb.greater_constraints) or
         is_101_coefficient(pb.less_constraints)) and
        is_integer_variable(pb.vars.values)) {
        double kappa, delta, theta;
        int limit;

        std::tie(kappa, delta, theta, limit) = get_parameters(ctx);

        return bx::generalized_wedelin(kappa, delta, theta, limit, pb);
    }

    ctx->info("no_solver_available");

    throw bx::solver_failure(bx::solver_error_tag::no_solver_available);
}

bx::result
optimize(std::shared_ptr<bx::context> ctx, bx::problem& pb)
{
    auto thread = get_thread_number(ctx);
    baryonyx_private::preprocess(ctx, pb);

    if (pb.greater_constraints.empty() and pb.less_constraints.empty() and
        is_boolean_coefficient(pb.equal_constraints) and
        is_boolean_variable(pb.vars.values))
        return bx::itm::inequalities_1coeff_wedelin_optimize(ctx, pb, thread);

    if (is_101_coefficient(pb.equal_constraints) and
        is_101_coefficient(pb.less_constraints) and
        is_101_coefficient(pb.greater_constraints) and
        is_boolean_variable(pb.vars.values))
        return bx::itm::inequalities_1coeff_wedelin_optimize(ctx, pb, thread);

    if ((is_101_coefficient(pb.equal_constraints) or
         is_101_coefficient(pb.greater_constraints) or
         is_101_coefficient(pb.less_constraints)) and
        is_integer_variable(pb.vars.values)) {
        double kappa, delta, theta;
        int limit;

        std::tie(kappa, delta, theta, limit) = get_parameters(ctx);

        return bx::generalized_wedelin(kappa, delta, theta, limit, pb);
    }

    ctx->info("no_solver_available");

    throw bx::solver_failure(bx::solver_error_tag::no_solver_available);
}
}
