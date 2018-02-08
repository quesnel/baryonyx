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

#include <baryonyx/core-compare>
#include <baryonyx/core>

#include <fstream>
#include <iterator>
#include <thread>

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

template<typename variableT>
static bool
is_boolean_variable(const variableT& vars) noexcept
{
    for (const auto& elem : vars)
        if (elem.type != bx::variable_type::binary)
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
is_boolean_coefficient(const constraintsT& csts) noexcept
{
    for (const auto& cst : csts)
        for (const auto& elem : cst.elements)
            if (elem.factor < 0 or elem.factor > 1)
                return false;

    return true;
}

template<typename constraintsT>
static bool
is_101_coefficient(const constraintsT& csts)
{
    for (const auto& cst : csts)
        for (const auto& elem : cst.elements)
            if (elem.factor < -1 or elem.factor > 1)
                return false;

    return true;
}

static bool
is_101_coefficient(const bx::problem& pb) noexcept
{
    return is_101_coefficient(pb.equal_constraints) and
           is_101_coefficient(pb.less_constraints) and
           is_101_coefficient(pb.greater_constraints);
}

namespace baryonyx_private {

bx::result
solve(const std::shared_ptr<bx::context>& ctx, bx::problem& pb)
{
    baryonyx_private::preprocess(ctx, pb);

    if (is_boolean_variable(pb.vars.values)) {
        if (pb.greater_constraints.empty() and pb.less_constraints.empty() and
            is_boolean_coefficient(pb.equal_constraints))
            return bx::itm::inequalities_Zcoeff_wedelin_solve(ctx, pb);

        if (is_101_coefficient(pb))
            return bx::itm::inequalities_Zcoeff_wedelin_solve(ctx, pb);

        return bx::itm::inequalities_Zcoeff_wedelin_solve(ctx, pb);
    }

    error(ctx, "no solver available for integer variable");

    throw bx::solver_failure(bx::solver_error_tag::no_solver_available);
}

bx::result
optimize(std::shared_ptr<bx::context> ctx, bx::problem& pb)
{
    auto th = get_thread_number(ctx);
    baryonyx_private::preprocess(ctx, pb);

    if (is_boolean_variable(pb.vars.values)) {
        if (pb.greater_constraints.empty() and pb.less_constraints.empty() and
            is_boolean_coefficient(pb.equal_constraints))
            return bx::itm::inequalities_Zcoeff_wedelin_optimize(ctx, pb, th);

        if (is_101_coefficient(pb))
            return bx::itm::inequalities_Zcoeff_wedelin_optimize(ctx, pb, th);

        return bx::itm::inequalities_Zcoeff_wedelin_optimize(ctx, pb, th);
    }

    error(ctx, "no optimizer available for integer variable");

    throw bx::solver_failure(bx::solver_error_tag::no_solver_available);
}
}
