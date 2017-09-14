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

#include <iterator>

#include "generalized-wedelin.hpp"
#include "inequalities-1coeff.hpp"
#include "mitm.hpp"
#include "utils.hpp"
#include "wedelin.hpp"

namespace baryonyx {

/**
 * Get number of thread to use in optimizer from parameters list. If an
 * error occured, this function returns 1.
 *
 * @param ctx An `baryonyx::context` where parameter "thread" is trying to be
 * read.
 *
 * @return An integer >= 1 if the value exist and can be convert to
 * positive integer or 1 if an error occurred.
 */
inline long int
get_thread_number(std::shared_ptr<baryonyx::context>& ctx) noexcept
{
    auto t = ctx->get_integer_parameter("thread",
                                        std::thread::hardware_concurrency());

    return t <= 0 ? 1 : baryonyx::numeric_cast<long int>(t);
}

std::tuple<double, double, double, long>
get_parameters(std::shared_ptr<baryonyx::context>& ctx) noexcept
{
    auto kappa = ctx->get_real_parameter("kappa", 0.001);
    auto theta = ctx->get_real_parameter("theta", 0.001);
    auto delta = ctx->get_real_parameter("delta", 0.001);
    auto limit = ctx->get_real_parameter("limit", 1000l);

    return std::make_tuple(kappa, delta, theta, limit);
}

template <typename variableT>
bool
is_boolean_variable(const variableT& vars)
{
    for (const auto& elem : vars)
        if (elem.type != variable_type::binary)
            return false;

    return true;
}

template <typename constraintsT>
bool
is_boolean_coefficient(const constraintsT csts)
{
    for (const auto& cst : csts)
        for (const auto& elem : cst.elements)
            if (elem.factor < 0 or elem.factor > 1)
                return false;

    return true;
}

template <typename variableT>
bool
is_integer_variable(const variableT& vars)
{
    for (const auto& elem : vars)
        if (elem.type != variable_type::general)
            return false;

    return true;
}

template <typename constraintsT>
bool
is_101_coefficient(const constraintsT csts)
{
    for (const auto& cst : csts)
        for (const auto& elem : cst.elements)
            if (elem.factor < -1 or elem.factor > 1)
                return false;

    return true;
}

template <typename functionT>
functionT
cleanup_function_element(functionT fct, long int& nb)
{
    if (fct.size() <= 1)
        return fct;

    std::sort(fct.begin(), fct.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.variable_index < rhs.variable_index;
    });

    functionT ret;
    auto prev = fct.begin();
    auto it = prev + 1;
    auto end = fct.end();

    ret.emplace_back(*prev);

    while (it != end) {
        if (it->variable_index == prev->variable_index) {
            Expects(ret.back().variable_index == it->variable_index);
            ret.back().factor += it->factor;
        } else {
            prev = it;
            ret.emplace_back(*it);
        }
        ++it;
    }

    {
        nb += numeric_cast<long int>(fct.size() - ret.size());
    }

    return ret;
}

void
cleanup_problem(std::shared_ptr<baryonyx::context> ctx, problem& pb)
{
    ctx->info("* cleaning problem definition:\n");

    long int nb_function_clean{ 0 };

    for (auto& elem : pb.equal_constraints)
        elem.elements =
          cleanup_function_element(elem.elements, nb_function_clean);
    for (auto& elem : pb.greater_constraints)
        elem.elements =
          cleanup_function_element(elem.elements, nb_function_clean);
    for (auto& elem : pb.greater_equal_constraints)
        elem.elements =
          cleanup_function_element(elem.elements, nb_function_clean);
    for (auto& elem : pb.less_constraints)
        elem.elements =
          cleanup_function_element(elem.elements, nb_function_clean);
    for (auto& elem : pb.less_equal_constraints)
        elem.elements =
          cleanup_function_element(elem.elements, nb_function_clean);

    ctx->info(" -> cleaned functions: %ld\n", nb_function_clean);
}

result
mitm_solve(std::shared_ptr<baryonyx::context> ctx, problem& pb)
{
    cleanup_problem(ctx, pb);

    if (pb.greater_constraints.empty() and
        pb.greater_equal_constraints.empty() and
        pb.less_constraints.empty() and pb.less_equal_constraints.empty() and
        is_boolean_coefficient(pb.equal_constraints) and
        is_boolean_variable(pb.vars.values)) {

        return baryonyx::inequalities_1coeff_wedelin_solve(ctx, pb);
    }

    if ((not pb.equal_constraints.empty() or
         not pb.greater_equal_constraints.empty() or
         not pb.less_equal_constraints.empty())) {

        for (int i{ 0 },
             e(numeric_cast<int>(pb.less_equal_constraints.size()));
             i != e;
             ++i)
            for (const auto& elem : pb.less_equal_constraints[i].elements)
                if (elem.factor < -1 or elem.factor > 1)
                    ctx->warning(
                      "Error at constraints %d (factor=%d)\n", i, elem.factor);

        if (pb.greater_constraints.empty() and pb.less_constraints.empty() and
            is_101_coefficient(pb.equal_constraints) and
            is_101_coefficient(pb.greater_equal_constraints) and
            is_101_coefficient(pb.less_equal_constraints) and
            is_boolean_variable(pb.vars.values)) {

            return baryonyx::inequalities_1coeff_wedelin_solve(ctx, pb);
        }
    }

    if ((is_101_coefficient(pb.equal_constraints) or
         is_101_coefficient(pb.greater_constraints) or
         is_101_coefficient(pb.greater_equal_constraints) or
         is_101_coefficient(pb.less_constraints) or
         is_101_coefficient(pb.less_equal_constraints)) and
        is_integer_variable(pb.vars.values)) {
        double kappa, delta, theta;
        long limit;

        std::tie(kappa, delta, theta, limit) = get_parameters(ctx);

        return baryonyx::generalized_wedelin(kappa, delta, theta, limit, pb);
    }

    throw baryonyx::solver_failure(solver_error_tag::no_solver_available);
}

result
mitm_optimize(std::shared_ptr<baryonyx::context> ctx, problem& pb)
{
    auto thread = get_thread_number(ctx);
    cleanup_problem(ctx, pb);

    if (pb.greater_constraints.empty() and
        pb.greater_equal_constraints.empty() and
        pb.less_constraints.empty() and pb.less_equal_constraints.empty() and
        is_boolean_coefficient(pb.equal_constraints) and
        is_boolean_variable(pb.vars.values)) {

        return baryonyx::inequalities_1coeff_wedelin_optimize(ctx, pb, thread);
    }

    if ((not pb.equal_constraints.empty() or
         not pb.greater_equal_constraints.empty() or
         not pb.less_equal_constraints.empty())) {

        for (int i{ 0 },
             e(numeric_cast<int>(pb.less_equal_constraints.size()));
             i != e;
             ++i)
            for (const auto& elem : pb.less_equal_constraints[i].elements)
                if (elem.factor < -1 or elem.factor > 1)
                    ctx->warning("Error %d (factor=%d)\n", i, elem.factor);

        if (pb.greater_constraints.empty() and pb.less_constraints.empty() and
            is_101_coefficient(pb.equal_constraints) and
            is_101_coefficient(pb.greater_equal_constraints) and
            is_101_coefficient(pb.less_equal_constraints) and
            is_boolean_variable(pb.vars.values)) {

            return baryonyx::inequalities_1coeff_wedelin_optimize(
              ctx, pb, thread);
        }
    }

    if ((is_101_coefficient(pb.equal_constraints) or
         is_101_coefficient(pb.greater_constraints) or
         is_101_coefficient(pb.greater_equal_constraints) or
         is_101_coefficient(pb.less_constraints) or
         is_101_coefficient(pb.less_equal_constraints)) and
        is_integer_variable(pb.vars.values)) {
        double kappa, delta, theta;
        long limit;

        std::tie(kappa, delta, theta, limit) = get_parameters(ctx);

        return baryonyx::generalized_wedelin(kappa, delta, theta, limit, pb);
    }

    throw baryonyx::solver_failure(solver_error_tag::no_solver_available);
}

} // namespace baryonyx
