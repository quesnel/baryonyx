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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_HPP

#include <baryonyx/core-compare>
#include <baryonyx/core-out>

#include "utils.hpp"

#include <algorithm>
#include <utility>

namespace baryonyx {
namespace itm {

enum class init_policy_type
{
    bastert = 0,
    random,
    best
};

inline const char*
init_policy_type_to_string(init_policy_type type) noexcept
{
    static const char* ret[] = { "bastert", "random", "best" };

    return ret[static_cast<int>(type)];
}

inline init_policy_type
get_init_policy_type(const context_ptr& ctx) noexcept
{
    auto str = context_get_string_parameter(ctx, "init-policy", "bastert");

    if (str == "random")
        return init_policy_type::random;

    if (str == "best")
        return init_policy_type::best;

    return init_policy_type::bastert;
}

enum class floating_point_type
{
    float_type = 0,
    double_type,
    longdouble_type
};

inline const char*
floating_point_type_to_string(floating_point_type type) noexcept
{
    static const char* ret[] = {
        "float",
        "double",
        "longdouble",
    };

    return ret[static_cast<int>(type)];
}

inline floating_point_type
get_floating_point_type(const context_ptr& ctx) noexcept
{
    auto str =
      context_get_string_parameter(ctx, "floating-point-type", "double");

    if (str == "float")
        return floating_point_type::float_type;
    if (str == "longdouble")
        return floating_point_type::longdouble_type;

    return floating_point_type::double_type;
}

enum class constraint_order
{
    none = 0,
    reversing,
    random_sorting,
    infeasibility_decr,
    infeasibility_incr,
};

inline const char*
constraint_order_to_string(constraint_order type) noexcept
{
    static const char* ret[] = {
        "none",
        "reversing",
        "random-sorting",
        "infeasibility-decr",
        "infeasibility-incr",
    };

    return ret[static_cast<int>(type)];
}

inline constraint_order
get_constraint_order(const context_ptr& ctx) noexcept
{
    auto str = context_get_string_parameter(ctx, "constraint-order", "none");

    if (str == "reversing")
        return constraint_order::reversing;
    if (str == "random-sorting")
        return constraint_order::random_sorting;
    if (str == "infeasibility-decr")
        return constraint_order::infeasibility_decr;
    if (str == "infeasibility-incr")
        return constraint_order::infeasibility_incr;

    return constraint_order::none;
}

struct parameters
{
    parameters(const context_ptr& ctx)
      : preprocessing(
          context_get_string_parameter(ctx, "preprocessing", "none"))
      , norm(context_get_string_parameter(ctx, "norm", "inf"))
      , time_limit(context_get_real_parameter(ctx, "time-limit", -1.0))
      , theta(context_get_real_parameter(ctx, "theta", 0.5))
      , delta(context_get_real_parameter(ctx, "delta", -1.0))
      , kappa_min(context_get_real_parameter(ctx, "kappa-min", 0.0))
      , kappa_step(context_get_real_parameter(ctx, "kappa-step", 1.e-3))
      , kappa_max(context_get_real_parameter(ctx, "kappa-max", 0.6))
      , alpha(context_get_real_parameter(ctx, "alpha", 1.0))
      , reverse_solution(
          context_get_real_parameter(ctx, "reverse-solution", -0.5))
      , pushing_k_factor(
          context_get_real_parameter(ctx, "pushing-k-factor", 0.9))
      , pushing_objective_amplifier(
          context_get_real_parameter(ctx, "pushing-objective-amplifier", 5))
      , init_random(context_get_real_parameter(ctx, "init-random", 0.5))
      , pushes_limit(context_get_integer_parameter(ctx, "pushes-limit", 100))
      , pushing_iteration_limit(
          context_get_integer_parameter(ctx, "pushing-iteration-limit", 50))
      , limit(context_get_integer_parameter(ctx, "limit", 1000))
      , w(context_get_integer_parameter(ctx, "w", 20))
      , print_level(context_get_integer_parameter(ctx, "print-level", 0))
      , order(get_constraint_order(ctx))
      , float_type(get_floating_point_type(ctx))
      , init_policy(get_init_policy_type(ctx))
    {
        if (limit < 0)
            limit = std::numeric_limits<int>::max();

        info(ctx,
             " * Global parameters:\n"
             "  - limit: {}\n"
             "  - time-limit: {:.10g}\n"
             "  - floating-point-type: {}\n"
             "  - print-level: {}\n",
             limit,
             time_limit,
             floating_point_type_to_string(float_type),
             print_level);

        info(ctx,
             " * In The Middle parameters:\n"
             "  - preprocessing: {}\n"
             "  - constraint-order: {}\n"
             "  - theta: {:.10g}\n"
             "  - delta: {:.10g}\n"
             "  - kappa: {:.10g} {:.10g} {:.10g}\n"
             "  - alpha: {:.10g}\n"
             "  - w: {}\n"
             "  - norm: {}\n",
             preprocessing,
             constraint_order_to_string(order),
             theta,
             delta,
             kappa_min,
             kappa_step,
             kappa_max,
             alpha,
             w,
             norm);

        info(ctx,
             " * Pushes system parameters:\n"
             "  - pushes-limit: {}\n"
             "  - pushing-objective-amplifier: {:.10g}\n"
             "  - pushing-iteration-limit: {}\n"
             "  - pushing-k-factor: {:.10g}\n",
             pushes_limit,
             pushing_objective_amplifier,
             pushing_iteration_limit,
             pushing_k_factor);

        info(ctx,
             " * Initialization parameters:\n"
             "  - init-policy: {}\n"
             "  - init-random: {:.10g}\n",
             init_policy_type_to_string(init_policy),
             init_random);
    }

    std::string preprocessing;
    std::string norm;
    double time_limit;
    double theta;
    double delta;
    double kappa_min;
    double kappa_step;
    double kappa_max;
    double alpha;
    double reverse_solution;
    double pushing_k_factor;
    double pushing_objective_amplifier;
    double init_random;
    int pushes_limit;
    int pushing_iteration_limit;
    int limit;
    int w;
    int print_level;
    constraint_order order;
    floating_point_type float_type;
    init_policy_type init_policy;
};

struct merged_constraint
{
    merged_constraint(std::vector<function_element> elements_,
                      int min_,
                      int max_,
                      int id_)
      : elements(std::move(elements_))
      , min(min_)
      , max(max_)
      , id(id_)
    {}

    std::vector<function_element> elements;
    int min;
    int max;
    int id;
};

std::vector<merged_constraint>
make_merged_constraints(const context_ptr& ctx,
                        const problem& pb,
                        const parameters& p);

result
solve(const baryonyx::context_ptr& ctx, problem& pb);

result
optimize(const baryonyx::context_ptr& ctx, problem& pb, int thread);
}
}

#endif
