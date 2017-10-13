/* Copyright (C) 2017 INRA
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublipnse, and/or sell copies of the Software, and to
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

namespace baryonyx {
namespace itm {

struct maximize_tag
{
};

struct minimize_tag
{
};

result
inequalities_1coeff_wedelin_solve(std::shared_ptr<context> ctx, problem& pb);

result
inequalities_1coeff_wedelin_optimize(std::shared_ptr<context> ctx,
                                     problem& pb,
                                     int thread);

enum class floating_point_type
{
    float_type,
    double_type,
    longdouble_type
};

inline const char*
floating_point_type_to_string(floating_point_type type) noexcept
{
    static const char* ret[] = {
        "float", "double", "longdouble",
    };

    switch (type) {
    case floating_point_type::float_type:
        return ret[0];
    case floating_point_type::double_type:
        return ret[1];
    case floating_point_type::longdouble_type:
        return ret[2];
    default:
        return ret[0];
    };
}

inline floating_point_type
get_floating_point_type(std::shared_ptr<context> ctx) noexcept
{
    auto str = ctx->get_string_parameter("floating-point-type", "double");

    if (str == "float")
        return floating_point_type::float_type;
    if (str == "longdouble")
        return floating_point_type::longdouble_type;

    return floating_point_type::double_type;
}

enum class constraint_order
{
    none,
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

    switch (type) {
    case constraint_order::none:
        return ret[0];
    case constraint_order::reversing:
        return ret[1];
    case constraint_order::random_sorting:
        return ret[2];
    case constraint_order::infeasibility_decr:
        return ret[3];
    case constraint_order::infeasibility_incr:
        return ret[4];
    default:
        return ret[0];
    }
}

inline constraint_order
get_constraint_order(std::shared_ptr<context> ctx) noexcept
{
    auto str = ctx->get_string_parameter("constraint-order", "none");

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
    parameters(std::shared_ptr<context> ctx)
      : preprocessing(ctx->get_string_parameter("preprocessing", "none"))
      , norm(ctx->get_string_parameter("norm", "inf"))
      , time_limit(ctx->get_real_parameter("time-limit", -1.0))
      , theta(ctx->get_real_parameter("theta", 0.5))
      , delta(ctx->get_real_parameter("delta", 0.01))
      , kappa_min(ctx->get_real_parameter("kappa-min", 0.0))
      , kappa_step(ctx->get_real_parameter("kappa-step", 1.e-3))
      , kappa_max(ctx->get_real_parameter("kappa-max", 0.6))
      , alpha(ctx->get_real_parameter("alpha", 1.0))
      , pushing_k_factor(ctx->get_real_parameter("pushing-k-factor", 0.9))
      , pushing_objective_amplifier(
          ctx->get_real_parameter("pushing-objective-amplifier", 5))
      , pushes_limit(ctx->get_integer_parameter("pushes-limit", 10))
      , pushing_iteration_limit(
          ctx->get_integer_parameter("pushing-iteration-limit", 20))
      , limit(ctx->get_integer_parameter("limit", 1000))
      , w(ctx->get_integer_parameter("w", 20))
      , serialize(ctx->get_integer_parameter("serialize", 0))
      , order(get_constraint_order(ctx))
      , float_type(get_floating_point_type(ctx))
    {
        if (limit < 0)
            limit = std::numeric_limits<int>::max();

        ctx->info("solver: inequalities_1coeff_wedelin\n"
                  "solver parameters:\n"
                  "  - preprocessing: %s\n"
                  "  - constraint-order: %s\n"
                  "  - time-limit: %.10g\n"
                  "  - theta: %.10g\n"
                  "  - delta: %.10g\n"
                  "  - limit: %d\n"
                  "  - kappa: %.10g %.10g %.10g\n"
                  "  - alpha: %.10g\n"
                  "  - w: %d\n"
                  "  - norm: %s\n"
                  "  - serialise: %d\n"
                  "  - float-point-type: %s\n",
                  preprocessing.c_str(),
                  constraint_order_to_string(order),
                  time_limit,
                  theta,
                  delta,
                  limit,
                  kappa_min,
                  kappa_step,
                  kappa_max,
                  alpha,
                  w,
                  norm.c_str(),
                  serialize,
                  floating_point_type_to_string(float_type));

        if (ctx->optimize())
            ctx->info("optimizer parameters:\n"
                      "  - pushed limit: %d\n"
                      "  - pushing objective amplifier: %.10g\n"
                      "  - pushing iteration limit: %d\n"
                      "  - pushing k factor: %.10g\n",
                      pushes_limit,
                      pushing_objective_amplifier,
                      pushing_iteration_limit,
                      pushing_k_factor);
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
    double pushing_k_factor;
    double pushing_objective_amplifier;
    int pushes_limit;
    int pushing_iteration_limit;
    int limit;
    int w;
    int serialize;
    constraint_order order;
    floating_point_type float_type;
};

struct merged_constraint
{
    merged_constraint(const std::vector<function_element>& elements_,
                      int min_,
                      int max_,
                      int id_)
      : elements(elements_)
      , min(min_)
      , max(max_)
      , id(id_)
    {
    }

    std::vector<function_element> elements;
    int min;
    int max;
    int id;
};

std::vector<merged_constraint>
make_merged_constraints(std::shared_ptr<context> ctx,
                        const problem& pb,
                        const parameters& p);
}
}

#endif
