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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_SOLVER_COMMON_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_SOLVER_COMMON_HPP

#include <baryonyx/core-compare>

#include <fmt/format.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <future>
#include <iterator>
#include <random>
#include <set>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include "itm-common.hpp"
#include "observer.hpp"
#include "private.hpp"
#include "sparse-matrix.hpp"
#include "utils.hpp"

#include <cassert>

namespace baryonyx {
namespace itm {

template<typename SolverT,
         typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT,
         typename observerT>
struct solver_functor
{
    using floatingpoint_type = floatingpointT;
    using mode_type = modeT;
    using constraint_order_type = constraintOrderT;
    using random_type = randomT;

    std::chrono::time_point<std::chrono::steady_clock> m_begin;
    std::chrono::time_point<std::chrono::steady_clock> m_end;

    const context_ptr& m_ctx;
    randomT& m_rng;

    result m_best;

    solver_functor(const context_ptr& ctx,
                   randomT& rng,
                   const std::vector<std::string>& variable_names,
                   const affected_variables& affected_vars)
      : m_ctx(ctx)
      , m_rng(rng)
    {
        m_best.affected_vars = affected_vars;
        m_best.variable_name = variable_names;
    }

    result operator()(const std::vector<merged_constraint>& constraints,
                      int variables,
                      const std::unique_ptr<floatingpointT[]>& original_costs,
                      const std::unique_ptr<floatingpointT[]>& norm_costs,
                      double cost_constant)
    {
        m_begin = std::chrono::steady_clock::now();
        m_end = m_begin;

        int i = 0;
        int pushed = -1;
        int best_remaining = -1;
        int pushing_iteration = m_ctx->parameters.pushing_iteration_limit;
        const auto& p = m_ctx->parameters;

        const auto kappa_min = static_cast<floatingpoint_type>(p.kappa_min);
        const auto kappa_step = static_cast<floatingpoint_type>(p.kappa_step);
        const auto kappa_max = static_cast<floatingpoint_type>(p.kappa_max);
        const auto alpha = static_cast<floatingpoint_type>(p.alpha);
        const auto theta = static_cast<floatingpoint_type>(p.theta);
        const auto delta = p.delta < 0
                             ? compute_delta<floatingpoint_type>(
                                 m_ctx, norm_costs, theta, variables)
                             : static_cast<floatingpoint_type>(p.delta);

        const auto pushing_k_factor =
          static_cast<floatingpoint_type>(p.pushing_k_factor);
        const auto pushing_objective_amplifier =
          static_cast<floatingpoint_type>(p.pushing_objective_amplifier);

        auto kappa = kappa_min;

        SolverT slv(m_rng,
                    length(constraints),
                    variables,
                    norm_costs,
                    constraints,
                    p.init_policy,
                    p.init_random);

        constraint_order_type compute(m_ctx, slv, m_rng);

        auto max_cost = max_cost_init(original_costs, variables, mode_type());
        bounds_printer<floatingpointT, modeT> bound_print(max_cost);
        observerT obs(slv, "img", p.limit);

        info(m_ctx, "* solver starts:\n");

        m_best.variables = slv.m;
        m_best.constraints = slv.n;

        bool start_pushing = false;

        for (;;) {
            int remaining = compute.run(slv, kappa, delta, theta);
            obs.make_observation();

            if (best_remaining == -1 or remaining < best_remaining) {
                best_remaining = remaining;
                m_best.duration = compute_duration(m_begin, m_end);
                m_best.loop = i;
                m_best.remaining_constraints = remaining;

                if (remaining == 0) {
                    m_best.status = baryonyx::result_status::success;
                    store_if_better(
                      slv.results(original_costs, cost_constant), slv.x, i);
                    start_pushing = true;
                } else {

                    // bound_print(slv, m_ctx, m_best);

                    info(m_ctx,
                         "  - violated constraints: {}/{} at {}s (loop: {})\n",
                         remaining,
                         slv.m,
                         compute_duration(m_begin, m_end),
                         i);
                }
            }

#ifndef BARYONYX_FULL_OPTIMIZATION
            print_solver(slv, m_ctx, m_best.variable_name, p.print_level);
#endif

            if (start_pushing) {
                ++pushing_iteration;

                if (pushed == -1)
                    info(m_ctx, "  - start push system:\n");

                if (pushing_iteration >= p.pushing_iteration_limit) {
                    pushed++;
                    pushing_iteration = 0;

                    remaining =
                      compute.push_and_run(slv,
                                           pushing_k_factor * kappa,
                                           delta,
                                           theta,
                                           pushing_objective_amplifier);

                    if (remaining == 0)
                        store_if_better(
                          slv.results(original_costs, cost_constant),
                          slv.x,
                          i);
                }

                if (pushed > p.pushes_limit) {
                    info(
                      m_ctx,
                      "    - Push system limit reached. Solution found: {}\n",
                      best_solution_value(m_best));
                    return m_best;
                }
            }

            if (i > p.w)
                kappa += kappa_step *
                         std::pow(static_cast<floatingpointT>(remaining) /
                                    static_cast<floatingpointT>(slv.m),
                                  alpha);

            ++i;
            if (p.limit > 0 and i > p.limit) {
                info(m_ctx, "  - Loop limit reached: {}\n", i);
                if (pushed == -1)
                    m_best.status = result_status::limit_reached;

                // if (context_get_integer_parameter(m_ctx, "print-level", 0) >
                // 0)
                //     print_missing_constraint(m_ctx,
                //                              slv.ap,
                //                              slv.A,
                //                              m_best.variable_value,
                //                              slv.b,
                //                              m_variable_names);

                return m_best;
            }

            if (kappa > kappa_max) {
                info(m_ctx, "  - Kappa max reached: {:+.6f}\n", kappa);
                if (pushed == -1)
                    m_best.status = result_status::kappa_max_reached;

                // if (context_get_integer_parameter(m_ctx, "print-level", 0) >
                // 0)
                //     print_missing_constraint(m_ctx,
                //                              slv.ap,
                //                              slv.A,
                //                              m_best.variable_value,
                //                              slv.b,
                //                              m_variable_names);

                return m_best;
            }

            m_end = std::chrono::steady_clock::now();
            if (is_time_limit(p.time_limit, m_begin, m_end)) {
                info(m_ctx, "  - Time limit reached: {} {:+.6f}\n", i, kappa);
                if (pushed == -1)
                    m_best.status = result_status::time_limit_reached;

                // if (context_get_integer_parameter(m_ctx, "print-level", 0) >
                // 0)
                //     print_missing_constraint(m_ctx,
                //                              slv.ap,
                //                              slv.A,
                //                              m_best.variable_value,
                //                              slv.b,
                //                              m_variable_names);

                return m_best;
            }
        }
    }

private:
    void store_if_better(double current, const x_type& x, int i)
    {
        if (m_best.solutions.empty() or
            is_better_solution(
              current, m_best.solutions.back().value, modeT())) {
            m_best.solutions.emplace_back(x, current);
            m_best.duration = compute_duration(m_begin, m_end);
            m_best.loop = i;

            info(m_ctx,
                 "  - Solution found: {} (i={} t={}s)\n",
                 current,
                 m_best.loop,
                 m_best.duration);
        } else {
            m_best.solutions.emplace_back(x, current);

            std::swap(*(m_best.solutions.rbegin()),
                      *(m_best.solutions.rbegin() + 1));
        }
    }
};

template<typename Enum>
constexpr std::size_t
enum_value(Enum e)
{
    return static_cast<std::size_t>(e);
};

template<typename Enum>
struct enum_size
{};

template<>
struct enum_size<solver_parameters::floating_point_type>
{
    static constexpr std::size_t size{ 3 };
};

template<>
struct enum_size<objective_function_type>
{
    static constexpr std::size_t size{ 2 };
};

template<>
struct enum_size<solver_parameters::constraint_order>
{
    static constexpr std::size_t size{ 5 };
};

template<typename Enum1, typename Enum2, typename Enum3>
constexpr std::size_t
combine_enums(Enum1 e1, Enum2 e2, Enum3 e3)
{
    return enum_value(e1) + enum_size<Enum1>::size * enum_value(e2) +
           enum_size<Enum1>::size * enum_size<Enum2>::size * enum_value(e3);
}

template<typename SolverT, typename constraintOrderT, typename randomT>
inline result
solve_problem(const context_ptr& ctx, const problem& pb);

/**
 * Template function to select the good template parameter for the underlying
 * itm::solver.
 */
template<template<typename floatingpointT, typename modeT, typename randomT>
         class Solver,
         typename random_type>
inline result
select_solver_parameters(const context_ptr& ctx, const problem& pb)
{
    switch (combine_enums(
      ctx->parameters.float_type, pb.type, ctx->parameters.order)) {
    case combine_enums(solver_parameters::floating_point_type::float_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::none):
        return solve_problem<Solver<float, minimize_tag, random_type>,
                             compute_none<float, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::double_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::none):
        return solve_problem<Solver<double, minimize_tag, random_type>,
                             compute_none<double, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::longdouble_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::none):
        return solve_problem<Solver<long double, minimize_tag, random_type>,
                             compute_none<long double, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::float_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::none):
        return solve_problem<Solver<float, maximize_tag, random_type>,
                             compute_none<float, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::double_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::none):
        return solve_problem<Solver<double, maximize_tag, random_type>,
                             compute_none<double, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::longdouble_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::none):
        return solve_problem<Solver<long double, minimize_tag, random_type>,
                             compute_none<long double, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::float_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::reversing):
        return solve_problem<Solver<float, minimize_tag, random_type>,
                             compute_reversing<float, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::double_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::reversing):
        return solve_problem<Solver<double, minimize_tag, random_type>,
                             compute_reversing<double, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::longdouble_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::reversing):
        return solve_problem<Solver<long double, minimize_tag, random_type>,
                             compute_reversing<long double, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::float_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::reversing):
        return solve_problem<Solver<float, maximize_tag, random_type>,
                             compute_reversing<float, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::double_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::reversing):
        return solve_problem<Solver<double, maximize_tag, random_type>,
                             compute_reversing<double, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::longdouble_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::reversing):
        return solve_problem<Solver<long double, minimize_tag, random_type>,
                             compute_reversing<long double, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::float_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::random_sorting):
        return solve_problem<Solver<float, minimize_tag, random_type>,
                             compute_random<float, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::double_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::random_sorting):
        return solve_problem<Solver<double, minimize_tag, random_type>,
                             compute_random<double, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::longdouble_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::random_sorting):
        return solve_problem<Solver<long double, minimize_tag, random_type>,
                             compute_random<long double, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::float_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::random_sorting):
        return solve_problem<Solver<float, maximize_tag, random_type>,
                             compute_random<float, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::double_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::random_sorting):
        return solve_problem<Solver<double, maximize_tag, random_type>,
                             compute_random<double, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::longdouble_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::random_sorting):
        return solve_problem<Solver<long double, minimize_tag, random_type>,
                             compute_random<long double, random_type>,
                             random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::float_type,
      baryonyx::objective_function_type::minimize,
      solver_parameters::constraint_order::infeasibility_decr):
        return solve_problem<Solver<float, minimize_tag, random_type>,
                             compute_infeasibility<float,
                                                   random_type,
                                                   compute_infeasibility_decr>,
                             random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::double_type,
      baryonyx::objective_function_type::minimize,
      solver_parameters::constraint_order::infeasibility_decr):
        return solve_problem<Solver<double, minimize_tag, random_type>,
                             compute_infeasibility<double,
                                                   random_type,
                                                   compute_infeasibility_decr>,
                             random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::longdouble_type,
      baryonyx::objective_function_type::minimize,
      solver_parameters::constraint_order::infeasibility_decr):
        return solve_problem<Solver<long double, minimize_tag, random_type>,
                             compute_infeasibility<long double,
                                                   random_type,
                                                   compute_infeasibility_decr>,
                             random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::float_type,
      baryonyx::objective_function_type::maximize,
      solver_parameters::constraint_order::infeasibility_decr):
        return solve_problem<Solver<float, maximize_tag, random_type>,
                             compute_infeasibility<float,
                                                   random_type,
                                                   compute_infeasibility_decr>,
                             random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::double_type,
      baryonyx::objective_function_type::maximize,
      solver_parameters::constraint_order::infeasibility_decr):
        return solve_problem<Solver<double, maximize_tag, random_type>,
                             compute_infeasibility<double,
                                                   random_type,
                                                   compute_infeasibility_decr>,
                             random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::longdouble_type,
      baryonyx::objective_function_type::maximize,
      solver_parameters::constraint_order::infeasibility_decr):
        return solve_problem<Solver<long double, minimize_tag, random_type>,
                             compute_infeasibility<long double,
                                                   random_type,
                                                   compute_infeasibility_decr>,
                             random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::float_type,
      baryonyx::objective_function_type::minimize,
      solver_parameters::constraint_order::infeasibility_incr):
        return solve_problem<Solver<float, minimize_tag, random_type>,
                             compute_infeasibility<float,
                                                   random_type,
                                                   compute_infeasibility_incr>,
                             random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::double_type,
      baryonyx::objective_function_type::minimize,
      solver_parameters::constraint_order::infeasibility_incr):
        return solve_problem<Solver<double, minimize_tag, random_type>,
                             compute_infeasibility<double,
                                                   random_type,
                                                   compute_infeasibility_incr>,
                             random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::longdouble_type,
      baryonyx::objective_function_type::minimize,
      solver_parameters::constraint_order::infeasibility_incr):
        return solve_problem<Solver<long double, minimize_tag, random_type>,
                             compute_infeasibility<long double,
                                                   random_type,
                                                   compute_infeasibility_incr>,
                             random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::float_type,
      baryonyx::objective_function_type::maximize,
      solver_parameters::constraint_order::infeasibility_incr):
        return solve_problem<Solver<float, maximize_tag, random_type>,
                             compute_infeasibility<float,
                                                   random_type,
                                                   compute_infeasibility_incr>,
                             random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::double_type,
      baryonyx::objective_function_type::maximize,
      solver_parameters::constraint_order::infeasibility_incr):
        return solve_problem<Solver<double, maximize_tag, random_type>,
                             compute_infeasibility<double,
                                                   random_type,
                                                   compute_infeasibility_incr>,
                             random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::longdouble_type,
      baryonyx::objective_function_type::maximize,
      solver_parameters::constraint_order::infeasibility_incr):
        return solve_problem<Solver<long double, minimize_tag, random_type>,
                             compute_infeasibility<long double,
                                                   random_type,
                                                   compute_infeasibility_incr>,
                             random_type>(ctx, pb);
    }

    return result(result_status::internal_error);
}

template<typename SolverT, typename constraintOrderT, typename randomT>
inline result
optimize_problem(const context_ptr& ctx, const problem& pb);

/**
 * Template function to select the good template parameter for the underlying
 * itm::solver.
 */
template<template<typename floatingpointT, typename modeT, typename randomT>
         class Solver,
         typename random_type>
auto
select_optimizer_parameters(const context_ptr& ctx, const problem& pb)
{
    switch (combine_enums(
      ctx->parameters.float_type, pb.type, ctx->parameters.order)) {
    case combine_enums(solver_parameters::floating_point_type::float_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::none):
        return optimize_problem<Solver<float, minimize_tag, random_type>,
                                compute_none<float, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::double_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::none):
        return optimize_problem<Solver<double, minimize_tag, random_type>,
                                compute_none<double, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::longdouble_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::none):
        return optimize_problem<Solver<long double, minimize_tag, random_type>,
                                compute_none<long double, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::float_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::none):
        return optimize_problem<Solver<float, maximize_tag, random_type>,
                                compute_none<float, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::double_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::none):
        return optimize_problem<Solver<double, maximize_tag, random_type>,
                                compute_none<double, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::longdouble_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::none):
        return optimize_problem<Solver<long double, minimize_tag, random_type>,
                                compute_none<long double, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::float_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::reversing):
        return optimize_problem<Solver<float, minimize_tag, random_type>,
                                compute_reversing<float, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::double_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::reversing):
        return optimize_problem<Solver<double, minimize_tag, random_type>,
                                compute_reversing<double, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::longdouble_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::reversing):
        return optimize_problem<Solver<long double, minimize_tag, random_type>,
                                compute_reversing<long double, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::float_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::reversing):
        return optimize_problem<Solver<float, maximize_tag, random_type>,
                                compute_reversing<float, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::double_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::reversing):
        return optimize_problem<Solver<double, maximize_tag, random_type>,
                                compute_reversing<double, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::longdouble_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::reversing):
        return optimize_problem<Solver<long double, minimize_tag, random_type>,
                                compute_reversing<long double, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::float_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::random_sorting):
        return optimize_problem<Solver<float, minimize_tag, random_type>,
                                compute_random<float, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::double_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::random_sorting):
        return optimize_problem<Solver<double, minimize_tag, random_type>,
                                compute_random<double, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::longdouble_type,
                       baryonyx::objective_function_type::minimize,
                       solver_parameters::constraint_order::random_sorting):
        return optimize_problem<Solver<long double, minimize_tag, random_type>,
                                compute_random<long double, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::float_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::random_sorting):
        return optimize_problem<Solver<float, maximize_tag, random_type>,
                                compute_random<float, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::double_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::random_sorting):
        return optimize_problem<Solver<double, maximize_tag, random_type>,
                                compute_random<double, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(solver_parameters::floating_point_type::longdouble_type,
                       baryonyx::objective_function_type::maximize,
                       solver_parameters::constraint_order::random_sorting):
        return optimize_problem<Solver<long double, minimize_tag, random_type>,
                                compute_random<long double, random_type>,
                                random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::float_type,
      baryonyx::objective_function_type::minimize,
      solver_parameters::constraint_order::infeasibility_decr):
        return optimize_problem<
          Solver<float, minimize_tag, random_type>,
          compute_infeasibility<float,
                                random_type,
                                compute_infeasibility_decr>,
          random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::double_type,
      baryonyx::objective_function_type::minimize,
      solver_parameters::constraint_order::infeasibility_decr):
        return optimize_problem<
          Solver<double, minimize_tag, random_type>,
          compute_infeasibility<double,
                                random_type,
                                compute_infeasibility_decr>,
          random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::longdouble_type,
      baryonyx::objective_function_type::minimize,
      solver_parameters::constraint_order::infeasibility_decr):
        return optimize_problem<
          Solver<long double, minimize_tag, random_type>,
          compute_infeasibility<long double,
                                random_type,
                                compute_infeasibility_decr>,
          random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::float_type,
      baryonyx::objective_function_type::maximize,
      solver_parameters::constraint_order::infeasibility_decr):
        return optimize_problem<
          Solver<float, maximize_tag, random_type>,
          compute_infeasibility<float,
                                random_type,
                                compute_infeasibility_decr>,
          random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::double_type,
      baryonyx::objective_function_type::maximize,
      solver_parameters::constraint_order::infeasibility_decr):
        return optimize_problem<
          Solver<double, maximize_tag, random_type>,
          compute_infeasibility<double,
                                random_type,
                                compute_infeasibility_decr>,
          random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::longdouble_type,
      baryonyx::objective_function_type::maximize,
      solver_parameters::constraint_order::infeasibility_decr):
        return optimize_problem<
          Solver<long double, minimize_tag, random_type>,
          compute_infeasibility<long double,
                                random_type,
                                compute_infeasibility_decr>,
          random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::float_type,
      baryonyx::objective_function_type::minimize,
      solver_parameters::constraint_order::infeasibility_incr):
        return optimize_problem<
          Solver<float, minimize_tag, random_type>,
          compute_infeasibility<float,
                                random_type,
                                compute_infeasibility_incr>,
          random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::double_type,
      baryonyx::objective_function_type::minimize,
      solver_parameters::constraint_order::infeasibility_incr):
        return optimize_problem<
          Solver<double, minimize_tag, random_type>,
          compute_infeasibility<double,
                                random_type,
                                compute_infeasibility_incr>,
          random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::longdouble_type,
      baryonyx::objective_function_type::minimize,
      solver_parameters::constraint_order::infeasibility_incr):
        return optimize_problem<
          Solver<long double, minimize_tag, random_type>,
          compute_infeasibility<long double,
                                random_type,
                                compute_infeasibility_incr>,
          random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::float_type,
      baryonyx::objective_function_type::maximize,
      solver_parameters::constraint_order::infeasibility_incr):
        return optimize_problem<
          Solver<float, maximize_tag, random_type>,
          compute_infeasibility<float,
                                random_type,
                                compute_infeasibility_incr>,
          random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::double_type,
      baryonyx::objective_function_type::maximize,
      solver_parameters::constraint_order::infeasibility_incr):
        return optimize_problem<
          Solver<double, maximize_tag, random_type>,
          compute_infeasibility<double,
                                random_type,
                                compute_infeasibility_incr>,
          random_type>(ctx, pb);

    case combine_enums(
      solver_parameters::floating_point_type::longdouble_type,
      baryonyx::objective_function_type::maximize,
      solver_parameters::constraint_order::infeasibility_incr):
        return optimize_problem<
          Solver<long double, minimize_tag, random_type>,
          compute_infeasibility<long double,
                                random_type,
                                compute_infeasibility_incr>,
          random_type>(ctx, pb);
    }

    return result(result_status::internal_error);
}

result
solve_equalities_01coeff(const context_ptr& ctx, const problem& pb);

result
solve_equalities_101coeff(const context_ptr& ctx, const problem& pb);

result
solve_inequalities_01coeff(const context_ptr& ctx, const problem& pb);

result
solve_inequalities_101coeff(const context_ptr& ctx, const problem& pb);

result
solve_inequalities_101coeff_buffered(const context_ptr& ctx,
                                     const problem& pb);

result
solve_inequalities_Zcoeff(const context_ptr& ctx, const problem& pb);

template<typename floatingpointT, typename modeT>
struct best_solution_recorder
{
    const context_ptr& m_ctx;
    std::mutex m_mutex;
    result m_best;

    best_solution_recorder(const context_ptr& ctx)
      : m_ctx(ctx)
    {}

    bool try_update(const result& current) noexcept
    {
        try {
            std::lock_guard<std::mutex> lock(m_mutex);

            if (current.status != result_status::success)
                return false;

            if (m_best.status != result_status::success or
                is_better_solution(current, m_best, modeT())) {

                info(m_ctx,
                     "  - Solution found: {} (i={} t={}s)\n",
                     best_solution_value(current),
                     current.loop,
                     current.duration);

                m_best = current;

                return true;
            }
        } catch (const std::exception& e) {
            error(m_ctx, "sync optimization error: {}", e.what());
        }

        return false;
    }
};

template<typename SolverT,
         typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT>
struct optimize_functor
{
    using floatingpoint_type = floatingpointT;
    using mode_type = modeT;
    using constraint_order_type = constraintOrderT;
    using random_type = randomT;

    std::chrono::time_point<std::chrono::steady_clock> m_begin;
    std::chrono::time_point<std::chrono::steady_clock> m_end;

    std::set<solution> m_all_solutions;

    const context_ptr& m_ctx;
    randomT m_rng;
    int m_thread_id;

    result m_best;

    optimize_functor(const context_ptr& ctx,
                     unsigned thread_id,
                     typename random_type::result_type seed,
                     const std::vector<std::string>& variable_names,
                     const affected_variables& affected_vars)
      : m_ctx(ctx)
      , m_rng(seed)
      , m_thread_id(thread_id)
    {
        m_best.affected_vars = affected_vars;
        m_best.variable_name = variable_names;
    }

    result operator()(
      best_solution_recorder<floatingpointT, modeT>& best_recorder,
      const std::vector<merged_constraint>& constraints,
      int variables,
      const std::unique_ptr<floatingpointT[]>& original_costs,
      const std::unique_ptr<floatingpointT[]>& norm_costs,
      double cost_constant)
    {
        m_begin = std::chrono::steady_clock::now();
        m_end = m_begin;

        int i = 0;
        int pushed = -1;
        int pushing_iteration = 0;
        const auto& p = m_ctx->parameters;

        auto init_policy = p.init_policy;
        const auto kappa_min = static_cast<floatingpoint_type>(p.kappa_min);
        const auto kappa_step = static_cast<floatingpoint_type>(p.kappa_step);
        const auto kappa_max = static_cast<floatingpoint_type>(p.kappa_max);
        const auto alpha = static_cast<floatingpoint_type>(p.alpha);
        const auto theta = static_cast<floatingpoint_type>(p.theta);
        const auto delta = p.delta < 0
                             ? compute_delta<floatingpoint_type>(
                                 m_ctx, norm_costs, theta, variables)
                             : static_cast<floatingpoint_type>(p.delta);

        const auto pushing_k_factor =
          static_cast<floatingpoint_type>(p.pushing_k_factor);
        const auto pushing_objective_amplifier =
          static_cast<floatingpoint_type>(p.pushing_objective_amplifier);

        auto kappa = kappa_min;

        SolverT slv(m_rng,
                    length(constraints),
                    variables,
                    norm_costs,
                    constraints,
                    p.init_policy,
                    p.init_random);

        constraint_order_type compute(m_ctx, slv, m_rng);

        for (; not is_time_limit(p.time_limit, m_begin, m_end);
             m_end = std::chrono::steady_clock::now(), ++i) {

            int remaining = compute.run(slv, kappa, delta, theta);
            if (remaining == 0)
                store_if_better(slv.results(original_costs, cost_constant),
                                slv.x,
                                i,
                                best_recorder);

            if (i > p.w)
                kappa += kappa_step *
                         std::pow(static_cast<floatingpoint_type>(remaining) /
                                    static_cast<floatingpoint_type>(slv.m),
                                  alpha);

            if ((p.limit > 0 and i >= p.limit) or kappa > kappa_max or
                pushed > p.pushes_limit) {
                if (m_best.solutions.empty())
                    init_policy =
                      init_solver(slv, x_type(), init_policy, p.init_random);
                else
                    init_policy =
                      init_solver(slv,
                                  m_best.solutions.back().variables,
                                  init_policy,
                                  p.init_random);

                i = 0;
                kappa = static_cast<floatingpoint_type>(kappa_min);
                pushed = -1;
                pushing_iteration = 0;

                continue;
            }

            if (pushed >= 0) {
                ++pushing_iteration;

                if (pushing_iteration >= p.pushing_iteration_limit) {
                    pushed++;
                    pushing_iteration = 0;

                    remaining =
                      compute.push_and_run(slv,
                                           pushing_k_factor * kappa,
                                           delta,
                                           theta,
                                           pushing_objective_amplifier);

                    if (remaining == 0)
                        store_if_better(
                          slv.results(original_costs, cost_constant),
                          slv.x,
                          i,
                          best_recorder);
                }
            }
        }

        std::copy(m_all_solutions.begin(),
                  m_all_solutions.end(),
                  std::back_inserter(m_best.solutions));

        return m_best;
    }

private:
    void store_if_better(
      double current,
      const x_type& x,
      int i,
      best_solution_recorder<floatingpointT, modeT>& best_recorder)
    {
        if (m_best.solutions.empty() or
            is_better_solution(
              current, m_best.solutions.back().value, modeT())) {
            m_best.status = baryonyx::result_status::success;

            // Store only the best solution, other solutions are stored into
            // the @c std::set to avoid duplicated solutions.

            if (m_best.solutions.empty())
                m_best.solutions.emplace_back(x, current);
            else
                m_best.solutions[0] = { x, current };

            m_best.duration = compute_duration(m_begin, m_end);
            m_best.loop = i;

            best_recorder.try_update(m_best);
        }

        m_all_solutions.emplace(x, current);
    }
};

template<typename SolverT, typename constraintOrderT, typename randomT>
inline result
solve_problem(const context_ptr& ctx, const problem& pb)
{
    info(ctx, "Solver initializing\n");

    using floatingpoint_type = typename SolverT::floatingpoint_type;
    using mode_type = typename SolverT::mode_type;

    result ret;
    auto affected_vars = std::move(pb.affected_vars);

    auto constraints{ make_merged_constraints(ctx, pb) };
    if (not constraints.empty() and not pb.vars.values.empty()) {
        randomT rng(init_random_generator_seed<randomT>(ctx));

        auto variables = numeric_cast<int>(pb.vars.values.size());
        auto cost =
          make_objective_function<floatingpoint_type>(pb.objective, variables);
        auto norm_costs = normalize_costs<floatingpoint_type, randomT>(
          ctx, cost, rng, variables);
        auto cost_constant = pb.objective.value;
        auto names = std::move(pb.vars.names);

        switch (ctx->parameters.observer) {
        case solver_parameters::observer_type::pnm: {
            using obs = pnm_observer<SolverT, floatingpoint_type>;

            solver_functor<SolverT,
                           floatingpoint_type,
                           mode_type,
                           constraintOrderT,
                           randomT,
                           obs>
              slv(ctx, rng, names, affected_vars);

            ret = slv(constraints, variables, cost, norm_costs, cost_constant);
        } break;
        case solver_parameters::observer_type::file: {
            using obs = file_observer<SolverT, floatingpoint_type>;

            solver_functor<SolverT,
                           floatingpoint_type,
                           mode_type,
                           constraintOrderT,
                           randomT,
                           obs>
              slv(ctx, rng, names, affected_vars);

            ret = slv(constraints, variables, cost, norm_costs, cost_constant);
        } break;
        default: {
            using obs = none_observer<SolverT, floatingpoint_type>;

            solver_functor<SolverT,
                           floatingpoint_type,
                           mode_type,
                           constraintOrderT,
                           randomT,
                           obs>
              slv(ctx, rng, names, affected_vars);

            ret = slv(constraints, variables, cost, norm_costs, cost_constant);
            break;
        }
        }

        ret.method = "inequalities_Zcoeff solver";
        ret.variable_name = std::move(names);
    } else {
        ret.status = result_status::success;
    }
    ret.affected_vars = std::move(affected_vars);

    return ret;
}

template<typename SolverT, typename constraintOrderT, typename randomT>
inline result
optimize_problem(const context_ptr& ctx, const problem& pb)
{
    using floatingpoint_type = typename SolverT::floatingpoint_type;
    using mode_type = typename SolverT::mode_type;

    info(ctx, "Optimizer initializing\n");
    result ret;
    auto affected_vars = std::move(pb.affected_vars);

    auto constraints{ make_merged_constraints(ctx, pb) };
    if (not constraints.empty() and not pb.vars.values.empty()) {
        randomT rng(init_random_generator_seed<randomT>(ctx));

        auto variables = numeric_cast<int>(pb.vars.values.size());
        auto cost =
          make_objective_function<floatingpoint_type>(pb.objective, variables);
        auto norm_costs = normalize_costs<floatingpoint_type, randomT>(
          ctx, cost, rng, variables);
        auto cost_constant = pb.objective.value;
        auto names = std::move(pb.vars.names);

        const auto thread = static_cast<unsigned>(ctx->parameters.thread);
        std::vector<std::thread> pool(ctx->parameters.thread);
        pool.clear();
        std::vector<std::future<result>> results(thread);
        results.clear();

        best_solution_recorder<floatingpoint_type, mode_type> result(ctx);

        if (thread == 1u)
            info(ctx, "optimizer starts with one thread\n");
        else
            info(ctx, "Optimizer starts with {} threads\n", thread);

        auto seeds = generate_seed(rng, thread);

        for (unsigned i = 0; i != thread; ++i) {
            std::packaged_task<baryonyx::result()> task(
              std::bind(optimize_functor<SolverT,
                                         floatingpoint_type,
                                         mode_type,
                                         constraintOrderT,
                                         randomT>(
                          ctx, i, seeds[i], names, affected_vars),
                        std::ref(result),
                        std::ref(constraints),
                        variables,
                        std::ref(cost),
                        std::ref(norm_costs),
                        cost_constant));

            results.emplace_back(task.get_future());

            pool.emplace_back(std::thread(std::move(task)));
        }

        for (auto& t : pool)
            t.join();

        std::set<solution> all_solutions;

        for (unsigned i = 0; i != thread; ++i) {
            auto current = results[i].get();
            if (current.status == result_status::success) {
                all_solutions.insert(current.solutions.begin(),
                                     current.solutions.end());

                if (ret.solutions.empty() or
                    is_better_solution(current.solutions.back().value,
                                       ret.solutions.back().value,
                                       mode_type()))
                    ret = current;
            }
        }

        ret.solutions.clear();

        std::copy(all_solutions.begin(),
                  all_solutions.end(),
                  std::back_inserter(ret.solutions));
    }

    return ret;
}

result
optimize_equalities_01coeff(const context_ptr& ctx, const problem& pb);

result
optimize_equalities_101coeff(const context_ptr& ctx, const problem& pb);

result
optimize_inequalities_01coeff(const context_ptr& ctx, const problem& pb);

result
optimize_inequalities_101coeff(const context_ptr& ctx, const problem& pb);

result
optimize_inequalities_101coeff_buffered(const context_ptr& ctx,
                                        const problem& pb);

result
optimize_inequalities_Zcoeff(const context_ptr& ctx, const problem& pb);

result
solve(const context_ptr& ctx, const problem& pb);

result
optimize(const context_ptr& ctx, const problem& pb);

/**
 * @brief Auto-tune baryonyx solver parameters to optimize the problem.
 *
 * @details The @c nlopt library or a manual test is used to found the best
 *     parameters @c baryonyx::optimize function for the specified problem.
 *
 * @param ctx A context.
 * @param pb Problem definition.
 * @return A representation of the result.
 * @throw @c baryonyx::solver_error
 */
baryonyx::result
automatic_optimizer(const baryonyx::context_ptr& ctx,
                    const baryonyx::problem& pb);

} // namespace itm
} // namespace baryonyx

#endif
