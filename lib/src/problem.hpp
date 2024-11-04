/* Copyright (C) 2016-2021 INRAE
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software && associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, && to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice && this permission notice shall be included
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

#ifndef ORG_VLEPROJECT_BARYONYX_LIB_PROBLEM_HPP
#define ORG_VLEPROJECT_BARYONYX_LIB_PROBLEM_HPP

#include <baryonyx/core-compare>
#include <baryonyx/core>

#include <fmt/format.h>

#include <algorithm>

namespace baryonyx {

/**
 * @brief Internal representation of a problem with preprocessed data.
 *
 * @details A equivalent class to @c raw_problem but with preprocessed
 *     structures to store affected variables and duplicated constraints from
 *     the original problem.
 *
 */
struct problem
{
    problem() = default;

    problem(const raw_problem& pb)
      : strings(pb.strings)
      , objective(pb.objective)
      , equal_constraints(pb.equal_constraints)
      , greater_constraints(pb.greater_constraints)
      , less_constraints(pb.less_constraints)
      , vars(pb.vars)
      , type(pb.type)
    {
        problem_type = which_problem_type();
    }

    problem(raw_problem&& pb)
      : strings(pb.strings)
      , objective(pb.objective)
      , equal_constraints(pb.equal_constraints)
      , greater_constraints(pb.greater_constraints)
      , less_constraints(pb.less_constraints)
      , vars(pb.vars)
      , type(pb.type)
    {
        problem_type = which_problem_type();
    }

    template<typename Constraints>
    constexpr static int coefficient_type(const Constraints& csts) noexcept
    {
        int ret = 0;

        for (const auto& cst : csts) {
            for (const auto& elem : cst.elements) {
                if (elem.factor == -1 && ret == 0)
                    ret = 1;
                else if (elem.factor <= -2 || elem.factor >= 2)
                    return 2;
            }
        }

        return ret;
    }

    int coefficient_type() const noexcept
    {
        int eq = coefficient_type(equal_constraints);
        int le = coefficient_type(less_constraints);
        int ge = coefficient_type(greater_constraints);

        return std::max(eq, std::max(le, ge));
    }

    problem_solver_type which_problem_type() const noexcept
    {
        auto coefficient = coefficient_type();

        if (greater_constraints.empty() && less_constraints.empty()) {
            switch (coefficient) {
            case 0:
                return problem_solver_type::equalities_01;
            case 1:
                return problem_solver_type::equalities_101;
            default:
                return problem_solver_type::equalities_Z;
            }
        } else {
            switch (coefficient) {
            case 0:
                return problem_solver_type::inequalities_01;
            case 1:
                return problem_solver_type::inequalities_101;
            default:
                return problem_solver_type::inequalities_Z;
            }
        }

        return problem_solver_type::inequalities_Z;
    }

    string_buffer_ptr strings;

    objective_function objective;

    std::vector<constraint> equal_constraints;
    std::vector<constraint> greater_constraints;
    std::vector<constraint> less_constraints;

    variables vars;
    affected_variables affected_vars;

    objective_function_type type = { objective_function_type::maximize };
    problem_solver_type problem_type = { problem_solver_type::equalities_01 };

    void swap(problem& other)
    {
        std::swap(other.objective, objective);
        std::swap(other.equal_constraints, equal_constraints);
        std::swap(other.greater_constraints, greater_constraints);
        std::swap(other.less_constraints, less_constraints);
        std::swap(other.vars, vars);
        std::swap(other.affected_vars, affected_vars);
        std::swap(other.type, type);
        std::swap(other.problem_type, problem_type);
    }
};

inline bool
operator==(const problem& lhs, const problem& rhs) noexcept
{
    return lhs.objective == rhs.objective &&
           lhs.equal_constraints == rhs.equal_constraints &&
           lhs.greater_constraints == rhs.greater_constraints &&
           lhs.less_constraints == rhs.less_constraints &&
           lhs.vars == rhs.vars && lhs.affected_vars == rhs.affected_vars &&
           lhs.type == rhs.type && lhs.problem_type == rhs.problem_type;
}

inline bool
operator!=(const problem& lhs, const problem& rhs) noexcept
{
    return !(lhs == rhs);
}

void
clear(raw_problem& pb);

bool
is_valid_solution(const problem& pb,
                  const std::vector<std::int8_t>& variable_value);

double
compute_solution(const problem& pb,
                 const std::vector<std::int8_t>& variable_value);

problem
unpreprocess(const context& ctx, const raw_problem& pb_);

problem
preprocess(const context& ctx, const raw_problem& pb_);

raw_problem
make_problem(const context& ctx, std::istream& is) noexcept;

/**
 * @brief Affect a variable to the @c pb problem.
 * @details Build a new @c problem from an original
 *  @c problem but a variable is assigned a specific value. All
 *      constraints, variables, objective function && affected variables are
 *      recursively updated.
 *
 * @param ctx Current context
 * @param pb Original problem
 * @param variable_index Which variable to assign.
 * @param variable_value The value of the binary variable.
 * @return A newly @c problem.
 *
 * @note Implemented in the private @c preprocessor-2.cpp file.
 */
problem
affect(const context& ctx,
       const problem& pb,
       int variable_index,
       bool variable_value);

/**
 * @brief Split the @c pb problem in two new problem.
 * @details Build two @c problem from an original
 *  @c problem but a variable is assigned to respectively 0 or 1. All
 *      constraints, variables, objective function && affected variables are
 *      recursively updated.
 *
 * @param ctx Current context
 * @param pb Original problem
 * @param variable_index_to_affect Which variable to assign respectively to 0
 *     or 1.
 * @return Two newly @c problem.
 *
 * @note Implemented in the private @c preprocessor-2.cpp file.
 */
std::tuple<problem, problem>
split(const context& ctx, const problem& pb, int variable_index_to_affect);

} // namespace baryonyx

template<typename FormatContext, typename Problem, typename Function>
void
write_function_element(FormatContext& ctx, const Problem& p, const Function& f)
{
    for (auto& elem : f) {
        if (elem.factor < 0)
            fmt::format_to(ctx.out(),
                           " {} {}",
                           elem.factor,
                           p.vars.names[elem.variable_index]);
        else if (elem.factor > 0)
            fmt::format_to(ctx.out(),
                           "+ {} {}",
                           elem.factor,
                           p.vars.names[elem.variable_index]);

        fmt::format_to(ctx.out(), " ");
    }
}

template<typename FormatContext, typename Problem, typename Function>
void
write_quadratic_element(FormatContext& ctx,
                        const Problem& p,
                        const Function& f)
{
    if (f.empty())
        return;

    fmt::format_to(ctx.out(), "+ [");

    for (auto& elem : f) {
        if (elem.variable_index_a == elem.variable_index_b)
            fmt::format_to(ctx.out(),
                           " {} {} ^2",
                           2.0 * elem.factor,
                           p.vars.names[elem.variable_index_a]);
        else
            fmt::format_to(ctx.out(),
                           " {} {} * {}",
                           2.0 * elem.factor,
                           p.vars.names[elem.variable_index_a],
                           p.vars.names[elem.variable_index_b]);
    }

    fmt::format_to(ctx.out(), " ] /2");
}

template<typename FormatContext, typename Problem, typename Constraint>
void
write_constraint(FormatContext& ctx,
                 const Problem& p,
                 const Constraint& cst,
                 const char* separator)
{
    if (!cst.label.empty())
        fmt::format_to(ctx.out(), "{}: ", cst.label);

    write_function_element(ctx, p, cst.elements);

    fmt::format_to(ctx.out(), "{}{}\n", separator, cst.value);
}

template<typename FormatContext, typename Problem>
void
write_constraints(FormatContext& ctx, const Problem& pb)
{
    for (std::size_t i = 0, e = pb.equal_constraints.size(); i != e; ++i)
        write_constraint(ctx, pb, pb.equal_constraints[i], " = ");

    for (std::size_t i = 0, e = pb.greater_constraints.size(); i != e; ++i)
        write_constraint(ctx, pb, pb.greater_constraints[i], " >= ");

    for (std::size_t i = 0, e = pb.less_constraints.size(); i != e; ++i)
        write_constraint(ctx, pb, pb.less_constraints[i], " <= ");
}

template<typename FormatContext, typename Problem>
void
write_bounds(FormatContext& ctx, const Problem& p)
{
    for (std::size_t i{ 0 }, e{ p.vars.names.size() }; i != e; ++i) {
        if (p.vars.values[i].min != 0)
            fmt::format_to(
              ctx.out(), "{} >= {}\n", p.vars.names[i], p.vars.values[i].min);

        if (p.vars.values[i].max != std::numeric_limits<int>::max())
            fmt::format_to(
              ctx.out(), "{} <= {}\n", p.vars.names[i], p.vars.values[i].max);
    }
}

template<typename FormatContext, typename Problem>
constexpr void
write_problem(FormatContext& ctx, const Problem& p)
{
    if (p.vars.names.empty())
        return;

    fmt::format_to(ctx.out(),
                   "{}\n",
                   p.type == baryonyx::objective_function_type::maximize
                     ? "maximize\n"
                     : "minimize\n");

    write_function_element(ctx, p, p.objective.elements);
    write_quadratic_element(ctx, p, p.objective.qelements);

    if (p.objective.value < 0)
        fmt::format_to(ctx.out(), "{}", p.objective.value);
    else if (p.objective.value > 0)
        fmt::format_to(ctx.out(), " + {}", p.objective.value);

    fmt::format_to(ctx.out(), "\nsubject to\n");
    write_constraints(ctx, p);

    fmt::format_to(ctx.out(), "bounds\n");
    write_bounds(ctx, p);

    bool have_binary = false;
    bool have_general = false;

    for (std::size_t i{ 0 }, e{ p.vars.names.size() }; i != e; ++i) {
        if (p.vars.values[i].type == baryonyx::variable_type::binary) {
            have_binary = true;
            if (have_general == true)
                break;
        } else if (p.vars.values[i].type == baryonyx::variable_type::general) {
            have_general = true;
            if (have_binary == true)
                break;
        }
    }

    if (have_binary) {
        fmt::format_to(ctx.out(), "binary\n");
        for (std::size_t i{ 0 }, e{ p.vars.names.size() }; i != e; ++i)
            if (p.vars.values[i].type == baryonyx::variable_type::binary)
                fmt::format_to(ctx.out(), " {}\n", p.vars.names[i]);
    }

    if (have_general) {
        fmt::format_to(ctx.out(), "general\n");
        for (std::size_t i{ 0 }, e{ p.vars.names.size() }; i != e; ++i)
            if (p.vars.values[i].type == baryonyx::variable_type::general)
                fmt::format_to(ctx.out(), " {}\n", p.vars.names[i]);
    }

    fmt::format_to(ctx.out(), "end\n");
};

template<>
struct fmt::formatter<baryonyx::problem>
{
    constexpr auto parse(format_parse_context& ctx)
      -> format_parse_context::iterator
    {
        return ctx.begin();
    }

    constexpr auto format(const baryonyx::problem& pb,
                          format_context& ctx) const
      -> format_context::iterator
    {
        write_problem(ctx, pb);
        return ctx.out();
    }
};

template<>
struct fmt::formatter<baryonyx::raw_problem>
{
    constexpr auto parse(format_parse_context& ctx)
    {
        return ctx.begin();
    }

    constexpr auto format(const baryonyx::raw_problem& pb,
                          format_context& ctx) const
      -> format_context::iterator
    {
        write_problem(ctx, pb);
        return ctx.out();
    }
};

#endif
