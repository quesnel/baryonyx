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

#ifndef ORG_VLEPROJECT_BARYONYX_LIB_PROBLEM_HPP
#define ORG_VLEPROJECT_BARYONYX_LIB_PROBLEM_HPP

#include <baryonyx/core-compare>
#include <baryonyx/core>

namespace baryonyx {

/**
 * @brief Internal representation of a problem with preprocessed data.
 *
 * @details A equivalent class to @c problem but with preprocessed
 * structures to store affected variables and duplicated constraints from the
 * originial problem.
 *
 */
struct problem
{
    problem() = default;

    problem(const raw_problem& pb)
      : objective(pb.objective)
      , equal_constraints(pb.equal_constraints)
      , greater_constraints(pb.greater_constraints)
      , less_constraints(pb.less_constraints)
      , vars(pb.vars)
      , type(pb.type)
    {}

    problem(raw_problem&& pb)
      : objective(pb.objective)
      , equal_constraints(pb.equal_constraints)
      , greater_constraints(pb.greater_constraints)
      , less_constraints(pb.less_constraints)
      , vars(pb.vars)
      , type(pb.type)
      , problem_type(which_problem_type())
    {}

    template<typename Constraints>
    int coefficient_type(const Constraints& csts) noexcept
    {
        int ret = 0;

        for (const auto& cst : csts) {
            for (const auto& elem : cst.elements) {
                if (elem.factor == -1 and ret == 0)
                    ret = 1;
                else if (elem.factor <= -2 or elem.factor >= 2)
                    return 2;
            }
        }

        return ret;
    }

    int coefficient_type() noexcept
    {
        int ret = 0;

        ret = coefficient_type(equal_constraints);
        if (ret != 2) {
            ret = coefficient_type(less_constraints);
            if (ret != 2) {
                ret = coefficient_type(less_constraints);
            }
        }

        return ret;
    }

    problem_solver_type which_problem_type() noexcept
    {
        auto coefficient = coefficient_type();

        if (greater_constraints.empty() and less_constraints.empty()) {
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

    objective_function objective;

    std::vector<constraint> equal_constraints;
    std::vector<constraint> greater_constraints;
    std::vector<constraint> less_constraints;

    variables vars;
    affected_variables affected_vars;

    objective_function_type type = { objective_function_type::maximize };
    problem_solver_type problem_type = { problem_solver_type::equalities_01 };
};

inline bool
operator==(const problem& lhs, const problem& rhs) noexcept
{
    return lhs.objective == rhs.objective and
           lhs.equal_constraints == rhs.equal_constraints and
           lhs.greater_constraints == rhs.greater_constraints and
           lhs.less_constraints == rhs.less_constraints and
           lhs.vars == rhs.vars and lhs.affected_vars == rhs.affected_vars and
           lhs.type == rhs.type and lhs.problem_type == rhs.problem_type;
}

inline bool
operator!=(const problem& lhs, const problem& rhs) noexcept
{
    return not(lhs == rhs);
}

bool
is_valid_solution(const problem& pb, const std::vector<bool>& variable_value);

double
compute_solution(const problem& pb, const std::vector<bool>& variable_value);

// bool
// is_valid_solution(const problem& pb, const result& r);

// double
// compute_solution(const problem& pb, const result& r);

problem
unpreprocess(const context_ptr& ctx, const raw_problem& pb_);

problem
preprocess(const context_ptr& ctx, const raw_problem& pb_);

raw_problem
read_problem(std::istream& is);

result
read_result(std::istream& is);

void
check_consistency(const raw_problem& pb);

void
check_consistency(const problem& pb);

/**
 * @brief Select the correct solver (only itm is available) and solve the
 *     problem.
 *
 * @param ctx Current context.
 * @param pb Preprocessed problem.
 *
 * @note Implemented in the private @c select.cpp file.
 */
result
solver_select(const context_ptr& ctx, const problem& pb);

/**
 * @brief Select the correct optimizez (only itm is available) and solve the
 *     problem.
 *
 * @param ctx Current context.
 * @param pb Preprocessed problem.
 *
 * @note Implemented in the private @c select.cpp file.
 */
result
optimizer_select(const context_ptr& ctx, const problem& pb);

/**
 * @brief Affect a variable to the @c pb problem.
 * @details Build a new @c problem from an original
 *  @c problem but a variable is assigned a specific value. All
 *      constraints, variables, objective function and affected variables are
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
affect(const context_ptr& ctx,
       const problem& pb,
       int variable_index,
       bool variable_value);

/**
 * @brief Split the @c pb problem in two new problem.
 * @details Build two @c problem from an original
 *  @c problem but a variable is assigned to respectively 0 or 1. All
 *      constraints, variables, objective function and affected variables are
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
split(const context_ptr& ctx, const problem& pb, int variable_index_to_affect);

} // namespace baryonyx

#endif
