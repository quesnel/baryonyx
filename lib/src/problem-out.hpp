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

#ifndef ORG_VLEPROJECT_BARYONYX_LIB_PROBLEM_OUT_HPP
#define ORG_VLEPROJECT_BARYONYX_LIB_PROBLEM_OUT_HPP

#include <array>
#include <chrono>
#include <iomanip>
#include <iterator>
#include <numeric>
#include <ostream>

#include <baryonyx/core>

#include "problem.hpp"

namespace baryonyx {

namespace detail {

template<typename Problem, typename Function>
void
write_function_element(std::ostream& os, const Problem& p, const Function& f)
{
    for (auto& elem : f) {
        os << ((elem.factor < 0) ? "- " : "+ ");
        if (elem.factor != 1)
            os << std::abs(elem.factor) << ' ';

        os << p.vars.names[elem.variable_index] << ' ';
    }
}

template<typename Problem, typename Constraint>
void
write_constraint(std::ostream& os,
                 const Problem& p,
                 const Constraint& cst,
                 const char* separator)
{
    if (!cst.label.empty())
        os << cst.label << ": ";

    write_function_element(os, p, cst.elements);
    os << separator << cst.value << '\n';
}

template<typename Problem>
void
write_constraints(std::ostream& os, const Problem& pb)
{
    for (std::size_t i = 0, e = pb.equal_constraints.size(); i != e; ++i)
        write_constraint(os, pb, pb.equal_constraints[i], " = ");

    for (std::size_t i = 0, e = pb.greater_constraints.size(); i != e; ++i)
        write_constraint(os, pb, pb.greater_constraints[i], " >= ");

    for (std::size_t i = 0, e = pb.less_constraints.size(); i != e; ++i)
        write_constraint(os, pb, pb.less_constraints[i], " <= ");
}

template<typename Problem>
void
write_bounds(std::ostream& os, const Problem& p)
{
    for (std::size_t i{ 0 }, e{ p.vars.names.size() }; i != e; ++i) {
        if (p.vars.values[i].min != 0)
            os << p.vars.names[i] << " >= " << p.vars.values[i].min << '\n';

        if (p.vars.values[i].max != std::numeric_limits<int>::max())
            os << p.vars.names[i] << " <= " << p.vars.values[i].max << '\n';
    }
}

template<typename Problem>
void
write_problem(std::ostream& os, const Problem& p)
{
    if (p.vars.names.empty())
        return;

    if (p.type == objective_function_type::maximize)
        os << "maximize\n";
    else
        os << "minimize\n";

    detail::write_function_element(os, p, p.objective.elements);

    if (p.objective.value < 0)
        os << p.objective.value;
    else if (p.objective.value > 0)
        os << " + " << p.objective.value;

    os << "\nsubject to\n";
    detail::write_constraints(os, p);

    os << "bounds\n";
    detail::write_bounds(os, p);

    bool have_binary = false;
    bool have_general = false;

    for (std::size_t i{ 0 }, e{ p.vars.names.size() }; i != e; ++i) {
        if (p.vars.values[i].type == variable_type::binary) {
            have_binary = true;
            if (have_general == true)
                break;
        } else if (p.vars.values[i].type == variable_type::general) {
            have_general = true;
            if (have_binary == true)
                break;
        }
    }

    if (have_binary) {
        os << "binary\n";
        for (std::size_t i{ 0 }, e{ p.vars.names.size() }; i != e; ++i)
            if (p.vars.values[i].type == variable_type::binary)
                os << ' ' << p.vars.names[i] << '\n';
    }

    if (have_general) {
        os << "general\n";
        for (std::size_t i{ 0 }, e{ p.vars.names.size() }; i != e; ++i)
            if (p.vars.values[i].type == variable_type::general)
                os << ' ' << p.vars.names[i] << '\n';
    }

    os << "end\n";
}

} // namespace detail

/**
 * Write @e lp problem into a stream.
 *
 */
inline std::ostream&
operator<<(std::ostream& os, const problem& p)
{
    if (os)
        detail::write_problem(os, p);

    return os;
}

/**
 * Write @e lp problem into a stream.
 *
 */
inline std::ostream&
operator<<(std::ostream& os, const raw_problem& p)
{
    if (os)
        detail::write_problem(os, p);

    return os;
}

/**
 * @brief Write baryonyx::result into a `dot sol` format.
 * @details The `dot sol` format use the same comment and variable names as the
 *     `dot lp` format. All meta data, constraint remaining, duration are
 *     stored into comment in `dot sol` format. Only variable affectation are
 *     full useful.
 *
 * @param os [out] class output std::ostream.
 * @param result [in] the baryonyx::result to write.
 */
inline std::ostream&
operator<<(std::ostream& os, const result& result)
{
    auto store = os.flags();

    os.flags(store | os.boolalpha);
    os << R"(\ solver................: )" << result.method << '\n'
       << R"(\ constraints...........: )" << result.constraints << '\n'
       << R"(\ variables.............: )" << result.variables << '\n'
       << R"(\ duration..............: )" << result.duration << "s\n"
       << R"(\ loop..................: )" << result.loop << '\n'
       << R"(\ status................: )";

    switch (result.status) {
    case result_status::internal_error:
        os << "internal error reported\n";
        break;
    case result_status::uninitialized:
        os << "uninitialized\n";
        break;
    case result_status::success:
        os << "solution found\n";

        if (result.solutions.empty()) // Baryonyx ensures solutions are not
            break;                    // empty.

        os << R"(\ value.................: )" << result.solutions.back().value
           << '\n'
           << R"(\ other value...........: )";

        for (const auto& elem : result.solutions)
            os << elem.value << ' ';
        os << '\n';

        os << "\\ variables.............: \n";

        for (std::size_t i{ 0 }, e{ result.affected_vars.names.size() };
             i != e;
             ++i)
            os << result.affected_vars.names[i] << '='
               << (result.affected_vars.values[i] ? 1 : 0) << '\n';

        for (std::size_t i{ 0 }, e{ result.variable_name.size() }; i != e; ++i)
            os << result.variable_name[i] << '='
               << (result.solutions.back().variables[i] ? 1 : 0) << '\n';
        break;
    case result_status::time_limit_reached:
        os << "time limit reached\n"
           << R"(\ remaining constraints.: )" << result.remaining_constraints
           << '\n';
        break;
    case result_status::kappa_max_reached:
        os << "kappa max reached\n"
           << R"(\ remaining constraints.: )" << result.remaining_constraints
           << '\n';
        break;
    case result_status::limit_reached:
        os << "limit reached\n"
           << R"(\ remaining constraints.: )" << result.remaining_constraints
           << '\n';
        break;
    }

    os.flags(store);

    return os;
}

struct resume
{
    template<typename Problem>
    resume(const Problem& pb, bool use_lp_format_ = true)
      : variables({})
      , constraints({})
      , minmax(compute_min_max_objective_function(pb.objective))
      , use_lp_format(use_lp_format_)
    {
        variables = std::accumulate(pb.vars.values.begin(),
                                    pb.vars.values.end(),
                                    variables,
                                    [](std::array<int, 3>& value, auto vv) {
                                        switch (vv.type) {
                                        case variable_type::real:
                                            value[0]++;
                                            break;
                                        case variable_type::binary:
                                            value[1]++;
                                            break;
                                        case variable_type::general:
                                            value[2]++;
                                            break;
                                        }

                                        return value;
                                    });

        constraints[0] = static_cast<int>(pb.equal_constraints.size());
        constraints[1] = static_cast<int>(pb.greater_constraints.size());
        constraints[2] = static_cast<int>(pb.less_constraints.size());

        problem_type = get_problem_type(pb);
    }

    std::string get_problem_type(const problem& pb) const
    {
        switch (pb.problem_type) {
        case baryonyx::problem_solver_type::equalities_01:
            return "equalities-01";
        case baryonyx::problem_solver_type::equalities_101:
            return "equalities-101";
        case baryonyx::problem_solver_type::equalities_Z:
            return "equalities-Z";
        case baryonyx::problem_solver_type::inequalities_01:
            return "inequalities-01";
        case baryonyx::problem_solver_type::inequalities_101:
            return "inequalities-101";
        case baryonyx::problem_solver_type::inequalities_Z:
            return "inequalities-Z";
        default:
            return {};
        }
    }

    std::string get_problem_type(const raw_problem&) const
    {
        return {};
    }

    std::array<int, 3> variables;
    std::array<int, 3> constraints;
    std::tuple<double, double> minmax;
    std::string problem_type;
    bool use_lp_format;
};

inline std::ostream&
operator<<(std::ostream& os, const resume& pb)
{
    auto store = os.flags();

    os << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    if (pb.use_lp_format) {
        os << "\\ Problem statistics:\n"
           << R"(\  type: )" << pb.problem_type << '\n'
           << R"(\  nb variables: )"
           << std::accumulate(pb.variables.begin(), pb.variables.end(), 0)
           << '\n'
           << R"(\   ..... real: )" << pb.variables[0] << '\n'
           << R"(\   ... binary: )" << pb.variables[1] << '\n'
           << R"(\   .. general: )" << pb.variables[2] << '\n'
           << R"(\  nb constraints: )"
           << std::accumulate(pb.constraints.begin(), pb.constraints.end(), 0)
           << '\n'
           << R"(\   ........ =  : )" << pb.constraints[0] << '\n'
           << R"(\   ........ >= : )" << pb.constraints[1] << '\n'
           << R"(\   ........ <= : )" << pb.constraints[2] << '\n'
           << R"(\  minimal value.: )" << std::get<0>(pb.minmax) << '\n'
           << R"(\  maximal value.: )" << std::get<1>(pb.minmax) << '\n';
    } else {
        os << "Problem statistics:\n"
           << "  * type: " << pb.problem_type << '\n'
           << "  * variables: "
           << std::accumulate(pb.variables.begin(), pb.variables.end(), 0)
           << '\n'
           << "    - real: " << pb.variables[0] << '\n'
           << "    - binary: " << pb.variables[1] << '\n'
           << "    - general: " << pb.variables[2] << '\n'
           << "  * constraints: "
           << std::accumulate(pb.constraints.begin(), pb.constraints.end(), 0)
           << '\n'
           << "    - constraint =  : " << pb.constraints[0] << '\n'
           << "    - constraint >= : " << pb.constraints[1] << '\n'
           << "    - constraint <= : " << pb.constraints[2] << '\n'
           << "  * objective:\n"
           << "    - minimal value.: " << std::get<0>(pb.minmax) << '\n'
           << "    - maximal value.: " << std::get<1>(pb.minmax) << '\n';
    }

    os.flags(store);

    return os;
}

inline std::ostream&
operator<<(std::ostream& os, const affected_variables& var)
{
    std::size_t i = 0, e = var.names.size();

    for (; i != e; ++i)
        os << var.names[i] << ": " << (var.values[i] ? 1 : 0) << '\n';

    return os;
}

struct best_solution_writer
{
    const baryonyx::result& res;

    best_solution_writer(const baryonyx::result& res_)
      : res(res_)
    {}
};

inline std::ostream&
operator<<(std::ostream& os, const best_solution_writer& writer)
{
    std::size_t i = 0, e = writer.res.variable_name.size();

    for (; i != e; ++i)
        os << writer.res.variable_name[i] << ": "
           << (writer.res.solutions.back().variables[i] ? 1 : 0) << '\n';

    return os;
}
} // namespace baryonyx

#endif
