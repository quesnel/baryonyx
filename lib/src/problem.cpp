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

#include "problem.hpp"
#include "utils.hpp"

#include <deque>
#include <istream>
#include <limits>
#include <ostream>
#include <unordered_map>

#include <fmt/format.h>

namespace {

template<typename Problem, typename Function>
void
write_function_element(std::ostream& os, const Problem& p, const Function& f)
{
    for (auto& elem : f) {
        if (elem.factor < 0)
            os << "- ";
        else
            os << "+ ";

        if (elem.factor != 1)
            os << std::abs(elem.factor) << ' ';

        os << p.vars.names[elem.variable_index] << ' ';
    }
}

template<typename Problem, typename Function>
void
write_quadratic_element(std::ostream& os, const Problem& p, const Function& f)
{
    if (f.empty())
        return;

    os << "+ [";

    for (auto& elem : f) {
        if (elem.variable_index_a == elem.variable_index_b)
            os << ' ' << (2.0 * elem.factor) << ' '
               << p.vars.names[elem.variable_index_a] << " ^2";
        else
            os << ' ' << (2.0 * elem.factor) << ' '
               << p.vars.names[elem.variable_index_a] << " * "
               << p.vars.names[elem.variable_index_b];
    }

    os << " ] /2";
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

    if (p.type == baryonyx::objective_function_type::maximize)
        os << "maximize\n";
    else
        os << "minimize\n";

    ::write_function_element(os, p, p.objective.elements);
    ::write_quadratic_element(os, p, p.objective.qelements);

    if (p.objective.value < 0)
        os << p.objective.value;
    else if (p.objective.value > 0)
        os << " + " << p.objective.value;

    os << "\nsubject to\n";
    ::write_constraints(os, p);

    os << "bounds\n";
    ::write_bounds(os, p);

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
        os << "binary\n";
        for (std::size_t i{ 0 }, e{ p.vars.names.size() }; i != e; ++i)
            if (p.vars.values[i].type == baryonyx::variable_type::binary)
                os << ' ' << p.vars.names[i] << '\n';
    }

    if (have_general) {
        os << "general\n";
        for (std::size_t i{ 0 }, e{ p.vars.names.size() }; i != e; ++i)
            if (p.vars.values[i].type == baryonyx::variable_type::general)
                os << ' ' << p.vars.names[i] << '\n';
    }

    os << "end\n";
}

} // anonymous namespace

namespace baryonyx {

/**
 * Write @e lp problem into a stream.
 *
 */
std::ostream&
operator<<(std::ostream& os, const problem& p)
{
    if (os)
        ::write_problem(os, p);

    return os;
}

/**
 * Write @e lp problem into a stream.
 *
 */
std::ostream&
operator<<(std::ostream& os, const raw_problem& p)
{
    if (os)
        ::write_problem(os, p);

    return os;
}

void
clear(raw_problem& pb)
{
    std::vector<objective_function_element>().swap(pb.objective.elements);
    std::vector<objective_quadratic_element>().swap(pb.objective.qelements);

    std::vector<constraint>().swap(pb.equal_constraints);
    std::vector<constraint>().swap(pb.greater_constraints);
    std::vector<constraint>().swap(pb.less_constraints);

    std::vector<std::string_view>().swap(pb.vars.names);
    std::vector<variable_value>().swap(pb.vars.values);

    pb.type = objective_function_type::maximize;
}

} // baryonyx namespace
