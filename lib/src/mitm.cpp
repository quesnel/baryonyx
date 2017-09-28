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

#include "generalized-wedelin.hpp"
#include "inequalities-1coeff.hpp"
#include "lpformat-io.hpp"
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
inline int
get_thread_number(std::shared_ptr<baryonyx::context>& ctx) noexcept
{
    auto t = ctx->get_integer_parameter("thread",
                                        std::thread::hardware_concurrency());

    return t <= 0 ? 1 : baryonyx::numeric_cast<int>(t);
}

std::tuple<double, double, double, int>
get_parameters(std::shared_ptr<baryonyx::context>& ctx) noexcept
{
    auto kappa = ctx->get_real_parameter("kappa", 0.001);
    auto theta = ctx->get_real_parameter("theta", 0.001);
    auto delta = ctx->get_real_parameter("delta", 0.001);
    auto limit = ctx->get_real_parameter("limit", 1000l);

    return std::make_tuple(kappa, delta, theta, limit);
}

template<typename variableT>
bool
is_boolean_variable(const variableT& vars)
{
    for (const auto& elem : vars)
        if (elem.type != variable_type::binary)
            return false;

    return true;
}

template<typename constraintsT>
bool
is_boolean_coefficient(const constraintsT csts)
{
    for (const auto& cst : csts)
        for (const auto& elem : cst.elements)
            if (elem.factor < 0 or elem.factor > 1)
                return false;

    return true;
}

template<typename variableT>
bool
is_integer_variable(const variableT& vars)
{
    for (const auto& elem : vars)
        if (elem.type != variable_type::general)
            return false;

    return true;
}

template<typename constraintsT>
bool
is_101_coefficient(const constraintsT csts)
{
    for (const auto& cst : csts)
        for (const auto& elem : cst.elements)
            if (elem.factor < -1 or elem.factor > 1)
                return false;

    return true;
}

template<typename functionT>
functionT
cleanup_function_element(functionT& fct, int& nb)
{
    std::size_t fct_size{ fct.size() };

    {
        auto it = std::remove_if(fct.begin(), fct.end(), [](const auto& elem) {
            return elem.factor == 0;
        });
        nb += std::distance(it, fct.end());
        fct.erase(it, fct.end());
    }

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
        nb += numeric_cast<int>(fct_size - ret.size());
    }

    return ret;
}

int
remove_variable(objective_function& obj, int variable_id, int variable_value)
{
    auto it = std::remove_if(obj.elements.begin(),
                             obj.elements.end(),
                             [variable_id](const auto& elem) {
                                 return elem.variable_index == variable_id;
                             });

    int ret = std::distance(it, obj.elements.end());

    for (auto jt = it; jt != obj.elements.end(); ++jt)
        obj.value += (jt->factor * variable_value);

    obj.elements.erase(it, obj.elements.end());

    for (auto& elem : obj.elements)
        if (elem.variable_index > variable_id)
            --elem.variable_index;

    return ret;
}

int
remove_variable(constraint& obj, int variable_id, int variable_value)
{
    auto it = std::remove_if(obj.elements.begin(),
                             obj.elements.end(),
                             [variable_id](const auto& elem) {
                                 return elem.variable_index == variable_id;
                             });

    int ret = std::distance(it, obj.elements.end());

    for (auto jt = it; jt != obj.elements.end(); ++jt)
        obj.value -= (jt->factor * variable_value);

    obj.elements.erase(it, obj.elements.end());

    for (auto& elem : obj.elements)
        if (elem.variable_index > variable_id)
            --elem.variable_index;

    return ret;
}

void
remove_variable(problem& pb,
                int variable_id,
                int variable_value,
                int& clean_nb,
                int& constraint_nb)
{
    pb.affected_vars.push_back(pb.vars.names[variable_id], variable_value);

    pb.vars.names.erase(pb.vars.names.begin() + variable_id);
    pb.vars.values.erase(pb.vars.values.begin() + variable_id);

    clean_nb += remove_variable(pb.objective, variable_id, variable_value);

    for (auto& elem : pb.equal_constraints)
        clean_nb += remove_variable(elem, variable_id, variable_value);

    for (auto& elem : pb.less_constraints)
        clean_nb += remove_variable(elem, variable_id, variable_value);

    for (auto& elem : pb.greater_constraints)
        clean_nb += remove_variable(elem, variable_id, variable_value);

    constraint_nb +=
      (pb.equal_constraints.size() + pb.less_constraints.size() +
       pb.greater_constraints.size());

    pb.equal_constraints.erase(
      std::remove_if(pb.equal_constraints.begin(),
                     pb.equal_constraints.end(),
                     [](const auto& v) { return v.elements.empty(); }),
      pb.equal_constraints.end());

    pb.less_constraints.erase(
      std::remove_if(pb.less_constraints.begin(),
                     pb.less_constraints.end(),
                     [](const auto& v) { return v.elements.empty(); }),
      pb.less_constraints.end());

    pb.greater_constraints.erase(
      std::remove_if(pb.greater_constraints.begin(),
                     pb.greater_constraints.end(),
                     [](const auto& v) { return v.elements.empty(); }),
      pb.greater_constraints.end());

    constraint_nb -=
      (pb.equal_constraints.size() + pb.less_constraints.size() +
       pb.greater_constraints.size());
}

static inline size_t
hash_constraint(const std::vector<function_element>& e) noexcept
{
    std::size_t seed{ e.size() };

    for (auto& elem : e)
        seed ^= elem.variable_index + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    return seed;
}

struct hashed_constraints
{
    hashed_constraints(const std::vector<function_element>& elems, int index_)
      : hash(hash_constraint(elems))
      , index(index_)
    {
    }

    std::size_t hash;
    int index;
};

// Tries to remove duplicated constraints (i.e. all function elements in the
// vector must be equal) with or without the same value. To improve
// calculation, we use a mapping vector between an hash of the function
// elements vector and the index into the constraint vector. We merge
// duplicated constraints and build an index vector with constraints preserved.
void
remove_duplicated_constraints(std::shared_ptr<baryonyx::context> ctx,
                              std::vector<constraint>& cst,
                              operator_type type,
                              int& constraint_nb)
{
    std::vector<hashed_constraints> t;
    t.reserve(cst.size());

    std::vector<int> order;
    order.reserve(cst.size());

    for (int i{ 0 }, e{ numeric_cast<int>(cst.size()) }; i != e; ++i)
        t.emplace_back(cst[i].elements, i);

    std::sort(t.begin(), t.end(), [](const auto& rhs, const auto& lhs) {
        return rhs.hash < lhs.hash;
    });

    std::size_t i{ 0 }, e{ t.size() };

    while (i != e) {
        order.emplace_back(t[i].index);
        auto j = i + 1;

        while (j != e and t[i].hash == t[j].hash and
               cst[t[i].index].elements == cst[t[j].index].elements) {
            switch (type) {
            case operator_type::equal:
                Expects(cst[t[i].index].value == cst[t[j].index].value);

                ctx->debug("      = constraints %d = %d must be removed\n",
                           t[j].index,
                           t[i].index);
                break;
            case operator_type::greater:
                cst[t[i].index].value =
                  std::max(cst[t[i].index].value, cst[t[j].index].value);

                ctx->debug("      >= constraints %d = %d must be removed\n",
                           t[j].index,
                           t[i].index);
                break;
            case operator_type::less:
                cst[t[i].index].value =
                  std::min(cst[t[i].index].value, cst[t[j].index].value);

                ctx->debug("      <= constraints %d = %d must be removed\n",
                           t[j].index,
                           t[i].index);
                break;
            default:
                break;
            }

            ++j;
            ++constraint_nb;
        }

        i = j;
    }

    std::sort(order.begin(), order.end(), std::less<int>());

    std::vector<constraint> result(order.size());
    for (std::size_t i{ 0 }, e{ order.size() }; i != e; ++i)
        result[i] = std::move(cst[order[i]]);

    result.swap(cst);
}

static void
remove_duplicated_constraints(std::shared_ptr<baryonyx::context> ctx,
                              problem& pb,
                              int& constraint)
{
    remove_duplicated_constraints(
      ctx, pb.equal_constraints, operator_type::equal, constraint);
    remove_duplicated_constraints(
      ctx, pb.less_constraints, operator_type::less, constraint);
    remove_duplicated_constraints(
      ctx, pb.greater_constraints, operator_type::greater, constraint);
}

static void
remove_assigned_variables(std::shared_ptr<baryonyx::context> ctx,
                          problem& pb,
                          int& variable,
                          int& clean,
                          int& constraint)
{
    int i{ 0 }, e{ numeric_cast<int>(pb.vars.values.size()) };

    while (i != e) {
        if (pb.vars.values[i].min == pb.vars.values[i].max) {
            ctx->debug("      bound pass: variable %s = %d must be removed\n",
                       pb.vars.names[i].c_str(),
                       pb.vars.values[i].max);

            remove_variable(pb, i, pb.vars.values[i].max, clean, constraint);

            ++variable;
            i = 0;
            e = numeric_cast<int>(pb.vars.values.size());
        } else {
            ++i;
        }
    }

    bool modify = false;
    int pass{ 0 };

    do {
        modify = false;

        {
            i = 0;
            e = numeric_cast<int>(pb.equal_constraints.size());

            while (i != e) {
                if (pb.equal_constraints[i].elements.size() == 1) {
                    ctx->debug(
                      "      equal constraint (pass %d): variable %s = "
                      "%d must be removed\n",
                      pass,
                      pb.vars
                        .names[pb.equal_constraints[i]
                                 .elements.front()
                                 .variable_index]
                        .c_str(),
                      pb.equal_constraints[i].value);

                    remove_variable(
                      pb,
                      pb.equal_constraints[i].elements.front().variable_index,
                      pb.equal_constraints[i].value,
                      clean,
                      constraint);

                    modify = true;
                    ++variable;
                    i = 0;
                } else {
                    ++i;
                }
            }
        }

        ++pass;
    } while (modify);
}

void
cleanup_problem(std::shared_ptr<baryonyx::context> ctx, problem& pb)
{
    ctx->info("preprocessor:\n");

    {
        ctx->info("  - cleaning functions (merge variables, remove zero) "
                  "coefficient:\n");
        int clean{ 0 };

        for (auto& elem : pb.equal_constraints)
            elem.elements = cleanup_function_element(elem.elements, clean);
        for (auto& elem : pb.greater_constraints)
            elem.elements = cleanup_function_element(elem.elements, clean);
        for (auto& elem : pb.less_constraints)
            elem.elements = cleanup_function_element(elem.elements, clean);

        cleanup_function_element(pb.objective.elements, clean);

        ctx->info("    `-> %d function element(s) merged.\n", clean);
    }

    {
        ctx->info("  - cleaning already assigned variables:\n");
        int variable{ 0 }, clean{ 0 }, constraint{ 0 };

        remove_assigned_variables(ctx, pb, variable, clean, constraint);
        remove_duplicated_constraints(ctx, pb, constraint);

        ctx->info("    `-> %d variable(s) - %d function element(s) - %d "
                  "constraint(s) removed.\n",
                  variable,
                  clean,
                  constraint);
    }

#ifndef BARYONYX_FULL_OPTIMIZATION
    {
        ctx->info("  - write preprocessed.lp: ");
        std::ofstream ofs("preprocessed.lp");
        if (ofs.is_open()) {
            details::problem_writer pw(pb, ofs);
            ctx->info("writing %s\n", ((pw) ? "done" : "fail"));
        } else {
            ctx->info("opening failed\n");
        }
    }
#endif
}

result
mitm_solve(std::shared_ptr<baryonyx::context> ctx, problem& pb)
{
    cleanup_problem(ctx, pb);

    if (pb.greater_constraints.empty() and pb.less_constraints.empty() and
        is_boolean_coefficient(pb.equal_constraints) and
        is_boolean_variable(pb.vars.values))
        return baryonyx::inequalities_1coeff_wedelin_solve(ctx, pb);

    if (is_101_coefficient(pb.equal_constraints) and
        is_101_coefficient(pb.less_constraints) and
        is_101_coefficient(pb.greater_constraints) and
        is_boolean_variable(pb.vars.values))
        return baryonyx::inequalities_1coeff_wedelin_solve(ctx, pb);

    if ((is_101_coefficient(pb.equal_constraints) or
         is_101_coefficient(pb.greater_constraints) or
         is_101_coefficient(pb.less_constraints)) and
        is_integer_variable(pb.vars.values)) {
        double kappa, delta, theta;
        int limit;

        std::tie(kappa, delta, theta, limit) = get_parameters(ctx);

        return baryonyx::generalized_wedelin(kappa, delta, theta, limit, pb);
    }

    ctx->info("no_solver_available");

    throw baryonyx::solver_failure(solver_error_tag::no_solver_available);
}

result
mitm_optimize(std::shared_ptr<baryonyx::context> ctx, problem& pb)
{
    auto thread = get_thread_number(ctx);
    cleanup_problem(ctx, pb);

    if (pb.greater_constraints.empty() and pb.less_constraints.empty() and
        is_boolean_coefficient(pb.equal_constraints) and
        is_boolean_variable(pb.vars.values))
        return baryonyx::inequalities_1coeff_wedelin_optimize(ctx, pb, thread);

    if (is_101_coefficient(pb.equal_constraints) and
        is_101_coefficient(pb.less_constraints) and
        is_101_coefficient(pb.greater_constraints) and
        is_boolean_variable(pb.vars.values))
        return baryonyx::inequalities_1coeff_wedelin_optimize(ctx, pb, thread);

    if ((is_101_coefficient(pb.equal_constraints) or
         is_101_coefficient(pb.greater_constraints) or
         is_101_coefficient(pb.less_constraints)) and
        is_integer_variable(pb.vars.values)) {
        double kappa, delta, theta;
        int limit;

        std::tie(kappa, delta, theta, limit) = get_parameters(ctx);

        return baryonyx::generalized_wedelin(kappa, delta, theta, limit, pb);
    }

    ctx->info("no_solver_available");

    throw baryonyx::solver_failure(solver_error_tag::no_solver_available);
}

} // namespace baryonyx
