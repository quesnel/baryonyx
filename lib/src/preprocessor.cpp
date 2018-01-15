/* Copyright (C) 2017 INRA
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

#include "private.hpp"
#include "utils.hpp"

#include <fstream>
#include <unordered_set>

#include <cassert>

//
// Remove empty constraints, where the function does not have element from the
// problem.
//
static void
remove_empty_constraints(baryonyx::problem& pb, int& constraint_nb)
{
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

//
// Remove function element equal to 0, merge function elements with the same
// variable index.
//
// For example:
// x + 0 y + z + z will now equal to x + 2 z
//
// TODO take into account function with -x +x ... : x must be removed.
//
template<typename functionT>
static functionT
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
            assert(ret.back().variable_index == it->variable_index &&
                   "parser error to assign variable index.");

            ret.back().factor += it->factor;
        } else {
            prev = it;
            ret.emplace_back(*it);
        }
        ++it;
    }

    nb += baryonyx::numeric_cast<int>(fct_size - ret.size());

    return ret;
}

//
// Remove the variable from the objective function and add the `variable_value
// * factor` to the constant of the objective function.
//
static void
remove_variable_from_objective(baryonyx::objective_function& obj,
                               int variable_id,
                               int variable_value)
{
    auto it = std::find_if(obj.elements.begin(),
                           obj.elements.end(),
                           [variable_id](const auto& elem) {
                               return elem.variable_index == variable_id;
                           });

    if (it != obj.elements.end()) {
        obj.value += (it->factor * variable_value);
        obj.elements.erase(it);
    }

    assert(std::find_if(obj.elements.begin(),
                        obj.elements.end(),
                        [variable_id](const auto& elem) {
                            return elem.variable_index == variable_id;
                        }) == obj.elements.end());

    for (auto& elem : obj.elements)
        if (elem.variable_index > variable_id)
            --elem.variable_index;
}

//
// Remove the variable from the constraint and add the `factor * value`  to the
// value of the constraints.
//
static void
remove_variable(baryonyx::constraint& cst, int variable_id, int variable_value)
{
    auto it = std::find_if(cst.elements.begin(),
                           cst.elements.end(),
                           [variable_id](const auto& elem) {
                               return elem.variable_index == variable_id;
                           });

    if (it != cst.elements.end()) {
        cst.value -= (it->factor * variable_value);
        cst.elements.erase(it);
    }

    assert(std::find_if(cst.elements.begin(),
                        cst.elements.end(),
                        [variable_id](const auto& elem) {
                            return elem.variable_index == variable_id;
                        }) == cst.elements.end());

    for (auto& elem : cst.elements)
        if (elem.variable_index > variable_id)
            --elem.variable_index;
}

//
// Remove the specified variable to the problem and add the variable and its
// value to the problem affected variables vector.
//
static void
remove_variable(baryonyx::problem& pb, int variable_id, int variable_value)
{
    printf("remove variable %d value %d\n", variable_id, variable_value);
    pb.affected_vars.push_back(pb.vars.names[variable_id], variable_value);

    pb.vars.names.erase(pb.vars.names.begin() + variable_id);
    pb.vars.values.erase(pb.vars.values.begin() + variable_id);

    remove_variable_from_objective(pb.objective, variable_id, variable_value);

    for (auto& elem : pb.equal_constraints)
        remove_variable(elem, variable_id, variable_value);

    for (auto& elem : pb.less_constraints)
        remove_variable(elem, variable_id, variable_value);

    for (auto& elem : pb.greater_constraints)
        remove_variable(elem, variable_id, variable_value);
}

static inline size_t
hash_constraint(const std::vector<baryonyx::function_element>& e) noexcept
{
    std::size_t seed{ e.size() };

    for (auto& elem : e)
        seed ^= elem.variable_index + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    return seed;
}

struct hashed_constraints
{
    hashed_constraints(const std::vector<baryonyx::function_element>& elems,
                       int index_)
      : hash(hash_constraint(elems))
      , index(index_)
    {}

    std::size_t hash;
    int index;
};

//
// Tries to remove duplicated constraints (i.e. all function elements in the
// vector must be equal) with or without the same value. To improve
// calculation, we use a mapping vector between an hash of the function
// elements vector and the index into the constraint vector. We merge
// duplicated constraints and build an index vector with constraints preserved.
//
static void
remove_duplicated_constraints(const std::shared_ptr<baryonyx::context>& ctx,
                              std::vector<baryonyx::constraint>& cst,
                              baryonyx::operator_type type,
                              int& constraint_nb)
{
    std::vector<hashed_constraints> t;
    t.reserve(cst.size());

    std::vector<int> order;
    order.reserve(cst.size());

    for (int i{ 0 }, e{ baryonyx::numeric_cast<int>(cst.size()) }; i != e; ++i)
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
            case baryonyx::operator_type::equal:
                baryonyx::Expects(cst[t[i].index].value ==
                                  cst[t[j].index].value);

                ctx->debug("      = constraints %s = %s must be removed\n",
                           cst[t[j].index].label.c_str(),
                           cst[t[i].index].label.c_str());
                break;
            case baryonyx::operator_type::greater:
                cst[t[i].index].value =
                  std::max(cst[t[i].index].value, cst[t[j].index].value);

                ctx->debug("      >= constraints %s = %s must be removed\n",
                           cst[t[j].index].label.c_str(),
                           cst[t[i].index].label.c_str());
                break;
            case baryonyx::operator_type::less:
                cst[t[i].index].value =
                  std::min(cst[t[i].index].value, cst[t[j].index].value);

                ctx->debug("      <= constraints %s = %s must be removed\n",
                           cst[t[j].index].label.c_str(),
                           cst[t[i].index].label.c_str());
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

    std::vector<baryonyx::constraint> result(order.size());
    for (std::size_t i{ 0 }, e{ order.size() }; i != e; ++i)
        result[i] = std::move(cst[order[i]]);

    result.swap(cst);
}

static void
remove_duplicated_constraints(const std::shared_ptr<baryonyx::context>& ctx,
                              baryonyx::problem& pb,
                              int& constraint)
{
    remove_duplicated_constraints(
      ctx, pb.equal_constraints, baryonyx::operator_type::equal, constraint);
    remove_duplicated_constraints(
      ctx, pb.less_constraints, baryonyx::operator_type::less, constraint);
    remove_duplicated_constraints(ctx,
                                  pb.greater_constraints,
                                  baryonyx::operator_type::greater,
                                  constraint);
}

static std::tuple<int, int>
get_factor_type(const std::vector<baryonyx::function_element>& elems)
{
    int min = 0, max = 0;

    for (const auto& elem : elems)
        if (elem.factor < 0)
            min = std::min(min, elem.factor);
        else
            max = std::max(max, elem.factor);

    return std::make_tuple(min, max);
}

static bool
try_remove_assigned_variable(const std::shared_ptr<baryonyx::context>& ctx,
                             baryonyx::problem& pb,
                             const baryonyx::constraint& cst,
                             baryonyx::operator_type type,
                             int& variable_nb)
{
    assert(cst.elements.size() == 1);

    int factor = cst.elements.front().factor;
    int variable_index{ -1 }, variable_value{ -1 };

    if (factor == 1) {
        switch (type) {
        case baryonyx::operator_type::undefined:
            break;
        case baryonyx::operator_type::equal:
            // 1 x = 0 or 1 x = 1
            if (cst.value == 0 or cst.value == 1) {
                variable_index = cst.elements.front().variable_index;
                variable_value = cst.value;
            } else {
                assert(false and "1x = -1: error");
            }
            break;
        case baryonyx::operator_type::greater:
            // 1 x >= 1
            if (cst.value == 1) {
                variable_index = cst.elements.front().variable_index;
                variable_value = 1;
            } else
                return true; // 1x >= 0 or 1x >= -1
            break;
        case baryonyx::operator_type::less:
            // 1 x <= 0
            if (cst.value == 0) {
                variable_index = cst.elements.front().variable_index;
                variable_value = 0;
            } else if (cst.value == 1)
                return true;
            else
                assert(false and "1x <= -1");
            break;
        }
    } else {
        switch (type) {
        case baryonyx::operator_type::undefined:
            break;
        case baryonyx::operator_type::equal:
            if (cst.value == 0) { // -x = 0 -> x = 0
                variable_index = cst.elements.front().variable_index;
                variable_value = 0;
            } else if (cst.value == -1) { // -x = -1 -> x = 1
                variable_index = cst.elements.front().variable_index;
                variable_value = 1;
            } else {
                assert(false and "-1x = 1");
            }
            break;
        case baryonyx::operator_type::greater:
            if (cst.value == 0) { // -x >= 0 -> x = 0
                variable_index = cst.elements.front().variable_index;
                variable_value = 0;
            } else if (cst.value == -1) // -x >= -1 -> x = 1
                return true;
            else
                assert(false and "-1x >= 1");
            break;
        case baryonyx::operator_type::less:
            // -x <= -1 -> x = 1
            if (cst.value >= 0)
                return true;

            variable_index = cst.elements.front().variable_index;
            variable_value = 1;
            break;
        }
    }

    if (variable_index < 0)
        return false;

    ctx->debug("      variable %s = %d must be removed\n",
               pb.vars.names[variable_index].c_str(),
               variable_value);

    remove_variable(pb, variable_index, variable_value);

    ++variable_nb;

    return true;
}

static bool
try_remove_assigned_variables_01(const std::shared_ptr<baryonyx::context>& ctx,
                                 baryonyx::problem& pb,
                                 const baryonyx::constraint& cst,
                                 baryonyx::operator_type type,
                                 int& variable)
{
    auto nb = baryonyx::numeric_cast<int>(cst.elements.size());
    int value{ -1 };
    bool found{ false };

    switch (type) {
    case baryonyx::operator_type::undefined:
        break;
    case baryonyx::operator_type::equal:
        if (cst.value == 0 or cst.value == nb) {
            value = cst.value / nb;
            found = true;
        }
        break;
    case baryonyx::operator_type::greater:
        if (cst.value == nb) {
            value = 1;
            found = true;
        } else if (cst.value == 0)
            return true;
        break;
    case baryonyx::operator_type::less:
        if (cst.value == 0) {
            value = 0;
            found = true;
        } else if (cst.value == nb)
            return true;
        break;
    }

    if (not found)
        return false;

    std::vector<int> id;
    id.reserve(nb);

    for (int i = 0; i != nb; ++i)
        id.push_back(cst.elements[i].variable_index);

    while (not id.empty()) {
        ctx->debug("      variable %s = %d must be removed (2)\n",
                   pb.vars.names[id.back()].c_str(),
                   value);

        remove_variable(pb, id.back(), value);

        for (auto& elem : id)
            if (elem > id.back())
                --elem;

        id.pop_back();
        ++variable;
    }

    return true;
}

static bool
try_remove_assigned_variables_101(
  const std::shared_ptr<baryonyx::context>& ctx,
  baryonyx::problem& pb,
  const baryonyx::constraint& cst,
  baryonyx::operator_type type,
  int& variable)
{
    int min = 0, max = 0;
    for (const auto& elem : cst.elements)
        if (elem.factor < 0)
            --min;
        else
            ++max;

    std::vector<std::pair<int, int>> id;

    switch (type) {
    case baryonyx::operator_type::undefined:
        break;
    case baryonyx::operator_type::equal:
        if (cst.value == min) {
            for (const auto& elem : cst.elements)
                id.emplace_back(std::make_pair(elem.variable_index,
                                               (elem.factor < 0) ? 1 : 0));
        } else if (cst.value == max) {
            for (const auto& elem : cst.elements)
                id.emplace_back(std::make_pair(elem.variable_index,
                                               (elem.factor > 0) ? 1 : 0));
        }
        break;
    case baryonyx::operator_type::greater:
        if (cst.value == max) { // -x + y >= 1, x = 0 and y = 1
            for (const auto& elem : cst.elements)
                id.emplace_back(std::make_pair(elem.variable_index,
                                               (elem.factor > 0) ? 1 : 0));
        } else if (cst.value == min) // -x + y >= -1, useless constraint.
            return true;
        break;
    case baryonyx::operator_type::less:
        if (cst.value == min) { // -x + y <= -1, x = 1, y = 0
            for (const auto& elem : cst.elements)
                id.emplace_back(std::make_pair(elem.variable_index,
                                               (elem.factor < 0) ? 1 : 0));
        } else if (cst.value == max) // -x + y <= 1 useless constraint.
            return true;
        break;
    }

    if (id.empty())
        return false;

    while (not id.empty()) {
        ctx->debug("      variable %s = %d must be removed (2)\n",
                   pb.vars.names[id.back().first].c_str(),
                   id.back().second);

        remove_variable(pb, id.back().first, id.back().second);

        for (auto& elem : id)
            if (elem.first > id.back().first)
                --elem.first;

        id.pop_back();
        ++variable;
    }

    return true;
}

//
// Try to remove the current constraint by affecting variables. For example,
// if a constraint is defined as: x + y = 0, with x and y >= 0 then, variable x
// and y are equals to 0 and the constraint can be removed.
//
static bool
try_remove_assigned_variables(const std::shared_ptr<baryonyx::context>& ctx,
                              baryonyx::problem& pb,
                              const baryonyx::constraint& cst,
                              baryonyx::operator_type type,
                              int& variable)
{

    if (cst.elements.size() == 1)
        return try_remove_assigned_variable(ctx, pb, cst, type, variable);

    auto factor = get_factor_type(cst.elements);
    if (std::get<0>(factor) == 0 and std::get<1>(factor) == 1)
        return try_remove_assigned_variables_01(ctx, pb, cst, type, variable);
    else if (std::get<0>(factor) == -1 and std::get<1>(factor) == 1)
        return try_remove_assigned_variables_101(ctx, pb, cst, type, variable);

    return false;
}

static void
remove_assigned_variables(const std::shared_ptr<baryonyx::context>& ctx,
                          baryonyx::problem& pb,
                          int& variable,
                          int& constraint_nb)
{
    bool modify = false;

    int i = 0;
    auto e = baryonyx::numeric_cast<int>(pb.vars.values.size());
    while (i != e) {
        if (pb.vars.values[i].min == pb.vars.values[i].max) {
            ctx->debug("      bound: variable %s = %d must be removed\n",
                       pb.vars.names[i].c_str(),
                       pb.vars.values[i].max);

            remove_variable(pb, i, pb.vars.values[i].max);

            modify = true;
            ++variable;
            i = 0;
            e = baryonyx::numeric_cast<int>(pb.vars.values.size());
        } else {
            ++i;
        }
    }

    if (modify)
        remove_empty_constraints(pb, constraint_nb);

    do {
        modify = false;

        i = 0;
        e = baryonyx::numeric_cast<int>(pb.equal_constraints.size());
        while (i != e) {
            if (try_remove_assigned_variables(ctx,
                                              pb,
                                              pb.equal_constraints[i],
                                              baryonyx::operator_type::equal,
                                              variable)) {

                ctx->debug("      = constraint: %s must be removed\n",
                           pb.equal_constraints[i].label.c_str());
                pb.equal_constraints.erase(pb.equal_constraints.begin() + i);
                modify = true;
                ++constraint_nb;
                i = 0;
                e = baryonyx::numeric_cast<int>(pb.equal_constraints.size());
            } else {
                ++i;
            }
        }

        i = 0;
        e = baryonyx::numeric_cast<int>(pb.less_constraints.size());
        while (i != e) {
            if (try_remove_assigned_variables(ctx,
                                              pb,
                                              pb.less_constraints[i],
                                              baryonyx::operator_type::less,
                                              variable)) {

                ctx->debug("      <= constraint: %s must be removed\n",
                           pb.less_constraints[i].label.c_str());
                pb.less_constraints.erase(pb.less_constraints.begin() + i);
                modify = true;
                ++constraint_nb;
                i = 0;
                e = baryonyx::numeric_cast<int>(pb.less_constraints.size());
            } else {
                ++i;
            }
        }

        i = 0;
        e = baryonyx::numeric_cast<int>(pb.greater_constraints.size());
        while (i != e) {
            if (try_remove_assigned_variables(ctx,
                                              pb,
                                              pb.greater_constraints[i],
                                              baryonyx::operator_type::greater,
                                              variable)) {

                ctx->debug("      >= constraint: %s must be removed\n",
                           pb.greater_constraints[i].label.c_str());
                pb.greater_constraints.erase(pb.greater_constraints.begin() +
                                             i);
                modify = true;
                ++constraint_nb;
                i = 0;
                e = baryonyx::numeric_cast<int>(pb.greater_constraints.size());
            } else {
                ++i;
            }
        }

        if (modify)
            remove_empty_constraints(pb, constraint_nb);

    } while (modify);
}

bool
is_use_variable(baryonyx::problem& pb, int i)
{
    for (auto& cst : pb.equal_constraints)
        for (auto& elem : cst.elements)
            if (elem.variable_index == i)
                return true;

    for (auto& cst : pb.less_constraints)
        for (auto& elem : cst.elements)
            if (elem.variable_index == i)
                return true;

    for (auto& cst : pb.greater_constraints)
        for (auto& elem : cst.elements)
            if (elem.variable_index == i)
                return true;

    return false;
}

void
remove_unused_variables(const std::shared_ptr<baryonyx::context>& ctx,
                        baryonyx::problem& pb,
                        int& variable)
{
    std::vector<int> unused_variable;

    {
        std::unordered_set<int> vars;
        vars.reserve(pb.vars.values.size());

        for (int i = 0, e = (int)pb.vars.values.size(); i != e; ++i)
            vars.insert(i);

        for (auto& cst : pb.equal_constraints) {
            for (auto& elem : cst.elements) {
                auto it = vars.find(elem.variable_index);
                if (it != vars.end())
                    vars.erase(it);
            }
        }

        for (auto& cst : pb.less_constraints) {
            for (auto& elem : cst.elements) {
                auto it = vars.find(elem.variable_index);
                if (it != vars.end())
                    vars.erase(it);
            }
        }

        for (auto& cst : pb.greater_constraints) {
            for (auto& elem : cst.elements) {
                auto it = vars.find(elem.variable_index);
                if (it != vars.end())
                    vars.erase(it);
            }
        }

        unused_variable.assign(vars.begin(), vars.end());
    }

    std::sort(
      unused_variable.begin(), unused_variable.end(), std::greater<int>());

    const bool max = pb.type == baryonyx::objective_function_type::maximize;

    for (auto var : unused_variable) {
        ctx->debug(
          "      unused variable %s %d\n", pb.vars.names[var].c_str(), var);

        //
        // Here, we are sur that variables are not in constraint vectors so we
        // only remove the variables from the objective function but, we need
        // to know the factor and compute the correct value of this variable
        // according to the cost.
        //

        auto it = std::find_if(
          pb.objective.elements.begin(),
          pb.objective.elements.end(),
          [var](const auto& elem) { return elem.variable_index == var; });

        int value = 0;

        if (it != pb.objective.elements.end()) {
            if (it->factor > 0)
                value =
                  max ? pb.vars.values[var].max : pb.vars.values[var].min;
            else if (it->factor < 0)
                value =
                  max ? pb.vars.values[var].min : pb.vars.values[var].max;
        }

        remove_variable(pb, var, value);
    }

    variable += (int)unused_variable.size();
}

namespace baryonyx_private {

void
preprocess(const std::shared_ptr<baryonyx::context>& ctx,
           baryonyx::problem& pb)
{
    ctx->info("preprocessing:\n");

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

        pb.objective.elements =
          cleanup_function_element(pb.objective.elements, clean);

        ctx->info("    `-> %d function element(s) merged.\n", clean);
    }

    {
        ctx->info("  - cleaning already assigned variables:\n");
        int variable{ 0 }, constraint{ 0 };

        remove_assigned_variables(ctx, pb, variable, constraint);
        remove_duplicated_constraints(ctx, pb, constraint);
        remove_unused_variables(ctx, pb, variable);

        ctx->info("    `-> %d variable(s) - %d constraint(s) removed.\n",
                  variable,
                  constraint);
    }

#ifndef BARYONYX_FULL_OPTIMIZATION
    {
        ctx->info("  - write preprocessed.lp: ");
        std::ofstream ofs("preprocessed.lp");
        if (ofs.is_open()) {
            if (baryonyx_private::write_problem(ofs, pb))
                ctx->info("writing done\n");
            else
                ctx->info("writing fail\n");
        } else {
            ctx->info("opening failed\n");
        }
    }
#endif
}
}
