/* Copyright (C) 2016-2021 INRAE
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
#include <baryonyx/core-out>

#include "debug.hpp"
#include "itm-common.hpp"
#include "memory.hpp"
#include "utils.hpp"

#include <set>
#include <unordered_map>

namespace baryonyx {
namespace itm {

struct merged_constraint_hash
{
    inline size_t operator()(
      const std::vector<baryonyx::function_element>& fct) const noexcept
    {
        std::size_t seed{ fct.size() };

        for (auto& elem : fct)
            seed ^=
              elem.variable_index + 0x9e3779b9 + (seed << 6) + (seed >> 2);

        return seed;
    }
};

//
// Fills the cache and uses it to search same function and merge equal, less
// and greater equal constraints. To be use in any order, the fill operation
// depends of the operator type passed in parameter.
//
template<typename cacheT, typename retT>
static void
fill_merged_constraints(const baryonyx::context& ctx,
                        cacheT& cache,
                        baryonyx::operator_type op,
                        const baryonyx::problem& pb,
                        retT& ret)
{
    switch (op) {
    case baryonyx::operator_type::equal:
        for (const auto& elem : pb.equal_constraints) {
            auto it = cache.find(elem.elements);
            if (it == cache.end()) {
                cache.emplace(elem.elements, ret.size());
                ret.emplace_back(elem.elements,
                                 baryonyx::numeric_cast<int>(elem.value),
                                 baryonyx::numeric_cast<int>(elem.value),
                                 elem.id);
            } else {
                if (ret[it->second].min <= elem.value &&
                    elem.value <= ret[it->second].max) {
                    ret[it->second].min = elem.value;
                    ret[it->second].max = elem.value;
                } else {
                    error(ctx,
                          "Constraint {}: is inconsistent with previous {}.\n",
                          elem.id,
                          ret[it->second].id);

                    throw problem_definition_failure(
                      "equal",
                      problem_definition_error_tag::multiple_constraint);
                }
            }
        }
        break;
    case baryonyx::operator_type::less:
        for (const auto& elem : pb.less_constraints) {
            auto it = cache.find(elem.elements);
            if (it == cache.end()) {
                cache.emplace(elem.elements, ret.size());
                ret.emplace_back(elem.elements,
                                 std::numeric_limits<int>::min(),
                                 baryonyx::numeric_cast<int>(elem.value),
                                 elem.id);
            } else {
                ret[it->second].max =
                  std::min(ret[it->second].max,
                           baryonyx::numeric_cast<int>(elem.value));
            }
        }
        break;
    case baryonyx::operator_type::greater:
        for (const auto& elem : pb.greater_constraints) {
            auto it = cache.find(elem.elements);
            if (it == cache.end()) {
                cache.emplace(elem.elements, ret.size());
                ret.emplace_back(elem.elements,
                                 elem.value,
                                 std::numeric_limits<int>::max(),
                                 elem.id);
            } else {
                ret[it->second].min =
                  std::max(ret[it->second].min,
                           baryonyx::numeric_cast<int>(elem.value));
            }
        }
        break;
    default:
        break;
    }
}

/** @brief Build a merged constaints vector without any sort i.e. it restores
 *    the initial order of the raw problem.
 */
static std::vector<merged_constraint>
make_unsorted_merged_constraints(const baryonyx::context& ctx,
                                 const baryonyx::problem& pb)
{
    std::unordered_map<std::vector<baryonyx::function_element>,
                       std::size_t,
                       merged_constraint_hash>
      cache;

    std::vector<merged_constraint> ret;
    ret.reserve(pb.equal_constraints.size() + pb.less_constraints.size() +
                pb.greater_constraints.size());

    fill_merged_constraints(
      ctx, cache, baryonyx::operator_type::equal, pb, ret);
    fill_merged_constraints(
      ctx, cache, baryonyx::operator_type::less, pb, ret);
    fill_merged_constraints(
      ctx, cache, baryonyx::operator_type::greater, pb, ret);

    std::sort(ret.begin(), ret.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.id < rhs.id;
    });

    return ret;
}

/** @brief Build a merged constaints vector and sort constraints
 *    using the type of constraints (equal, less, greater or any order).
 */
static std::vector<merged_constraint>
make_ordered_merged_constraints(const baryonyx::context& ctx,
                                const baryonyx::problem& pb)
{
    std::unordered_map<std::vector<baryonyx::function_element>,
                       std::size_t,
                       merged_constraint_hash>
      cache;

    std::vector<merged_constraint> ret;
    ret.reserve(pb.equal_constraints.size() + pb.less_constraints.size() +
                pb.greater_constraints.size());

    baryonyx::operator_type op[3] = { baryonyx::operator_type::less,
                                      baryonyx::operator_type::greater,
                                      baryonyx::operator_type::equal };

    switch (ctx.parameters.pre_order) {
    case baryonyx::solver_parameters::pre_constraint_order::less_greater_equal:
        fill_merged_constraints(ctx, cache, op[0], pb, ret);
        fill_merged_constraints(ctx, cache, op[1], pb, ret);
        fill_merged_constraints(ctx, cache, op[2], pb, ret);
        break;
    case baryonyx::solver_parameters::pre_constraint_order::less_equal_greater:
        fill_merged_constraints(ctx, cache, op[0], pb, ret);
        fill_merged_constraints(ctx, cache, op[2], pb, ret);
        fill_merged_constraints(ctx, cache, op[1], pb, ret);
        break;
    case baryonyx::solver_parameters::pre_constraint_order::greater_less_equal:
        fill_merged_constraints(ctx, cache, op[1], pb, ret);
        fill_merged_constraints(ctx, cache, op[0], pb, ret);
        fill_merged_constraints(ctx, cache, op[2], pb, ret);
        break;
    case baryonyx::solver_parameters::pre_constraint_order::greater_equal_less:
        fill_merged_constraints(ctx, cache, op[1], pb, ret);
        fill_merged_constraints(ctx, cache, op[2], pb, ret);
        fill_merged_constraints(ctx, cache, op[0], pb, ret);
        break;
    case baryonyx::solver_parameters::pre_constraint_order::equal_less_greater:
        fill_merged_constraints(ctx, cache, op[2], pb, ret);
        fill_merged_constraints(ctx, cache, op[0], pb, ret);
        fill_merged_constraints(ctx, cache, op[1], pb, ret);
        break;
    case baryonyx::solver_parameters::pre_constraint_order::equal_greater_less:
        fill_merged_constraints(ctx, cache, op[2], pb, ret);
        fill_merged_constraints(ctx, cache, op[1], pb, ret);
        fill_merged_constraints(ctx, cache, op[0], pb, ret);
        break;
    default:
        break;
    }

    return ret;
}

//
// Reorder constraints according to the special operator proposed into the
// preprocessing parameter. If parameter is badly defined, we fall back to the
// unsorted_merged_constraints.
//
static std::vector<merged_constraint>
make_special_merged_constraints(const baryonyx::context& ctx,
                                const baryonyx::problem& pb)
{
    auto csts = make_unsorted_merged_constraints(ctx, pb);

    std::vector<int> variable_constraint_degree(pb.vars.values.size(), 0);
    for (const auto& cst : csts)
        for (const auto& elem : cst.elements)
            ++variable_constraint_degree[elem.variable_index];

    std::vector<double> variable_cost(pb.vars.values.size(), 0.0);
    for (const auto& elem : pb.objective.elements)
        variable_cost[elem.variable_index] =
          elem.factor /
          static_cast<double>(variable_constraint_degree[elem.variable_index]);

    std::vector<std::pair<int, double>> constraints_ratio_min(csts.size());
    for (int i = 0, ie = length(csts); i != ie; ++i) {
        constraints_ratio_min[i] = {
            i, variable_cost[csts[i].elements[0].variable_index]
        };

        for (int j = 1, ej = length(csts[i].elements); j != ej; ++j) {
            constraints_ratio_min[i].second =
              std::min(constraints_ratio_min[i].second,
                       variable_cost[csts[i].elements[j].variable_index]);
        }
    }

    if (ctx.parameters.pre_order ==
        solver_parameters::pre_constraint_order::p1) {
        std::sort(constraints_ratio_min.begin(),
                  constraints_ratio_min.end(),
                  [](const auto& lhs, const auto& rhs) {
                      return lhs.second > rhs.second;
                  });
    } else {
        std::sort(constraints_ratio_min.begin(),
                  constraints_ratio_min.end(),
                  [](const auto& lhs, const auto& rhs) {
                      return lhs.second < rhs.second;
                  });
    }

    for (int i = 0, e = length(csts); i != e; ++i)
        csts[i].id = constraints_ratio_min[i].first;

    std::sort(csts.begin(), csts.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.id < rhs.id;
    });

    for (int i = 0, e = length(csts); i != e; ++i)
        std::sort(csts[i].elements.begin(),
                  csts[i].elements.end(),
                  [&variable_cost](const auto& lhs, const auto& rhs) {
                      return variable_cost[lhs.variable_index] >
                             variable_cost[rhs.variable_index];
                  });

    return csts;
}

static void
improve_memory_usage(std::vector<merged_constraint>& csts)
{
    for (auto& cst : csts)
        std::sort(cst.elements.begin(),
                  cst.elements.end(),
                  [](const auto& lhs, const auto& rhs) {
                      return lhs.variable_index < rhs.variable_index;
                  });
}

std::vector<merged_constraint>
make_merged_constraints(const context& ctx, const problem& pb)
{
    info(ctx,
         "  - merge constraint according to the `{}` algorithm\n",
         ctx.parameters.pre_order);

    auto original_nb = static_cast<int>(pb.equal_constraints.size() +
                                        pb.less_constraints.size() +
                                        pb.greater_constraints.size());

    std::vector<merged_constraint> ret;
    ret.reserve(original_nb);

    switch (ctx.parameters.pre_order) {
    case solver_parameters::pre_constraint_order::none:
        ret = make_unsorted_merged_constraints(ctx, pb);
        break;
    case solver_parameters::pre_constraint_order::memory:
        ret = make_unsorted_merged_constraints(ctx, pb);
        improve_memory_usage(ret);
        break;
    case solver_parameters::pre_constraint_order::less_greater_equal:
    case solver_parameters::pre_constraint_order::less_equal_greater:
    case solver_parameters::pre_constraint_order::greater_less_equal:
    case solver_parameters::pre_constraint_order::greater_equal_less:
    case solver_parameters::pre_constraint_order::equal_less_greater:
    case solver_parameters::pre_constraint_order::equal_greater_less:
        ret = make_ordered_merged_constraints(ctx, pb);
        improve_memory_usage(ret);
        break;
    case solver_parameters::pre_constraint_order::p1:
    case solver_parameters::pre_constraint_order::p2:
    case solver_parameters::pre_constraint_order::p3:
    case solver_parameters::pre_constraint_order::p4:
        ret = make_special_merged_constraints(ctx, pb);
        break;
    }

    info(ctx,
         "  - merged constraints removed: {}\n"
         "  - problem memory used: {}\n",
         original_nb - static_cast<int>(ret.size()),
         to_string(memory_consumed_size(memory_consumed(pb))));

    return ret;
}

} // namespace itm
} // namespace baryonyx
