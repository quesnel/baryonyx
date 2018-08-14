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

#include <baryonyx/core-compare>

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
fill_merged_constraints(const baryonyx::context_ptr& ctx,
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

//
// Restore original order: reorder merged constraints according to the
// merged_constraint::id (initially initialized in the lp format reading
// process).
//
static std::vector<merged_constraint>
make_unsorted_merged_constraints(const baryonyx::context_ptr& ctx,
                                 const baryonyx::problem& pb)
{
    info(ctx, "  - merge constraints without any sort\n");

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

//
// Reorder constraints according to the type of operator proposed into the
// preprocessing parameter. If parameter is badly defined, we fall back to the
// unsorted_merged_constraints.
//
static std::vector<merged_constraint>
make_ordered_merged_constraints(const baryonyx::context_ptr& ctx,
                                const baryonyx::problem& pb)
{
    bx_expects(static_cast<int>(ctx->parameters.pre_order) >= 5);

    info(ctx, "  - merge constraints with ordered sort\n");

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

    switch (ctx->parameters.pre_order) {
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
make_special_merged_constraints(const baryonyx::context_ptr& ctx,
                                const baryonyx::problem& pb)
{
    bx_expects(static_cast<int>(ctx->parameters.pre_order) >= 1);
    bx_expects(static_cast<int>(ctx->parameters.pre_order) <= 4);

    info(ctx,
         "  - merge constraint according to the `{}` algorithm\n",
         baryonyx::to_string(ctx->parameters.pre_order));

    auto ret = make_unsorted_merged_constraints(ctx, pb);

    std::vector<int> vars(pb.vars.values.size(), 0);
    std::vector<std::set<int>> linkvars(pb.vars.values.size());
    std::vector<std::set<int>> linkcst(pb.vars.values.size());

    for (auto it{ ret.begin() }, et{ ret.end() }; it != et; ++it) {
        for (std::size_t i{ 0 }, e{ it->elements.size() }; i != e; ++i) {

            linkcst[it->elements[i].variable_index].emplace(
              numeric_cast<int>(std::distance(it, ret.end())));

            for (std::size_t j{ 0 }; j != e; ++j) {
                if (i != j)
                    linkvars[it->elements[i].variable_index].emplace(
                      it->elements[j].variable_index);
            }
        }

        for (auto& f : it->elements) {
            vars[f.variable_index]++;
        }
    }

    std::vector<std::pair<merged_constraint, int>> tosort;

    if (ctx->parameters.pre_order ==
        solver_parameters::pre_constraint_order::variables_number) {
        for (auto& cst : ret) {
            tosort.emplace_back(cst, 0);
            for (auto& elem : cst.elements) {
                for (auto& s : linkcst[elem.variable_index])
                    tosort.back().second +=
                      static_cast<int>(linkvars[s].size());
            }
        }

        std::sort(
          tosort.begin(), tosort.end(), [](const auto& lhs, const auto& rhs) {
              return rhs.second < lhs.second;
          });

        ret.clear();
        std::transform(tosort.begin(),
                       tosort.end(),
                       std::back_inserter(ret),
                       [](const auto& elem) { return elem.first; });

        return ret;
    }

    if (ctx->parameters.pre_order ==
        solver_parameters::pre_constraint_order::variables_weight) {
        std::sort(
          ret.begin(), ret.end(), [vars](const auto& lhs, const auto& rhs) {
              int sumlhs{ 0 };
              int sumrhs{ 0 };

              for (auto& f : lhs.elements)
                  sumlhs += vars[f.variable_index];

              for (auto& f : rhs.elements)
                  sumrhs += vars[f.variable_index];

              return sumlhs < sumrhs;
          });

        return ret;
    }

    if (ctx->parameters.pre_order ==
        solver_parameters::pre_constraint_order::constraints_weight) {
        std::sort(ret.begin(),
                  ret.end(),
                  [linkvars, linkcst](const auto& lhs, const auto& rhs) {
                      std::size_t sumlhs{ 1 };
                      std::size_t sumrhs{ 1 };

                      for (auto& f : lhs.elements)
                          sumlhs *= linkcst[f.variable_index].size();

                      for (auto& f : rhs.elements)
                          sumrhs *= linkcst[f.variable_index].size();

                      return sumrhs > sumlhs;
                  });

        return ret;
    }

    if (ctx->parameters.pre_order ==
        solver_parameters::pre_constraint_order::implied) {
        // Algorithm to build a tosort vector according to constraints of type:
        // -x1 -x2 -x3 +x4 <= 0

        for (auto& cst : ret) {
            tosort.emplace_back(cst, 0);

            int nbneg{ 0 }, nbpos{ 0 };

            for (auto& elem : cst.elements)
                if (elem.factor < 0)
                    nbneg++;
                else
                    nbpos++;

            if (((nbneg > 1 && nbpos == 1) || (nbpos > 1 && nbneg == 1)) &&
                cst.min == cst.max) {
                tosort.back().second = 1;
            }
        }

        std::sort(
          tosort.begin(), tosort.end(), [](const auto& lhs, const auto& rhs) {
              return rhs.second < lhs.second;
          });

        ret.clear();
        std::transform(tosort.begin(),
                       tosort.end(),
                       std::back_inserter(ret),
                       [](const auto& elem) { return elem.first; });

        return ret;
    }

    return ret;
}

std::vector<merged_constraint>
make_merged_constraints(const context_ptr& ctx, const problem& pb)
{
    auto original_nb = static_cast<int>(pb.equal_constraints.size() +
                                        pb.less_constraints.size() +
                                        pb.greater_constraints.size());

    std::vector<merged_constraint> ret;
    ret.reserve(original_nb);

    if (ctx->parameters.pre_order ==
        solver_parameters::pre_constraint_order::none) {
        ret = make_unsorted_merged_constraints(ctx, pb);
    } else if (static_cast<int>(ctx->parameters.pre_order) >= 5) {
        ret = make_ordered_merged_constraints(ctx, pb);
    } else {
        ret = make_special_merged_constraints(ctx, pb);
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
