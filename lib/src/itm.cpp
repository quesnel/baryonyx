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

#include "itm.hpp"

#include <set>

#include <cassert>

namespace bx = baryonyx;

struct merged_constraint_hash
{
    inline size_t operator()(
      const std::vector<bx::function_element>& fct) const noexcept
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
fill_merged_constraints(const std::shared_ptr<bx::context>& ctx,
                        cacheT& cache,
                        bx::operator_type op,
                        const bx::problem& pb,
                        retT& ret)
{
    switch (op) {
    case bx::operator_type::equal:
        for (const auto& elem : pb.equal_constraints) {
            auto it = cache.find(elem.elements);
            if (it == cache.end()) {
                cache.emplace(elem.elements, ret.size());
                ret.emplace_back(elem.elements,
                                 bx::numeric_cast<int>(elem.value),
                                 bx::numeric_cast<int>(elem.value),
                                 elem.id);
            } else {
                if (ret[it->second].min <= elem.value and
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
    case bx::operator_type::less:
        for (const auto& elem : pb.less_constraints) {
            auto it = cache.find(elem.elements);
            if (it == cache.end()) {
                cache.emplace(elem.elements, ret.size());
                ret.emplace_back(elem.elements,
                                 std::numeric_limits<int>::min(),
                                 bx::numeric_cast<int>(elem.value),
                                 elem.id);
            } else {
                ret[it->second].max = std::min(
                  ret[it->second].max, bx::numeric_cast<int>(elem.value));
            }
        }
        break;
    case bx::operator_type::greater:
        for (const auto& elem : pb.greater_constraints) {
            auto it = cache.find(elem.elements);
            if (it == cache.end()) {
                cache.emplace(elem.elements, ret.size());
                ret.emplace_back(elem.elements,
                                 elem.value,
                                 std::numeric_limits<int>::max(),
                                 elem.id);
            } else {
                ret[it->second].min = std::max(
                  ret[it->second].min, bx::numeric_cast<int>(elem.value));
            }
        }
        break;
    default:
        break;
    }
}

static inline bx::operator_type
operator_from_string(const std::string& op) noexcept
{
    if (op == "equal")
        return bx::operator_type::equal;
    if (op == "less")
        return bx::operator_type::less;
    if (op == "greater")
        return bx::operator_type::greater;

    return bx::operator_type::undefined;
}

//
// Restore original order: reorder merged constraints according to the
// merged_constraint::id (initially initialized in the lp format reading
// process).
//
static std::vector<bx::itm::merged_constraint>
make_unsorted_merged_constraints(const std::shared_ptr<bx::context>& ctx,
                                 const bx::problem& pb)
{
    info(ctx, "  - merge constraints without any sort\n");

    std::unordered_map<std::vector<bx::function_element>,
                       std::size_t,
                       merged_constraint_hash>
      cache;

    std::vector<bx::itm::merged_constraint> ret;
    ret.reserve(pb.equal_constraints.size() + pb.less_constraints.size() +
                pb.greater_constraints.size());

    fill_merged_constraints(ctx, cache, bx::operator_type::equal, pb, ret);
    fill_merged_constraints(ctx, cache, bx::operator_type::less, pb, ret);
    fill_merged_constraints(ctx, cache, bx::operator_type::greater, pb, ret);

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
static std::vector<bx::itm::merged_constraint>
make_ordered_merged_constraints(const std::shared_ptr<bx::context>& ctx,
                                const bx::problem& pb,
                                const bx::itm::parameters& p)
{
    info(ctx, "  - merge constraints with ordered sort\n");

    auto first = p.preprocessing.find(',', 0);
    if (first != std::string::npos) {
        auto second = p.preprocessing.find(',', first + 1);
        if (second != std::string::npos) {
            auto op0 = p.preprocessing.substr(0, first);
            auto op1 = p.preprocessing.substr(first + 1, second - first - 1);
            auto op2 = p.preprocessing.substr(second + 1);

            bx::operator_type op[3] = { bx::operator_type::less,
                                        bx::operator_type::greater,
                                        bx::operator_type::equal };

            op[0] = operator_from_string(op0);
            op[1] = operator_from_string(op1);
            op[2] = operator_from_string(op2);

            if (op[0] != bx::operator_type::undefined and
                op[1] != bx::operator_type::undefined and
                op[2] != bx::operator_type::undefined and op[0] != op[1] and
                op[1] != op[2] and op[0] != op[2]) {

                std::unordered_map<std::vector<bx::function_element>,
                                   std::size_t,
                                   merged_constraint_hash>
                  cache;

                std::vector<bx::itm::merged_constraint> ret;
                ret.reserve(pb.equal_constraints.size() +
                            pb.less_constraints.size() +
                            pb.greater_constraints.size());

                fill_merged_constraints(ctx, cache, op[0], pb, ret);
                fill_merged_constraints(ctx, cache, op[1], pb, ret);
                fill_merged_constraints(ctx, cache, op[2], pb, ret);

                return ret;
            }
        }

        warning(ctx,
                "    Bad preprocessing order mode: `{}`. Falling back "
                "to the unsorted merged constraint\n",
                p.preprocessing);
    }

    return make_unsorted_merged_constraints(ctx, pb);
}

//
// Reorder constraints according to the special operator proposed into the
// preprocessing parameter. If parameter is badly defined, we fall back to the
// unsorted_merged_constraints.
//
static std::vector<bx::itm::merged_constraint>
make_special_merged_constraints(const std::shared_ptr<bx::context>& ctx,
                                const bx::problem& pb,
                                const bx::itm::parameters& p)
{
    info(ctx,
         "  - merge constraint according to the `{}` algorithm\n",
         p.preprocessing);

    auto ret = make_unsorted_merged_constraints(ctx, pb);

    std::vector<int> vars(pb.vars.values.size(), 0);
    std::vector<std::set<int>> linkvars(pb.vars.values.size());
    std::vector<std::set<int>> linkcst(pb.vars.values.size());

    for (auto it{ ret.begin() }, et{ ret.end() }; it != et; ++it) {
        for (std::size_t i{ 0 }, e{ it->elements.size() }; i != e; ++i) {

            linkcst[it->elements[i].variable_index].emplace(
              std::distance(it, ret.end()));

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

    std::vector<std::pair<bx::itm::merged_constraint, int>> tosort;

    if (p.preprocessing == "variables-number") {
        for (auto& cst : ret) {
            tosort.emplace_back(cst, 0);
            for (auto& elem : cst.elements) {
                for (auto& s : linkcst[elem.variable_index])
                    tosort.back().second += linkvars[s].size();
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

    if (p.preprocessing == "variables-weight") {
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

    if (p.preprocessing == "constraints-weight") {
        std::sort(ret.begin(),
                  ret.end(),
                  [linkvars, linkcst](const auto& lhs, const auto& rhs) {
                      int sumlhs{ 1 };
                      int sumrhs{ 1 };

                      for (auto& f : lhs.elements)
                          sumlhs *= linkcst[f.variable_index].size();

                      for (auto& f : rhs.elements)
                          sumrhs *= linkcst[f.variable_index].size();

                      return sumrhs > sumlhs;
                  });

        return ret;
    }

    if (p.preprocessing == "implied") {
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

            if (((nbneg > 1 and nbpos == 1) or (nbpos > 1 and nbneg == 1)) and
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

namespace baryonyx {
namespace itm {

std::vector<merged_constraint>
make_merged_constraints(const std::shared_ptr<context>& ctx,
                        const problem& pb,
                        const parameters& p)
{
    auto original_nb = static_cast<int>(pb.equal_constraints.size() +
                                        pb.less_constraints.size() +
                                        pb.greater_constraints.size());

    std::vector<merged_constraint> ret;
    ret.reserve(original_nb);

    if (p.preprocessing.empty() or p.preprocessing == "none") {
        ret = make_unsorted_merged_constraints(ctx, pb);
    } else if (p.preprocessing.find(',') != std::string::npos) {
        ret = make_ordered_merged_constraints(ctx, pb, p);
    } else {
        ret = make_special_merged_constraints(ctx, pb, p);
    }

    info(ctx,
         "  - merged constraints removed: {}\n",
         original_nb - static_cast<int>(ret.size()));

    return ret;
}
}
}
