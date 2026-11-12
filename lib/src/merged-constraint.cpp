/* Copyright (C) 2016-2026 INRAE
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

// The solvers assume each row can reach its bounds with 0-1 variables; a row
// that cannot would break that invariant at the first update instead of
// being reported.
static void
check_satisfiable(const context& ctx,
                  const std::vector<merged_constraint>& csts)
{
    for (const auto& cst : csts) {
        long long lower = 0;
        long long upper = 0;

        for (const auto& elem : cst.elements)
            (elem.factor < 0 ? lower : upper) += elem.factor;

        if (upper < cst.min || cst.max < lower) {
            error(ctx,
                  "Constraint {}: bounds [{}, {}] cannot be reached with 0-1 "
                  "variables (reachable range [{}, {}]).\n",
                  cst.id,
                  cst.min,
                  cst.max,
                  lower,
                  upper);

            throw solver_failure(solver_error_tag::unrealisable_constraint);
        }
    }
}

std::vector<merged_constraint>
make_merged_constraints(const context& ctx, const problem& pb)
{
    auto original_nb = static_cast<int>(pb.equal_constraints.size() +
                                        pb.less_constraints.size() +
                                        pb.greater_constraints.size());

    auto ret = make_unsorted_merged_constraints(ctx, pb);

    info(ctx,
         "  - merged constraints removed: {}\n"
         "  - problem memory used: {}\n",
         original_nb - static_cast<int>(ret.size()),
         to_string(memory_consumed_size(memory_consumed(pb))));

    check_satisfiable(ctx, ret);

    return ret;
}

} // namespace itm
} // namespace baryonyx
