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

#include <optional>
#include <set>
#include <unordered_set>

namespace baryonyx {
namespace itm {

/**
 * Positions of the rows already in the merged constraints, hashed and
 * compared through the elements of these rows: a key copying the elements
 * would allocate it again for each row of each solver run.
 */
class row_index
{
public:
    explicit row_index(const std::vector<merged_constraint>& rows,
                       const std::size_t capacity)
      : m_rows(0, hash{ &rows }, equal{ &rows })
    {
        m_rows.reserve(capacity);
    }

    /** @return the position of the row with these elements, if any. */
    std::optional<std::size_t> find(
      const std::vector<baryonyx::function_element>& elements) const
    {
        const auto it = m_rows.find(elements);

        if (it == m_rows.end())
            return std::nullopt;

        return *it;
    }

    void insert(const std::size_t position)
    {
        m_rows.insert(position);
    }

private:
    using elements_type = std::vector<baryonyx::function_element>;

    struct hash
    {
        using is_transparent = void;

        const std::vector<merged_constraint>* rows;

        std::size_t operator()(const elements_type& fct) const noexcept
        {
            std::size_t seed{ fct.size() };

            for (auto& elem : fct)
                seed ^= static_cast<std::size_t>(elem.variable_index) +
                        0x9e3779b9 + (seed << 6) + (seed >> 2);

            return seed;
        }

        std::size_t operator()(const std::size_t position) const noexcept
        {
            return (*this)((*rows)[position].elements);
        }
    };

    struct equal
    {
        using is_transparent = void;

        const std::vector<merged_constraint>* rows;

        const elements_type& get(const std::size_t position) const noexcept
        {
            return (*rows)[position].elements;
        }

        static const elements_type& get(const elements_type& elements) noexcept
        {
            return elements;
        }

        template<typename Lhs, typename Rhs>
        bool operator()(const Lhs& lhs, const Rhs& rhs) const noexcept
        {
            return get(lhs) == get(rhs);
        }
    };

    std::unordered_set<std::size_t, hash, equal> m_rows;
};

//
// Uses the index to search same function and merge equal, less and greater
// equal constraints. To be use in any order, the fill operation depends of
// the operator type passed in parameter.
//
template<typename retT>
static void
fill_merged_constraints(const baryonyx::context& ctx,
                        row_index& cache,
                        baryonyx::operator_type op,
                        const baryonyx::problem& pb,
                        retT& ret)
{
    const auto add = [&](const baryonyx::constraint& elem, int min, int max) {
        ret.emplace_back(elem.elements, min, max, elem.id);
        cache.insert(ret.size() - 1u);
    };

    switch (op) {
    case baryonyx::operator_type::equal:
        for (const auto& elem : pb.equal_constraints) {
            const auto found = cache.find(elem.elements);
            if (!found) {
                add(elem,
                    baryonyx::numeric_cast<int>(elem.value),
                    baryonyx::numeric_cast<int>(elem.value));
            } else {
                auto& row = ret[*found];
                if (row.min <= elem.value && elem.value <= row.max) {
                    row.min = elem.value;
                    row.max = elem.value;
                } else {
                    error(ctx,
                          "Constraint {}: is inconsistent with previous {}.\n",
                          elem.id,
                          row.id);

                    throw problem_definition_failure(
                      "equal",
                      problem_definition_error_tag::multiple_constraint);
                }
            }
        }
        break;
    case baryonyx::operator_type::less:
        for (const auto& elem : pb.less_constraints) {
            const auto found = cache.find(elem.elements);
            if (!found) {
                add(elem,
                    std::numeric_limits<int>::min(),
                    baryonyx::numeric_cast<int>(elem.value));
            } else {
                ret[*found].max = std::min(
                  ret[*found].max, baryonyx::numeric_cast<int>(elem.value));
            }
        }
        break;
    case baryonyx::operator_type::greater:
        for (const auto& elem : pb.greater_constraints) {
            const auto found = cache.find(elem.elements);
            if (!found) {
                add(elem, elem.value, std::numeric_limits<int>::max());
            } else {
                ret[*found].min = std::max(
                  ret[*found].min, baryonyx::numeric_cast<int>(elem.value));
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
    const auto rows = pb.equal_constraints.size() +
                      pb.less_constraints.size() +
                      pb.greater_constraints.size();

    std::vector<merged_constraint> ret;
    ret.reserve(rows);
    row_index cache(ret, rows);

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
