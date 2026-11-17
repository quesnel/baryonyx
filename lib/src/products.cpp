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

#include "products.hpp"
#include "private.hpp"
#include "utils.hpp"

#include <cstdint>
#include <unordered_map>
#include <utility>

namespace baryonyx {

namespace {

enum product_bound : std::uint8_t
{
    upper_bound = 1, // z <= x
    lower_bound = 2  // z >= x + y - 1
};

/** How a variable is linked to the rows that define it as a product. */
struct product_candidate
{
    int rows = 0;
    int product_rows = 0;
    int a = -1;
    int b = -1;
    std::uint8_t bounds = 0;
    bool invalid = false;

    void add_partner(const int partner) noexcept
    {
        if (a == partner || b == partner)
            return;

        if (a < 0)
            a = partner;
        else if (b < 0)
            b = partner;
        else
            invalid = true;
    }
};

// Reads a row in its "<=" form: z - x <= 0 bounds z by x, and
// x + y - z <= 1 bounds z from below by x + y - 1.
void
classify_row(const constraint& cst,
             const bool greater,
             std::vector<product_candidate>& candidates)
{
    const auto size = cst.elements.size();
    if (size != 2 && size != 3)
        return;

    const int sign = greater ? -1 : 1;
    const int rhs = sign * cst.value;
    int positives = 0;
    int negatives = 0;

    for (const auto& elem : cst.elements) {
        const int factor = sign * elem.factor;
        if (factor == 1)
            ++positives;
        else if (factor == -1)
            ++negatives;
        else
            return;
    }

    const auto find = [&](const int factor) {
        for (const auto& elem : cst.elements)
            if (sign * elem.factor == factor)
                return elem.variable_index;
        return -1;
    };

    if (size == 2 && positives == 1 && rhs == 0) {
        auto& z = candidates[static_cast<std::size_t>(find(1))];
        ++z.product_rows;
        z.bounds |= upper_bound;
        z.add_partner(find(-1));
    } else if (size == 3 && negatives == 1 && rhs == 1) {
        const int index = find(-1);
        auto& z = candidates[static_cast<std::size_t>(index)];
        ++z.product_rows;
        z.bounds |= lower_bound;
        for (const auto& elem : cst.elements)
            if (elem.variable_index != index)
                z.add_partner(elem.variable_index);
    }
}

} // anonymous namespace

std::optional<raw_problem>
linearized_products::reduce(const context& ctx, const raw_problem& pb)
{
    m_products.clear();

    const auto n = pb.vars.values.size();
    std::vector<product_candidate> candidates(n);
    std::vector<double> linear(n, 0.0);
    std::vector<char> quadratic(n, false);

    for (const auto* csts : { &pb.equal_constraints,
                              &pb.greater_constraints,
                              &pb.less_constraints })
        for (const auto& cst : *csts)
            for (const auto& elem : cst.elements)
                ++candidates[static_cast<std::size_t>(elem.variable_index)]
                    .rows;

    for (const auto& cst : pb.greater_constraints)
        classify_row(cst, true, candidates);
    for (const auto& cst : pb.less_constraints)
        classify_row(cst, false, candidates);

    const double sense =
      pb.type == objective_function_type::maximize ? -1.0 : 1.0;
    for (const auto& elem : pb.objective.elements)
        linear[static_cast<std::size_t>(elem.variable_index)] +=
          sense * elem.factor;
    for (const auto& elem : pb.objective.qelements) {
        quadratic[static_cast<std::size_t>(elem.variable_index_a)] = true;
        quadratic[static_cast<std::size_t>(elem.variable_index_b)] = true;
    }

    // With only z <= x and z <= y, the optimum takes z = x.y only if z does
    // not cost; with only z >= x + y - 1, only if z does not pay. Otherwise
    // the rows do not define a product and z is kept.
    std::vector<char> removed(n, false);
    std::vector<char> partner(n, false);
    for (std::size_t z = 0; z != n; ++z) {
        const auto& c = candidates[z];
        const auto& value = pb.vars.values[z];

        if (c.invalid || c.rows == 0 || c.product_rows != c.rows ||
            c.a < 0 || c.b < 0 || quadratic[z] || partner[z] ||
            value.type == variable_type::real || value.min != 0 ||
            value.max != 1)
            continue;

        if ((c.bounds == upper_bound && linear[z] > 0.0) ||
            (c.bounds == lower_bound && linear[z] < 0.0))
            continue;

        const auto a = static_cast<std::size_t>(c.a);
        const auto b = static_cast<std::size_t>(c.b);
        if (removed[a] || removed[b])
            continue;

        removed[z] = true;
        partner[a] = true;
        partner[b] = true;
        m_products.push_back(
          { pb.vars.names[z], pb.vars.names[a], pb.vars.names[b] });
    }

    if (m_products.empty())
        return std::nullopt;

    raw_problem ret;
    ret.strings = pb.strings;
    ret.type = pb.type;
    ret.status = pb.status;

    std::vector<int> index(n, -1);
    for (std::size_t i = 0; i != n; ++i) {
        if (removed[i])
            continue;

        index[i] = length(ret.vars.names);
        ret.vars.names.emplace_back(pb.vars.names[i]);
        ret.vars.values.emplace_back(pb.vars.values[i]);
    }

    const auto map = [&index](const int variable) {
        return index[static_cast<std::size_t>(variable)];
    };

    ret.objective.value = pb.objective.value;
    for (const auto& elem : pb.objective.elements)
        if (const int i = map(elem.variable_index); i >= 0)
            ret.objective.elements.emplace_back(elem.factor, i);

    // One term per pair of variables, as the parser builds them.
    std::unordered_map<std::uint64_t, std::size_t> terms;
    const auto add_term = [&](const double factor, int a, int b) {
        if (a > b)
            std::swap(a, b);
        const auto key = (static_cast<std::uint64_t>(a) << 32) |
                         static_cast<std::uint32_t>(b);
        const auto [it, inserted] =
          terms.try_emplace(key, ret.objective.qelements.size());
        if (inserted)
            ret.objective.qelements.emplace_back(factor, a, b);
        else
            ret.objective.qelements[it->second].factor += factor;
    };

    for (const auto& elem : pb.objective.qelements)
        add_term(
          elem.factor, map(elem.variable_index_a), map(elem.variable_index_b));

    for (const auto& elem : pb.objective.elements) {
        const auto z = static_cast<std::size_t>(elem.variable_index);
        if (removed[z] && elem.factor != 0.0)
            add_term(elem.factor, map(candidates[z].a), map(candidates[z].b));
    }

    int removed_rows = 0;
    const auto copy = [&](const std::vector<constraint>& from,
                          std::vector<constraint>& to) {
        for (const auto& cst : from) {
            bool keep = true;
            for (const auto& elem : cst.elements)
                keep = keep && map(elem.variable_index) >= 0;

            if (!keep) {
                ++removed_rows;
                continue;
            }

            auto& out = to.emplace_back();
            out.label = cst.label;
            out.value = cst.value;
            out.id = cst.id;
            out.elements.reserve(cst.elements.size());
            for (const auto& elem : cst.elements)
                out.elements.emplace_back(elem.factor,
                                          map(elem.variable_index));
        }
    };

    copy(pb.equal_constraints, ret.equal_constraints);
    copy(pb.greater_constraints, ret.greater_constraints);
    copy(pb.less_constraints, ret.less_constraints);

    info(ctx,
         "  - linearized products: {} variables and {} constraints replaced "
         "by quadratic terms\n",
         m_products.size(),
         removed_rows);

    return ret;
}

void
linearized_products::restore(result& r) const
{
    if (m_products.empty())
        return;

    std::unordered_map<std::string_view, std::size_t> position;
    for (std::size_t i = 0, e = r.variable_name.size(); i != e; ++i)
        position.emplace(r.variable_name[i], i);

    std::unordered_map<std::string_view, bool> fixed;
    for (std::size_t i = 0, e = r.affected_vars.names.size(); i != e; ++i)
        fixed.emplace(r.affected_vars.names[i], r.affected_vars.values[i]);

    bool complete = !r.solutions.empty();
    for (const auto& s : r.solutions)
        complete = complete && s.variables.size() == r.variable_name.size();

    for (const auto& p : m_products) {
        const auto pa = position.find(p.a);
        const auto pb = position.find(p.b);

        if (pa == position.end() && pb == position.end()) {
            const auto fa = fixed.find(p.a);
            const auto fb = fixed.find(p.b);
            if (fa != fixed.end() && fb != fixed.end())
                r.affected_vars.push_back(p.name, fa->second && fb->second);
            continue;
        }

        if (!complete)
            continue;

        const auto value = [&](const solution& s,
                               const auto& it,
                               std::string_view name) {
            if (it != position.end())
                return s.variables[it->second] != 0;

            const auto f = fixed.find(name);
            return f != fixed.end() && f->second;
        };

        r.variable_name.emplace_back(p.name);
        for (auto& s : r.solutions)
            s.variables.emplace_back(
              static_cast<std::int8_t>(value(s, pa, p.a) && value(s, pb, p.b)));
    }
}

} // namespace baryonyx
