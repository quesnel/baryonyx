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
#include <baryonyx/core>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <numeric>
#include <span>
#include <string_view>
#include <tuple>
#include <vector>

#include "debug.hpp"
#include "memory.hpp"
#include "private.hpp"
#include "problem.hpp"
#include "utils.hpp"

namespace {

namespace bx = baryonyx;

enum class constraint_kind
{
    equal,
    greater,
    less
};

constexpr std::array<constraint_kind, 3> constraint_kinds{
    constraint_kind::equal,
    constraint_kind::greater,
    constraint_kind::less
};

constexpr std::string_view constraint_kind_names[] = { "equal",
                                                       "greater",
                                                       "less" };

constexpr std::size_t
to_index(constraint_kind kind) noexcept
{
    return static_cast<std::size_t>(kind);
}

constexpr bool
is_satisfied(constraint_kind kind, long long lhs, long long rhs) noexcept
{
    switch (kind) {
    case constraint_kind::equal:
        return lhs == rhs;
    case constraint_kind::greater:
        return lhs >= rhs;
    case constraint_kind::less:
        return lhs <= rhs;
    }

    return false;
}

[[maybe_unused]] static bool
is_valid_variable_index(const std::vector<bx::constraint>& csts,
                        int size) noexcept
{
    for (const auto& cst : csts)
        for (auto elem : cst.elements)
            if (!(elem.variable_index >= 0 && elem.variable_index < size))
                return false;

    return true;
}

/// The constraints of one kind in which each variable appears, once per
/// element, stored in two flat arrays: a vector per variable would make one
/// allocation per variable for each preprocessing.
class variable_constraints
{
public:
    variable_constraints(const std::vector<bx::constraint>& csts,
                         const std::size_t variables)
      : m_start(variables + 1u, 0)
    {
        for (const auto& cst : csts)
            for (const auto& elem : cst.elements)
                ++m_start[static_cast<std::size_t>(elem.variable_index) + 1u];

        for (std::size_t j = 1; j < m_start.size(); ++j)
            m_start[j] += m_start[j - 1];

        m_constraints.resize(static_cast<std::size_t>(m_start.back()));
        auto next = m_start;
        for (int i = 0, e = bx::length(csts); i != e; ++i)
            for (const auto& elem : csts[static_cast<std::size_t>(i)].elements)
                m_constraints[static_cast<std::size_t>(
                  next[static_cast<std::size_t>(elem.variable_index)]++)] = i;
    }

    std::span<const int> operator()(const int variable) const noexcept
    {
        const auto first = m_start[static_cast<std::size_t>(variable)];
        const auto last = m_start[static_cast<std::size_t>(variable) + 1u];

        return { m_constraints.data() + first,
                 static_cast<std::size_t>(last - first) };
    }

private:
    std::vector<int> m_start;
    std::vector<int> m_constraints;
};

/// Propagates the 0-1 assignments implied by bounds and constraints. A
/// constraint is tracked by the number of its elements whose variable has not
/// been propagated yet; once at most one remains, the constraint is decided
/// exactly and removed. Any contradiction makes the problem infeasible and is
/// reported with a @c solver_failure.
template<typename Problem>
class preprocessor
{
private:
    static constexpr std::int8_t unassigned = -1;

    const bx::context& ctx;
    const Problem& pb;
    std::array<const std::vector<bx::constraint>*, 3> constraints;
    std::array<std::vector<int>, 3> remaining;
    std::array<variable_constraints, 3> cache;
    std::vector<std::int8_t> values;
    std::vector<int> to_propagate;
    std::array<std::vector<int>, 3> touched;
    std::array<std::vector<char>, 3> is_touched;
    int assigned = 0;

    const bx::constraint& get(constraint_kind kind, int index) const noexcept
    {
        return (*constraints[to_index(kind)])[index];
    }

    [[noreturn]] void unrealisable(constraint_kind kind, int index) const
    {
        const auto& cst = get(kind, index);

        error(ctx,
              "  - Preprocessor: {} constraint {} (id: {}) cannot be "
              "satisfied with 0-1 variables.\n",
              constraint_kind_names[to_index(kind)],
              cst.label,
              cst.id);

        throw bx::solver_failure(
          bx::solver_error_tag::unrealisable_constraint);
    }

    [[noreturn]] void unrealisable(int variable) const
    {
        error(ctx,
              "  - Preprocessor: variable {} cannot take a 0-1 value.\n",
              pb.vars.names[variable]);

        throw bx::solver_failure(
          bx::solver_error_tag::unrealisable_constraint);
    }

    [[nodiscard]] bool assign(int variable, bool value)
    {
        if (values[variable] == unassigned) {
            values[variable] = value;
            ++assigned;
            to_propagate.emplace_back(variable);

            debug(ctx,
                  "    - variable {} assigned to {}.\n",
                  pb.vars.names[variable],
                  value);

            return true;
        }

        return values[variable] == static_cast<std::int8_t>(value);
    }

    // At most one element remains unassigned: enumerate its two values.
    void decide(constraint_kind kind, int index)
    {
        const auto& cst = get(kind, index);
        remaining[to_index(kind)][index] = 0;

        long long rhs = cst.value;
        int free_variable = -1;
        int free_factor = 0;

        for (const auto& elem : cst.elements) {
            const auto v = values[elem.variable_index];

            if (v != unassigned) {
                rhs -= static_cast<long long>(elem.factor) * v;
            } else {
                bx_assert(free_variable == -1);
                free_variable = elem.variable_index;
                free_factor = elem.factor;
            }
        }

        const bool can_be_0 = is_satisfied(kind, 0, rhs);

        if (free_variable < 0) {
            if (!can_be_0)
                unrealisable(kind, index);
            return;
        }

        const bool can_be_1 = is_satisfied(kind, free_factor, rhs);

        if (!can_be_0 && !can_be_1)
            unrealisable(kind, index);

        if (can_be_0 != can_be_1 && !assign(free_variable, can_be_1))
            unrealisable(kind, index);

        debug(ctx,
              "      - {} constraint {} removed.\n",
              constraint_kind_names[to_index(kind)],
              cst.label);
    }

    // Compares the right-hand side with the smallest and largest reachable
    // left-hand side: a constraint always satisfied is removed, one satisfied
    // only at an extremum fixes all its variables.
    void tighten(constraint_kind kind, int index)
    {
        const auto& cst = get(kind, index);

        long long rhs = cst.value;
        long long lo = 0;
        long long hi = 0;
        int free_elements = 0;

        for (const auto& elem : cst.elements) {
            const auto v = values[elem.variable_index];

            if (v != unassigned) {
                rhs -= static_cast<long long>(elem.factor) * v;
            } else {
                ++free_elements;
                (elem.factor < 0 ? lo : hi) += elem.factor;
            }
        }

        if (free_elements <= 1) {
            decide(kind, index);
            return;
        }

        const bool reach_lo = is_satisfied(kind, lo, rhs);
        const bool reach_hi = is_satisfied(kind, hi, rhs);

        if (kind == constraint_kind::equal ? (rhs < lo || hi < rhs)
                                           : (!reach_lo && !reach_hi))
            unrealisable(kind, index);

        const bool always =
          kind == constraint_kind::equal ? lo == hi : (reach_lo && reach_hi);

        if (always) {
            remaining[to_index(kind)][index] = 0;

            debug(ctx,
                  "      - {} constraint {} removed.\n",
                  constraint_kind_names[to_index(kind)],
                  cst.label);
            return;
        }

        // A value of a variable that pushes the reachable left-hand side
        // beyond the right-hand side is excluded, so the variable takes the
        // other one (activity-based fixing). The decided values follow from
        // the same bounds, so they can be applied in any order.
        const bool upper = kind != constraint_kind::greater;
        const bool lower = kind != constraint_kind::less;

        for (const auto& elem : cst.elements) {
            if (values[elem.variable_index] != unassigned || elem.factor == 0)
                continue;

            const long long a = elem.factor;
            const bool excludes_1 = (upper && a > 0 && lo + a > rhs) ||
                                    (lower && a < 0 && hi + a < rhs);
            const bool excludes_0 = (upper && a < 0 && lo - a > rhs) ||
                                    (lower && a > 0 && hi - a < rhs);

            if (excludes_0 && excludes_1)
                unrealisable(kind, index);

            if ((excludes_0 || excludes_1) &&
                !assign(elem.variable_index, excludes_0))
                unrealisable(kind, index);
        }
    }

    void propagate()
    {
        while (!to_propagate.empty()) {
            const int variable = to_propagate.back();
            to_propagate.pop_back();

            for (auto kind : constraint_kinds) {
                auto& counters = remaining[to_index(kind)];

                for (int cst : cache[to_index(kind)](variable)) {
                    if (counters[cst] <= 0)
                        continue;

                    if (--counters[cst] <= 1) {
                        decide(kind, cst);
                    } else if (!is_touched[to_index(kind)][cst]) {
                        is_touched[to_index(kind)][cst] = true;
                        touched[to_index(kind)].emplace_back(cst);
                    }
                }
            }
        }
    }

    void initialize()
    {
        values.assign(pb.vars.values.size(), unassigned);
        assigned = 0;
        to_propagate.clear();

        for (auto kind : constraint_kinds) {
            const auto& csts = *constraints[to_index(kind)];
            auto& counters = remaining[to_index(kind)];

            counters.resize(csts.size());
            for (std::size_t i = 0, e = csts.size(); i != e; ++i)
                counters[i] = bx::length(csts[i].elements);
        }

        for (int i = 0, e = bx::length(pb.vars.values); i != e; ++i) {
            const auto lower = std::max(pb.vars.values[i].min, 0);
            const auto upper = std::min(pb.vars.values[i].max, 1);

            if (lower > upper)
                unrealisable(i);

            if (lower == upper && !assign(i, lower == 1))
                unrealisable(i);
        }

        for (auto kind : constraint_kinds) {
            auto& list = touched[to_index(kind)];
            const auto size = constraints[to_index(kind)]->size();

            is_touched[to_index(kind)].assign(size, true);
            list.resize(size);
            std::iota(list.begin(), list.end(), 0);
        }

        reduce();
    }

    // Tightens again every constraint that lost a variable since its last
    // check, until no new assignment appears.
    void reduce()
    {
        propagate();

        std::vector<int> list;
        bool again = true;

        while (again) {
            again = false;

            for (auto kind : constraint_kinds) {
                list.clear();
                list.swap(touched[to_index(kind)]);

                for (int cst : list) {
                    is_touched[to_index(kind)][cst] = false;

                    if (remaining[to_index(kind)][cst] > 0)
                        tighten(kind, cst);
                }

                again = again || !to_propagate.empty();
                propagate();
            }
        }
    }

    // Counts the remaining constraints that a variable going up (or down)
    // could violate.
    std::pair<int, int> locks(int variable) const noexcept
    {
        int up = 0;
        int down = 0;

        for (auto kind : constraint_kinds) {
            for (int cst : cache[to_index(kind)](variable)) {
                if (remaining[to_index(kind)][cst] <= 0)
                    continue;

                int factor = 0;
                for (const auto& elem : get(kind, cst).elements)
                    if (elem.variable_index == variable)
                        factor += elem.factor;

                if (factor == 0)
                    continue;

                const bool positive = factor > 0;
                if (kind != constraint_kind::greater)
                    ++(positive ? up : down);
                if (kind != constraint_kind::less)
                    ++(positive ? down : up);
            }
        }

        return { up, down };
    }

    // Dual fixing: a variable that no remaining constraint prevents from
    // going down (up) and whose cost does not favour the other value keeps an
    // optimal solution at 0 (1). Variables of quadratic terms are kept since
    // their best value depends on other variables.
    bool dual_fixing()
    {
        std::vector<double> linear(values.size(), 0.0);
        std::vector<char> quadratic(values.size(), false);

        const double sense =
          pb.type == bx::objective_function_type::maximize ? -1.0 : 1.0;

        for (const auto& elem : pb.objective.elements)
            linear[elem.variable_index] += sense * elem.factor;

        for (const auto& elem : pb.objective.qelements) {
            quadratic[elem.variable_index_a] = true;
            quadratic[elem.variable_index_b] = true;
        }

        bool fixed = false;

        for (int i = 0, e = bx::length(values); i != e; ++i) {
            if (values[i] != unassigned || quadratic[i])
                continue;

            const auto [up, down] = locks(i);

            if (down == 0 && linear[i] >= 0.0)
                fixed = assign(i, false) || fixed;
            else if (up == 0 && linear[i] <= 0.0)
                fixed = assign(i, true) || fixed;
        }

        return fixed;
    }

    // Alternates propagation and dual fixing until neither assigns anything.
    void simplify()
    {
        reduce();

        while (dual_fixing())
            reduce();
    }

public:
    explicit preprocessor(const bx::context& ctx_, const Problem& pb_)
      : ctx(ctx_)
      , pb(pb_)
      , constraints{ &pb.equal_constraints,
                     &pb.greater_constraints,
                     &pb.less_constraints }
      , cache{
          variable_constraints(pb.equal_constraints, pb.vars.values.size()),
          variable_constraints(pb.greater_constraints, pb.vars.values.size()),
          variable_constraints(pb.less_constraints, pb.vars.values.size())
      }
    {
    }

    bx::problem operator()(int variable_index, bool variable_value)
    {
        initialize();

        if (!assign(variable_index, variable_value))
            unrealisable(variable_index);

        simplify();

        info(ctx,
             "  - Preprocessor finish. Removed variables {} (size: {})\n",
             assigned,
             to_string(bx::memory_consumed_size(memory_consumed(pb))));

        return make_problem();
    }

    bx::problem operator()(std::span<const std::pair<int, bool>> assignments)
    {
        initialize();

        for (const auto& [variable, value] : assignments)
            if (!assign(variable, value))
                unrealisable(variable);

        simplify();

        return make_problem();
    }

    bx::problem operator()()
    {
        initialize();
        simplify();

        info(ctx,
             "  - Preprocessor finish. Removed variables {} (size: {})\n",
             assigned,
             to_string(bx::memory_consumed_size(memory_consumed(pb))));

        return make_problem();
    }

private:
    auto make_problem() const -> bx::problem
    {
        bx::problem copy;

        copy.strings = pb.strings;
        copy.type = pb.type;

        // Maps each origin variable to its index in the copy, or to -1 and its
        // assigned value.

        std::vector<std::pair<int, bool>> mapping(pb.vars.values.size(),
                                                  { -1, false });

        int c = 0;
        for (int i = 0, e = bx::length(pb.vars.values); i != e; ++i)
            mapping[i] = values[i] == unassigned
                           ? std::make_pair(c++, false)
                           : std::make_pair(-1, values[i] == 1);

        bx_assert(c == bx::length(pb.vars.values) - assigned);

        copy.objective = objective_function_exclude_copy(mapping);
        std::tie(copy.vars, copy.affected_vars) = variables_exclude_copy();

        for (auto kind : constraint_kinds) {
            auto& out = kind == constraint_kind::equal ? copy.equal_constraints
                        : kind == constraint_kind::greater
                          ? copy.greater_constraints
                          : copy.less_constraints;

            constraints_exclude_copy(mapping,
                                     remaining[to_index(kind)],
                                     *constraints[to_index(kind)],
                                     out);
        }

        bx_ensures(is_valid_variable_index(copy.equal_constraints,
                                           bx::length(copy.vars.values)) &&
                   is_valid_variable_index(copy.greater_constraints,
                                           bx::length(copy.vars.values)) &&
                   is_valid_variable_index(copy.less_constraints,
                                           bx::length(copy.vars.values)));

        copy.problem_type = copy.which_problem_type();

        info(ctx,
             "- Preprocessing finishes (size: {})\n",
             to_string(bx::memory_consumed_size(memory_consumed(copy))));

        if (ctx.parameters.debug) {
            if (std::FILE* ofs = std::fopen("preprocessed.lp", "w"); ofs) {
                info(ctx, "- Write file preprocessed.lp\n");
                fmt::print(ofs, "{}", copy);
                std::fclose(ofs);
            }
        }

        return copy;
    }

    auto objective_function_exclude_copy(
      const std::vector<std::pair<int, bool>>& mapping) const
      -> bx::objective_function
    {
        bx::objective_function ret;

        ret.value = pb.objective.value;

        for (const auto& elem : pb.objective.elements) {
            const auto [index, value] = mapping[elem.variable_index];

            if (index >= 0)
                ret.elements.emplace_back(elem.factor, index);
            else if (value)
                ret.value += elem.factor;
        }

        for (const auto& elem : pb.objective.qelements) {
            const auto [a, a_value] = mapping[elem.variable_index_a];
            const auto [b, b_value] = mapping[elem.variable_index_b];

            if (a >= 0 && b >= 0)
                ret.qelements.emplace_back(elem.factor, a, b);
            else if (a >= 0 && b_value)
                ret.elements.emplace_back(elem.factor, a);
            else if (b >= 0 && a_value)
                ret.elements.emplace_back(elem.factor, b);
            else if (a < 0 && b < 0 && a_value && b_value)
                ret.value += elem.factor;
        }

        return ret;
    }

    static auto affected_variable_copy(const bx::raw_problem& /*pb*/,
                                       size_t size) -> bx::affected_variables
    {
        bx::affected_variables ret;
        ret.names.reserve(size);
        ret.values.reserve(size);

        return ret;
    }

    static auto affected_variable_copy(const bx::problem& pb, size_t size)
      -> bx::affected_variables
    {
        bx::affected_variables ret;
        ret.names.reserve(pb.affected_vars.names.size() + size);
        ret.values.reserve(pb.affected_vars.values.size() + size);

        std::copy(pb.affected_vars.names.cbegin(),
                  pb.affected_vars.names.end(),
                  std::back_inserter(ret.names));

        std::copy(pb.affected_vars.values.cbegin(),
                  pb.affected_vars.values.end(),
                  std::back_inserter(ret.values));

        return ret;
    }

    auto variables_exclude_copy() const
      -> std::tuple<bx::variables, bx::affected_variables>
    {
        auto ret_aff_vars =
          affected_variable_copy(pb, static_cast<std::size_t>(assigned));

        bx::variables ret_vars;

        ret_vars.names.reserve(pb.vars.names.size() - assigned);
        ret_vars.values.reserve(pb.vars.names.size() - assigned);

        for (std::size_t i = 0, e = pb.vars.values.size(); i != e; ++i) {
            if (values[i] == unassigned) {
                ret_vars.names.emplace_back(pb.vars.names[i]);
                ret_vars.values.emplace_back(pb.vars.values[i]);
            } else {
                ret_aff_vars.names.emplace_back(pb.vars.names[i]);
                ret_aff_vars.values.emplace_back(values[i]);
            }
        }

        return std::make_tuple(ret_vars, ret_aff_vars);
    }

    // Decided constraints (counter at 0) are already satisfied by the
    // assignments and are not copied.
    static void constraints_exclude_copy(
      const std::vector<std::pair<int, bool>>& mapping,
      const std::vector<int>& constraints_size,
      const std::vector<bx::constraint>& constraints,
      std::vector<bx::constraint>& copy)
    {
        for (int i = 0, e = bx::length(constraints); i != e; ++i) {
            if (constraints_size[i] <= 0)
                continue;

            auto& cst = copy.emplace_back();
            cst.id = constraints[i].id;
            cst.label = constraints[i].label;
            cst.value = constraints[i].value;

            for (const auto& elem : constraints[i].elements) {
                const auto [index, value] = mapping[elem.variable_index];

                if (index >= 0)
                    cst.elements.emplace_back(elem.factor, index);
                else if (value)
                    cst.value -= elem.factor;
            }
        }
    }
};

// The solvers only handle 0-1 variables: an integer variable above 1 would
// silently be solved as binary, while undeclared variables, real in the LP
// format, have long been read as binary.
static void
check_variable_domains(const bx::context& ctx, const bx::raw_problem& pb)
{
    int real = 0;

    for (int i = 0, e = bx::length(pb.vars.values); i != e; ++i) {
        const auto& v = pb.vars.values[i];

        if (v.type == bx::variable_type::real) {
            ++real;
        } else if (v.max > 1) {
            error(ctx,
                  "  - Variable {} has bounds [{}, {}]: only 0-1 variables "
                  "are supported.\n",
                  pb.vars.names[i],
                  v.min,
                  v.max);

            throw bx::solver_failure(
              bx::solver_error_tag::not_binary_variable);
        }
    }

    if (real > 0)
        warning(ctx,
                "  - {} real variables are solved as 0-1 variables.\n",
                real);
}

} // namespace anonymous

namespace baryonyx {

std::tuple<problem, problem>
split(const context& ctx, const problem& pb, int variable_index_to_affect)
{
    bx_expects(variable_index_to_affect >= 0 &&
               variable_index_to_affect < length(pb.vars.values));

    info(ctx,
         "  - Preprocessor starts split of variable {} (size: {})\n",
         pb.vars.names[variable_index_to_affect],
         to_string(memory_consumed_size(memory_consumed(pb))));

    ::preprocessor<problem> pp(ctx, pb);

    return std::make_tuple(pp(variable_index_to_affect, true),
                           pp(variable_index_to_affect, false));
}

problem
affect(const context& ctx,
       const problem& pb,
       int variable_index,
       bool variable_value)
{
    bx_expects(variable_index >= 0 && variable_index < length(pb.vars.values));

    info(
      ctx,
      "  - Preprocessor starts affectation of variable {} to {} (size: {})\n",
      pb.vars.names[variable_index],
      variable_value,
      to_string(memory_consumed_size(memory_consumed(pb))));

    ::preprocessor<problem> pp(ctx, pb);

    return pp(variable_index, variable_value);
}

struct affecter::impl
{
    impl(const context& ctx, const problem& pb)
      : pp(ctx, pb)
    {
    }

    ::preprocessor<problem> pp;
};

affecter::affecter(const context& ctx, const problem& pb)
  : m_impl(std::make_unique<impl>(ctx, pb))
{
}

affecter::~affecter() noexcept = default;

problem
affecter::operator()(std::span<const std::pair<int, bool>> assignments)
{
    return m_impl->pp(assignments);
}

problem
affect(const context& ctx,
       const problem& pb,
       std::span<const std::pair<int, bool>> assignments)
{
    return affecter(ctx, pb)(assignments);
}

problem
preprocess(const context& ctx, const raw_problem& raw_pb)
{
    info(ctx,
         "- Preprocessing starts (size: {})\n",
         to_string(memory_consumed_size(memory_consumed(raw_pb))));

    check_variable_domains(ctx, raw_pb);

    ::preprocessor<raw_problem> pp(ctx, raw_pb);

    return pp();
}

problem
unpreprocess(const context& ctx, const raw_problem& raw_pb)
{
    info(ctx,
         "- Unprepossessing starts (size: {})\n",
         to_string(memory_consumed_size(memory_consumed(raw_pb))));

    check_variable_domains(ctx, raw_pb);

    return problem(raw_pb);
}

} // namespace baryonyx
