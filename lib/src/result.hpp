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

#ifndef ORG_VLEPROJECT_BARYONYX_LIB_PRIVATE_RESULT_HPP
#define ORG_VLEPROJECT_BARYONYX_LIB_PRIVATE_RESULT_HPP

#include <baryonyx/core>

#include "itm-common.hpp"

#include <iosfwd>

namespace baryonyx {

template<typename Mode>
struct raw_result
{
    raw_result() = default;

    raw_result(raw_result&&) = default;
    raw_result& operator=(raw_result&&) = default;

    raw_result(int variables)
      : x(variables)
    {}

    raw_result(const bit_array& x_,
               double value_,
               double duration_,
               long int loop_)
      : x(x_)
      , value(value_)
      , duration(duration_)
      , loop(loop_)
      , remaining_constraints(0)
    {}

    raw_result(const bit_array& x_,
               index remaining_constraints_,
               double duration_,
               long int loop_)
      : x(x_)
      , value(itm::bad_value<Mode, double>())
      , duration(duration_)
      , loop(loop_)
      , remaining_constraints(remaining_constraints_)
    {}

    bit_array x;
    double value = itm::bad_value<Mode, double>();
    double duration = 0.0;
    std::size_t hash = 0;
    long int loop = 0;
    index remaining_constraints = std::numeric_limits<index>::max();

    void make_hash() noexcept
    {
        hash = bit_array_hash()(x);
    }

    bool is_solution() const noexcept
    {
        return remaining_constraints == 0;
    }

    void swap(raw_result& other) noexcept
    {
        x.swap(other.x);
        std::swap(value, other.value);
        std::swap(duration, other.duration);
        std::swap(loop, other.loop);
        std::swap(remaining_constraints, other.remaining_constraints);
    }
};

template<typename Mode>
struct raw_result_compare
{
    using is_transparent = std::true_type;

    bool operator()(const raw_result<Mode>& lhs,
                    const raw_result<Mode>& rhs) const noexcept
    {
        if (lhs.remaining_constraints == 0 && rhs.remaining_constraints == 0)
            return itm::is_better_solution<Mode>(lhs.value, rhs.value);

        return lhs.remaining_constraints < rhs.remaining_constraints;
    }

    bool operator()(const raw_result<Mode>& lhs, double value) const noexcept
    {
        return itm::is_better_solution<Mode>(lhs.value, value);
    }

    bool operator()(double value, const raw_result<Mode>& rhs) const noexcept
    {
        return itm::is_better_solution<Mode>(value, rhs.value);
    }
};

template<typename Mode>
bool
operator==(const raw_result<Mode>& lhs, const raw_result<Mode>& rhs) noexcept
{
    return lhs.remaining_constraints == rhs.remaining_constraints &&
           lhs.x == rhs.x;
}

template<typename Mode>
void
convert(const raw_result<Mode>& source, solution& sol, int variables)
{
    sol.value = source.value;
    sol.variables.resize(variables);

    for (int i = 0; i != variables; ++i)
        sol.variables[i] = static_cast<var_value>(source.x[i]);
}

std::istream&
operator>>(std::istream& is, result& r);

struct best_solution_writer
{
    const result& res;

    best_solution_writer(const result& res_)
      : res(res_)
    {}
};

std::ostream&
operator<<(std::ostream& os, const best_solution_writer& writer);

} // namespace baryonyx

template<>
struct fmt::formatter<baryonyx::result>
{
    constexpr auto parse(format_parse_context& ctx)
    {
        return ctx.begin();
    }

    template<typename FormatContext>
    auto format(const baryonyx::result& result, FormatContext& ctx)
    {
        fmt::format_to(ctx.out(),
                            "\\ solver................: {}\n"
                            "\\ constraints...........: {}\n"
                            "\\ variables.............: {}\n"
                            "\\ duration..............: {}s\n"
                            "\\ loop..................: {}\n"
                            "\\ status................: {}\n",
                            result.method,
                            result.constraints,
                            result.variables,
                            result.duration,
                            result.loop,
                            result.status);

        if (result.status == baryonyx::result_status::success) {
            if (!result.solutions.empty()) {
                fmt::format_to(ctx.out(),
                               "\\ value.................: {:.10g}\n",
                               result.solutions.back().value);

                if (result.solutions.size() > 1) {
                    fmt::format_to(ctx.out(), "\\ other value...........: ");
                    for (const auto& elem : result.solutions)
                        fmt::format_to(ctx.out(), "{:.10g} ", elem.value);
                    fmt::format_to(ctx.out(), "\n");
                }

                fmt::format_to(ctx.out(), "\\ variables.............: \n");

                for (std::size_t i = 0, e = result.affected_vars.names.size();
                     i != e;
                     ++i)
                    fmt::format_to(ctx.out(),
                                   "{}={}\n",
                                   result.affected_vars.names[i],
                                   result.affected_vars.values[i] ? 1 : 0);

                for (std::size_t i = 0, e = result.variable_name.size();
                     i != e;
                     ++i)
                    fmt::format_to(ctx.out(),
                                   "{}={}\n",
                                   result.variable_name[i],
                                   result.solutions.back().variables[i] ? 1
                                                                        : 0);
            }
        } else if (result.status ==
                     baryonyx::result_status::time_limit_reached ||
                   result.status ==
                     baryonyx::result_status::kappa_max_reached ||
                   result.status == baryonyx::result_status::limit_reached) {
            fmt::format_to(ctx.out(),
                           "\\ remaining constraints.: {}\n",
                           result.remaining_constraints);
        }

        return ctx.out();
    }
};

#endif
