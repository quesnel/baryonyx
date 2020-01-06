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

#ifndef ORG_VLEPROJECT_BARYONYX_LIB_PRIVATE_RESUME_HPP
#define ORG_VLEPROJECT_BARYONYX_LIB_PRIVATE_RESUME_HPP

#include <array>
#include <numeric>
#include <string_view>

#include <fmt/format.h>

#include <baryonyx/core-utils>
#include <baryonyx/core>

#include "problem.hpp"

namespace baryonyx {

struct resume
{
    resume(const raw_problem& pb, bool use_lp_format_ = true)
      : variables({})
      , constraints({})
      , minmax(compute_min_max_objective_function(pb.objective))
      , use_lp_format(use_lp_format_)
    {
        variables = std::accumulate(pb.vars.values.begin(),
                                    pb.vars.values.end(),
                                    variables,
                                    [](std::array<int, 3>& value, auto vv) {
                                        switch (vv.type) {
                                        case variable_type::real:
                                            value[0]++;
                                            break;
                                        case variable_type::binary:
                                            value[1]++;
                                            break;
                                        case variable_type::general:
                                            value[2]++;
                                            break;
                                        }

                                        return value;
                                    });

        constraints[0] = static_cast<int>(pb.equal_constraints.size());
        constraints[1] = static_cast<int>(pb.greater_constraints.size());
        constraints[2] = static_cast<int>(pb.less_constraints.size());
    }

    std::array<int, 3> variables;
    std::array<int, 3> constraints;
    std::tuple<double, double> minmax;
    bool use_lp_format;
};

} // namespace baryonyx

template<>
struct fmt::formatter<baryonyx::resume>
{
    constexpr auto parse(format_parse_context& ctx)
    {
        return ctx.begin();
    }

    template<typename FormatContext>
    auto format(const baryonyx::resume& res, FormatContext& ctx)
    {
        const auto variables{ std::accumulate(
          res.variables.begin(), res.variables.end(), 0) };
        const auto constraints{ std::accumulate(
          res.constraints.begin(), res.constraints.end(), 0) };

        return format_to(ctx.out(),
                         "Problem statistics:\n"
                         "  * type: undefined\n"
                         "  * variables: {}\n"
                         "    - real: {}\n"
                         "    - binary: {}\n"
                         "    - general: {}\n"
                         "  * constraints: {}\n"
                         "    - constraint ==: {}\n"
                         "    - constraint >=: {}\n"
                         "    - constraint <=: {}\n"
                         "  * objective:\n"
                         "    - minimal value.: {:.10g}\n"
                         "    - maximal value.: {:.10g}\n",
                         variables,
                         res.variables[0],
                         res.variables[1],
                         res.variables[2],
                         constraints,
                         res.constraints[0],
                         res.constraints[1],
                         res.constraints[2],
                         std::get<0>(res.minmax),
                         std::get<1>(res.minmax));
    }
};

#endif
