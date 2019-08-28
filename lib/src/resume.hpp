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
#include <chrono>
#include <iomanip>
#include <iterator>
#include <numeric>
#include <ostream>
#include <string_view>

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

inline std::ostream&
operator<<(std::ostream& os, const resume& res)
{
    auto store = os.flags();

    os << std::setprecision(std::numeric_limits<double>::digits10 + 1);

    if (res.use_lp_format) {
        os << "\\ Problem statistics:\n"
           << R"(\  type: ) undefined\n'
           << R"(\  nb variables: )"
           << std::accumulate(res.variables.begin(), res.variables.end(), 0)
           << '\n'
           << R"(\   ..... real: )" << res.variables[0] << '\n'
           << R"(\   ... binary: )" << res.variables[1] << '\n'
           << R"(\   .. general: )" << res.variables[2] << '\n'
           << R"(\  nb constraints: )"
           << std::accumulate(res.constraints.begin(), res.constraints.end(), 0)
           << '\n'
           << R"(\   ........ =  : )" << res.constraints[0] << '\n'
           << R"(\   ........ >= : )" << res.constraints[1] << '\n'
           << R"(\   ........ <= : )" << res.constraints[2] << '\n'
           << R"(\  minimal value.: )" << std::get<0>(res.minmax) << '\n'
           << R"(\  maximal value.: )" << std::get<1>(res.minmax) << '\n';
    } else {
        os << "Problem statistics:\n"
           << "  * type: undefined\n"
           << "  * variables: "
           << std::accumulate(res.variables.begin(), res.variables.end(), 0)
           << '\n'
           << "    - real: " << res.variables[0] << '\n'
           << "    - binary: " << res.variables[1] << '\n'
           << "    - general: " << res.variables[2] << '\n'
           << "  * constraints: "
           << std::accumulate(res.constraints.begin(), res.constraints.end(), 0)
           << '\n'
           << "    - constraint =  : " << res.constraints[0] << '\n'
           << "    - constraint >= : " << res.constraints[1] << '\n'
           << "    - constraint <= : " << res.constraints[2] << '\n'
           << "  * objective:\n"
           << "    - minimal value.: " << std::get<0>(res.minmax) << '\n'
           << "    - maximal value.: " << std::get<1>(res.minmax) << '\n';
    }

    os.flags(store);

    return os;
}

inline std::ostream&
operator<<(std::ostream& os, const affected_variables& var)
{
    std::size_t i = 0, e = var.names.size();

    for (; i != e; ++i)
        os << var.names[i] << ": " << (var.values[i] ? 1 : 0) << '\n';

    return os;
}

} // namespace baryonyx

#endif
