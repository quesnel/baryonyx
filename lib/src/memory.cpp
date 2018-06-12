/* Copyright (C) 2018 INRA
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

#include "memory.hpp"

#include <numeric>

#include <fmt/format.h>

static std::size_t
memory_consumed(const baryonyx::variables& vv) noexcept
{
    auto ret = sizeof(baryonyx::variables);

    for (const auto& str : vv.names)
        ret += str.capacity();

    return ret + vv.values.capacity() * sizeof(baryonyx::variable_value);
}

static std::size_t
memory_consumed(const baryonyx::constraint& cst) noexcept
{
    auto ret = sizeof(baryonyx::constraint);

    return ret + cst.elements.capacity() * sizeof(baryonyx::function_element);
}

static std::size_t
memory_consumed(const baryonyx::objective_function& obj) noexcept
{
    auto ret = sizeof(baryonyx::objective_function);

    return ret + obj.elements.capacity() *
                   sizeof(baryonyx::objective_function_element);
}

static std::size_t
memory_consumed(const baryonyx::affected_variables& av) noexcept
{
    auto ret = sizeof(baryonyx::affected_variables);

    ret = std::accumulate(av.names.cbegin(),
                          av.names.cend(),
                          ret,
                          [](std::size_t size, const std::string& str) {
                              return str.capacity() + size;
                          });

    return ret + av.values.capacity();
}

static std::size_t
memory_consumed(const std::vector<baryonyx::constraint>& csts) noexcept
{
    return std::accumulate(
      csts.cbegin(), csts.cend(), 0, [](std::size_t size, const auto& elem) {
          return size + memory_consumed(elem);
      });
}

namespace baryonyx {

std::string
to_string(std::tuple<double, show_size_type> mem)
{
    switch (std::get<1>(mem)) {
    case show_size_type::GB:
        return fmt::format("{:.2f} GB", std::get<0>(mem));
    case show_size_type::MB:
        return fmt::format("{:.2f} MB", std::get<0>(mem));
    case show_size_type::KB:
        return fmt::format("{:.2f} KB", std::get<0>(mem));
    default:
        return fmt::format("{:.2f} B", std::get<0>(mem));
    }
}

std::size_t
memory_consumed(const problem& pb) noexcept
{
    auto ret = sizeof(problem);

    ret += ::memory_consumed(pb.objective);
    ret += ::memory_consumed(pb.equal_constraints);
    ret += ::memory_consumed(pb.greater_constraints);
    ret += ::memory_consumed(pb.less_constraints);
    ret += ::memory_consumed(pb.vars);
    ret += ::memory_consumed(pb.affected_vars);

    return ret;
}

} // namespace baryonyx
