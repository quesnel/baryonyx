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

#ifndef ORG_VLEPROJECT_BARYONYX_PRIVATE_MEMORY_HPP
#define ORG_VLEPROJECT_BARYONYX_PRIVATE_MEMORY_HPP

#include <baryonyx/core>

#include <tuple>

namespace baryonyx {

struct problem;

enum class show_size_type
{
    GB = 0,
    MB,
    KB,
    B
};

template<typename T>
inline std::tuple<double, show_size_type>
memory_consumed_size(T size) noexcept
{
    static_assert(
      (std::is_integral<T>::value || std::is_floating_point<T>::value),
      "Integer or real required.");

    auto kb = static_cast<double>(size) / 1024.;
    auto mb = kb / 1024.;
    auto gb = mb / 1024.;

    auto ret = std::make_tuple(static_cast<double>(size), show_size_type::B);

    if (gb > 0.5) {
        std::get<0>(ret) = gb;
        std::get<1>(ret) = show_size_type::GB;
    } else if (mb > 0.5) {
        std::get<0>(ret) = mb;
        std::get<1>(ret) = show_size_type::MB;
    } else if (kb > 0.5) {
        std::get<0>(ret) = kb;
        std::get<1>(ret) = show_size_type::KB;
    }

    return ret;
}

std::string
to_string(std::tuple<double, show_size_type> memory_used);

/**
 * @brief Calculates an estimate of the consumed memory by the data structures.
 *
 * @details Return the size in bytes of the memory used by the @c
 *     baryonyx::problem data structures.
 */
std::size_t
memory_consumed(const baryonyx::raw_problem& pb) noexcept;

/**
 * @brief Calculates an estimate of the consumed memory by the data structures.
 *
 * @details Return the size in bytes of the memory used by the @c
 *     baryonyx::problem data structures.
 */
std::size_t
memory_consumed(const baryonyx::problem& pb) noexcept;

} // namespace baryonyx

#endif
