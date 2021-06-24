/* Copyright (C) 2018-2021 INRAE
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

#include <string>
#include <tuple>

namespace baryonyx {

struct raw_problem;
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

    const auto kb = static_cast<double>(size) / 1024.;
    const auto mb = kb / 1024.;
    const auto gb = mb / 1024.;

    if (gb > 0.5)
        return { gb, show_size_type::GB };

    if (mb > 0.5)
        return { mb, show_size_type::MB };

    if (kb > 0.5)
        return { kb, show_size_type::KB };

    return { static_cast<double>(size), show_size_type::B };
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
memory_consumed(const raw_problem& pb) noexcept;

/**
 * @brief Calculates an estimate of the consumed memory by the data structures.
 *
 * @details Return the size in bytes of the memory used by the @c
 *     baryonyx::problem data structures.
 */
std::size_t
memory_consumed(const problem& pb) noexcept;

} // namespace baryonyx

#endif
