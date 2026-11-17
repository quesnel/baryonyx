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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_PRODUCTS_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_PRODUCTS_HPP

#include <baryonyx/core>

#include <optional>
#include <string_view>
#include <vector>

namespace baryonyx {

/**
 * Variables z = x.y written with linear constraints (z <= x, z <= y,
 * z >= x + y - 1). @c reduce removes them with their constraints and moves
 * their cost to a quadratic term of the objective; @c restore puts them back
 * in the result with the value x.y.
 */
class linearized_products
{
public:
    /** The reduced problem, or nothing if @c pb has no such variable. */
    std::optional<raw_problem> reduce(const context& ctx,
                                      const raw_problem& pb);

    void restore(result& r) const;

    int size() const noexcept
    {
        return static_cast<int>(m_products.size());
    }

private:
    struct product
    {
        std::string_view name;
        std::string_view a;
        std::string_view b;
    };

    std::vector<product> m_products;
};

} // namespace baryonyx

#endif
