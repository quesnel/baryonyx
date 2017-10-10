/* Copyright (C) 2017 INRA
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublipnse, and/or sell copies of the Software, and to
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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_KNAPSACK_SOLVER_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_KNAPSACK_SOLVER_HPP

#include "fixed_2darray.hpp"

#include <vector>

#include <cassert>

namespace baryonyx {

template <typename T>
std::vector<bool>
knapsack_solver_dp(const std::vector<T>& values,
                   const std::vector<T>& weights,
                   T capacity)
{
    assert(values.size() == weights.size() &&
           "values and weights must be compatible");

    fixed_2darray<T> best(values.size() + 1, capacity + 1, 0);

    for (std::size_t i{ 1 }, ei{ values.size() + 1 }; i != ei; ++i) {
        for (std::size_t j{ 1 }, ej{ static_cast<std::size_t>(capacity + 1) };
             j != ej;
             ++j) {
            if ((std::size_t)weights[i - 1] > j) {
                best(i, j) = best(i - 1, j);
            } else {
                auto c1 = best(i - 1, j);
                auto c2 = best(i - 1, j - weights[i - 1]) + values[i - 1];
                best(i, j) = std::max(c1, c2);
            }
        }
    }

    std::size_t i{ values.size() };
    std::size_t j{ static_cast<std::size_t>(capacity) };
    std::vector<bool> x(i, false);
    int id{ 0 };

    for (; i > 0; --i, ++id) {
        if (best(i, j) != best(i - 1, j)) {
            x[id] = true;
            j -= weights[i - 1];
        }
    }

    return x;
}
}

#endif
