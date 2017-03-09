/* Copyright (C) 2016 INRA
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

#ifndef ORG_VLEPROJECT_LP_UTILS_HPP
#define ORG_VLEPROJECT_LP_UTILS_HPP

#include <cmath>
#include <lpcore>

namespace lp {

inline void
Expects(bool condition)
{
    if (not condition)
        throw precondition_error("precondition failure");
}

inline void
Ensures(bool condition)
{
    if (not condition)
        throw postcondition_error("postcondition failure");
}

inline void
Expects(bool condition, const char* s)
{
    if (not condition)
        throw precondition_error(s);
}

inline void
Ensures(bool condition, const char* s)
{
    if (not condition)
        throw postcondition_error(s);
}

template <typename T>
inline bool
is_essentially_equal(const T v1, const T v2, const T epsilon)
{
    static_assert(std::is_floating_point<T>::value);

    return fabs((v1) - (v2)) <=
           ((fabs(v1) > fabs(v2) ? fabs(v2) : fabs(v1)) * (epsilon));
}
}
#endif
