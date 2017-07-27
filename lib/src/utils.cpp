/* Copyright (C) 2017 INRA
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

#include "utils.hpp"

namespace baryonyx {

std::string
stringf(const char* format, ...) noexcept
{
    try {
        int n;
        int size = 256;
        std::string ret(size, '\0');

        for (;;) {
            va_list ap;
            va_start(ap, format);
            n = vsnprintf(&ret[0], size, format, ap);
            va_end(ap);

            if (n < 0)
                return {};

            if (n < size) {
                ret.erase(n, std::string::npos);
                return ret;
            }

            size = n + 1;
            ret.resize(size);
        }
    } catch (const std::bad_alloc& e) {
        fputs("baryonyx::stringf: not enough memory\n", stderr);
        return {};
    }
}

std::string
vstringf(const char* format, va_list ap) noexcept
{
    try {
        int n;
        int size = 256;
        std::string ret(size, '\0');

        for (;;) {
            n = vsnprintf(&ret[0], size, format, ap);

            if (n < 0)
                return {};

            if (n < size) {
                ret.erase(n, std::string::npos);
                return ret;
            }

            size = n + 1;
            ret.resize(size);
        }
    } catch (const std::bad_alloc& e) {
        fputs("baryonyx::vstringf: not enough memory\n", stderr);
        return {};
    }
}

} // namespace baryonyx
