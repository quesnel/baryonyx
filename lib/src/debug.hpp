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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_DEBUG_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_DEBUG_HPP

#include "private.hpp"

#include <exception>
#include <source_location>
#include <stdexcept>

#include <cstdio>

//
// In internal API, we prefer abort application when precondition,
// postcondition or assertion fail. For external API, user can select to use
// exception or return with message in terminal.
//

namespace baryonyx::details {

[[noreturn]] inline void
fail_fast(const char* type,
          const char* cond,
          const std::source_location location =
            std::source_location::current())
{
    std::fprintf(stderr,
                 "%s [%s] failure at %s:%u in %s\n",
                 type,
                 cond,
                 location.file_name(),
                 static_cast<unsigned>(location.line()),
                 location.function_name());

#ifndef BARYONYX_FAIL_FAST
    std::terminate();
#else
    throw std::logic_error(cond);
#endif
}

} // namespace baryonyx::details

#define BX_CONTRACT_CHECK(type, cond)                                         \
    do {                                                                      \
        if (!(cond)) [[unlikely]]                                             \
            baryonyx::details::fail_fast(type, #cond);                        \
    } while (false)

#ifdef BARYONYX_FULL_OPTIMIZATION
#define bx_expects(cond)
#define bx_ensures(cond)
#define bx_assert(cond)
#else
#define bx_expects(cond) BX_CONTRACT_CHECK("Precondition", cond)
#define bx_ensures(cond) BX_CONTRACT_CHECK("Postcondition", cond)
#define bx_assert(cond) BX_CONTRACT_CHECK("Assertion", cond)
#endif

#define bx_reach() baryonyx::details::fail_fast("Reached", "unreachable code")

#endif
