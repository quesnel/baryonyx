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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_DEBUG_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_DEBUG_HPP

#include <exception>

#include <cstdio>

#define bx_stringify_detail(x) #x
#define bx_stringify(x) bx_stringify_detail(x)

#if defined(__clang__) || defined(__GNUC__)
#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
#else
#define likely(x) (!!(x))
#define unlikely(x) (!!(x))
#endif

//
// In internal API, we prefer abort application when precondition,
// postcondition or assertion fail. For external API, user can select to use
// exception or return with message in terminal.
//

namespace baryonyx {
namespace details {

inline void
print(const char* type,
      const char* cond,
      const char* file,
      const char* line) noexcept
{
    fprintf(stderr, "%s [%s] failure at %s: %s\n", type, cond, file, line);
}

[[noreturn]] inline auto
fail_fast(const char* type,
          const char* cond,
          const char* file,
          const char* line) -> void
{
    print(type, cond, file, line);

#ifdef BARYONYX_FAIL_FAST
    std::terminate();
#else
    throw std::logic_error(cond);
#endif
}

} // namespace details
} // namespace baryonyx

#define BX_CONTRACT_FAIL(type)                                                \
    baryonyx::details::fail_fast(                                             \
      type, bx_stringify(cond), __FILE__, bx_stringify(__LINE__))

#define BX_CONTRACT_CHECK(type, cond)                                         \
    (likely(cond)                                                             \
       ? static_cast<void>(0)                                                 \
       : baryonyx::details::fail_fast(                                        \
           type, bx_stringify(cond), __FILE__, bx_stringify(__LINE__)))

#define BX_CONTRACT_CHECK_RETURN_VAL(type, cond, val)                         \
    do {                                                                      \
        if (unlikely(!(cond))) {                                              \
            baryonyx::details::print(                                         \
              type, bx_stringify(cond), __FILE__, bx_stringify(__LINE__));    \
            return val;                                                       \
        }                                                                     \
    } while (0)

#define BX_CONTRACT_CHECK_RETURN(type, cond)                                  \
    do {                                                                      \
        if (unlikely(!(cond))) {                                              \
            baryonyx::details::print(                                         \
              type, bx_stringify(cond), __FILE__, bx_stringify(__LINE__));    \
            return;                                                           \
        }                                                                     \
    } while (0)

#ifdef BARYONYX_FULL_OPTIMIZATION
#define bx_expects(cond)
#define bx_ensures(cond)
#define bx_assert(cond)
#define bx_reach()
#else
#define bx_expects(cond) BX_CONTRACT_CHECK("Precondition", cond)
#define bx_ensures(cond) BX_CONTRACT_CHECK("Postcondition", cond)
#define bx_assert(cond) BX_CONTRACT_CHECK("Assertion", cond)
#define bx_reach() BX_CONTRACT_FAIL("Reached")
#endif

#define bx_return_val_if_fail(cond, val)                                      \
    BX_CONTRACT_CHECK_RETURN_VAL("Precondition", cond, val)
#define bx_return_if_fail(cond) BX_CONTRACT_CHECK_RETURN("Precondition", cond)

#endif
