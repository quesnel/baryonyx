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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_PRIVATE_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_PRIVATE_HPP

#include <baryonyx/core>

static inline void
#if defined(__GNUC__)
  __attribute__((always_inline, format(printf, 2, 3)))
#endif
  lp_log_null(const std::shared_ptr<baryonyx::context>&, const char*, ...)
{
}

static inline void
#if defined(__GNUC__)
  __attribute__((always_inline, format(printf, 2, 3)))
#endif
  lp_log_null(const baryonyx::context*, const char*, ...)
{
}

#define lp_log_cond(ctx, prio, arg...)                                        \
    do {                                                                      \
        if (ctx->get_log_priority() >= prio) {                                \
            ctx->log(prio, __FILE__, __LINE__, __FUNCTION__, ##arg);          \
        }                                                                     \
    } while (0)

#ifndef BARYONYX_FULL_OPTIMIZATION
#ifndef BARYONYX_DISABLE_DEBUG
#define lp_debug(ctx, arg...) lp_log_cond(ctx, 7, ##arg)
#else
#define lp_debug(ctx, arg...) lp_log_null(ctx, ##arg)
#endif
#else
#define lp_debug(ctx, arg...) lp_log_null(ctx, ##arg)
#endif

#define BARYONYX_UNROLL_PTR(first, last, op)                                  \
    do {                                                                      \
        for (; (last) - (first) >= 4; (first) += 4) {                         \
            op(((first) + 0));                                                \
            op(((first) + 1));                                                \
            op(((first) + 2));                                                \
            op(((first) + 3));                                                \
        }                                                                     \
        for (; (first) != (last); ++(first))                                  \
            op((first));                                                      \
    } while (0);

namespace baryonyx {

struct maximize_tag
{
};

struct minimize_tag
{
};

} // namespace baryonyx

namespace baryonyx_private {

baryonyx::problem
read_problem(std::istream& is);

baryonyx::result
read_result(std::istream& is);

bool
write_problem(std::ostream& os, const baryonyx::problem& pb);

bool
check_consistency(const baryonyx::problem& pb);

void
preprocess(const std::shared_ptr<baryonyx::context>& ctx,
           baryonyx::problem& pb);

baryonyx::result
solve(const std::shared_ptr<baryonyx::context>& ctx, baryonyx::problem& pb);

baryonyx::result
optimize(std::shared_ptr<baryonyx::context> ctx, baryonyx::problem& pb);

} // namespace baryonyx

#endif
