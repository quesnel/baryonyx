/* Copyright (C) 2016-2018 INRA
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

#include <baryonyx/core-compare>
#include <baryonyx/core>

#include <fstream>
#include <iterator>
#include <thread>

#include "itm-common.hpp"
#include "itm-solver-common.hpp"
#include "private.hpp"
#include "utils.hpp"

//
// Get number of thread to use in optimizer from parameters list or
// from the standard thread API. If an error occurred, this function
// returns 1.
//
static inline int
get_thread_number(const baryonyx::context_ptr& ctx) noexcept
{
    int thread = ctx->parameters.thread;

    if (thread <= 0)
        thread = static_cast<int>(std::thread::hardware_concurrency());

    if (thread <= 0)
        return 1;

    return thread;
}

namespace baryonyx {

baryonyx::result
solver_select(const baryonyx::context_ptr& ctx, const baryonyx::problem& pb)
{
    return baryonyx::itm::solve(ctx, pb);
}

baryonyx::result
optimizer_select(const baryonyx::context_ptr& ctx, const baryonyx::problem& pb)
{
    auto th = get_thread_number(ctx);

    if (ctx->parameters.auto_tune ==
        solver_parameters::auto_tune_parameters::disabled)
        return baryonyx::itm::optimize(ctx, pb, th);
    else
        return baryonyx::itm::automatic_optimizer(ctx, pb, th);
}
} // namespace baryonyx
