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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_AUTOMATIC_OPTIMIZE_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_AUTOMATIC_OPTIMIZE_HPP

#include <baryonyx/core>

namespace baryonyx {

/**
 * @brief Auto-tune baryonyx solver parameters to optimize the problem.
 *
 * @details The @c nlopt library or a manual test is used to found the best
 *     parameters @c baryonyx::optimize function for the specified problem.
 *
 * @param ctx A context.
 * @param pb Problem definition.
 * @param thread Number of thread available for running optimization.
 * @return A representation of the result.
 * @throw @c baryonyx::solver_error
 */
result
automatic_optimize(const baryonyx::context_ptr& ctx,
                   baryonyx::problem& pb,
                   int thread);

} // namespace baryonyx

#endif
