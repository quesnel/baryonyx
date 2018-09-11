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

#include "debug.hpp"
#include "itm.hpp"
#include "private.hpp"
#include "problem.hpp"
#include "utils.hpp"

#include <fstream>
#include <list>

/**
 * @brief Return true if lhs is better to rhs.
 */
static bool
is_best(const baryonyx::result& lhs,
        const baryonyx::result& rhs,
        baryonyx::objective_function_type type)
{
    if (lhs) {
        if (rhs) {
            if (type == baryonyx::objective_function_type::maximize) {
                return lhs.solutions.front().value >
                       rhs.solutions.front().value;
            } else {
                return lhs.solutions.front().value <
                       rhs.solutions.front().value;
            }
        } else {
            return true;
        }
    } else {
        if (rhs) {
            return false;
        } else {
            return lhs.constraints < rhs.constraints;
        }
    }
}

template<typename Iterator>
auto
find_best_it(Iterator first,
             Iterator last,
             const baryonyx::result& current,
             baryonyx::objective_function_type type) -> Iterator
{
    for (; first != last; ++first)
        if (is_best(current, first->second, type))
            return first;

    return last;
};

static baryonyx::result
optimize(const baryonyx::context_ptr& ctx, const baryonyx::problem& pb)
{
    // NOTE 1 - add a solver parameter to limit the size of the list.

    auto old_log_priority = ctx->log_priority;
    ctx->log_priority = baryonyx::context::message_type::notice;

    auto best = baryonyx::itm::optimize(ctx, pb);

    if (best)
        baryonyx::notice(ctx,
                         "  - branch optimization found solution {:f}\n",
                         best.solutions.front().value);

    bx_ensures(best.annoying_variable >= 0 &&
               best.annoying_variable < baryonyx::length(pb.vars.values));

    baryonyx::notice(ctx,
                     "- branch-optimization splits on {}\n",
                     pb.vars.names[best.annoying_variable]);

    auto sp = baryonyx::split(ctx, pb, best.annoying_variable);

    std::list<std::pair<baryonyx::problem, baryonyx::result>> jobs;
    jobs.emplace_front(std::move(std::get<0>(sp)), best);
    jobs.emplace_back(std::move(std::get<1>(sp)), best);

    while (!jobs.empty()) {
        auto top = jobs.front().first;
        auto current = baryonyx::itm::optimize(ctx, top);

        baryonyx::notice(ctx,
                         "- branch-optimization splits on {}\n",
                         top.vars.names[current.annoying_variable]);

        auto sp = baryonyx::split(ctx, top, current.annoying_variable);

        jobs.pop_front();

        if ((!best && current) ||
            (best && current && is_best(current, best, top.type))) {
            baryonyx::notice(ctx,
                             "  - branch optimization found solution {:f}\n",
                             current.solutions.front().value);

            best = current;
            jobs.emplace_front(std::move(std::get<0>(sp)), current);
            jobs.emplace_front(std::move(std::get<1>(sp)), current);
        } else {
            auto it =
              find_best_it(jobs.begin(), jobs.end(), current, top.type);

            it = jobs.emplace(it, std::move(std::get<0>(sp)), current);
            jobs.emplace(it, std::move(std::get<1>(sp)), current);
        }
    }

    ctx->log_priority = old_log_priority;

    return best;
}

namespace baryonyx {
namespace itm {

result
branch_optimize(const baryonyx::context_ptr& ctx, const baryonyx::problem& pb)
{
    baryonyx::notice(ctx, "- branch-optimization starts\n");

    return ::optimize(ctx, pb);
}

} // namespace itm
} // namespace baryonyx
