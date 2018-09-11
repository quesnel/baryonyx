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

struct branch_result
{
    branch_result(const baryonyx::problem& pb)
      : problem(pb)
    {}

    branch_result(const baryonyx::problem& pb, const baryonyx::result& r)
      : problem(pb)
      , result(r)
    {}

    baryonyx::problem problem;
    baryonyx::result result;
};

struct branch_result_list
{
    void push_back(const baryonyx::context_ptr& ctx,
                   const baryonyx::problem& pb,
                   int var)
    {
        if (pb.vars.values.size() <= 1)
            return;

        baryonyx::notice(ctx, "- branch-optimization splits\n");

        auto ret = baryonyx::split(ctx, pb, var);

        data.emplace_back(std::move(std::get<0>(ret)));
        data.emplace_back(std::move(std::get<1>(ret)));
    }

    void push_back(const baryonyx::context_ptr& ctx)
    {
        if (data.front().problem.vars.values.size() <= 1)
            return;

        baryonyx::notice(ctx, "- branch-optimization splits\n");

        auto ret = baryonyx::split(
          ctx, data.front().problem, data.front().result.annoying_variable);

        std::ofstream ofs1("brp-001.lp");
        ofs1 << std::get<0>(ret);

        std::ofstream ofs2("brp-002.lp");
        ofs2 << std::get<1>(ret);

        data.emplace_back(std::move(std::get<0>(ret)));
        data.emplace_back(std::move(std::get<1>(ret)));
    }

    void pop_front()
    {
        data.pop_front();
    }

    void optimize(const baryonyx::context_ptr& ctx)
    {
        data.front().result =
          baryonyx::itm::optimize(ctx, data.front().problem);

        if (data.front().result.annoying_variable < 0 ||
            data.front().result.annoying_variable >=
              baryonyx::length(data.front().problem.vars.names))
            data.front().result.annoying_variable =
              baryonyx::length(data.front().problem.vars.names) / 2;
    }

    const branch_result& front() const
    {
        return data.front();
    }

    bool empty() const
    {
        return data.empty();
    }

    std::list<branch_result> data;
};

static bool
is_best(const branch_result& lhs, const branch_result& rhs)
{
    bx_expects(lhs.result && !lhs.result.solutions.empty());
    bx_expects(rhs.result && !rhs.result.solutions.empty());

    if (lhs.problem.type == baryonyx::objective_function_type::maximize) {
        return lhs.result.solutions.front().value <
               rhs.result.solutions.front().value;
    } else {
        return rhs.result.solutions.front().value <
               lhs.result.solutions.front().value;
    }
}

static baryonyx::result
optimize(const baryonyx::context_ptr& ctx, const baryonyx::problem& pb)
{
    // NOTE 1 - add a solver parameter to limit the size of the list.

    auto old_log_priority = ctx->log_priority;
    ctx->log_priority = baryonyx::context::message_type::notice;

    auto ret = baryonyx::itm::optimize(ctx, pb);

    if (ret)
        baryonyx::notice(ctx,
                         "  - branch optimization found solution {:f}\n",
                         ret.solutions.front().value);

    bx_ensures(ret.annoying_variable >= 0 &&
               ret.annoying_variable < baryonyx::length(pb.vars.values));

    branch_result best{ pb, ret };
    branch_result_list jobs;

    jobs.push_back(ctx, pb, ret.annoying_variable);

    while (!jobs.empty()) {
        jobs.optimize(ctx);

        if (!best.result) {
            if (jobs.front().result)
                best = jobs.front();

            jobs.push_back(ctx);
        } else {
            if (jobs.front().result) {
                if (is_best(jobs.front(), best)) {
                    baryonyx::notice(
                      ctx,
                      "  - branch optimization found solution {:f}\n",
                      jobs.front().result.solutions.front().value);

                    best = jobs.front();
                }

                jobs.push_back(ctx);
            }
        }

        jobs.pop_front();
    }

    ctx->log_priority = old_log_priority;

    return best.result;
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
