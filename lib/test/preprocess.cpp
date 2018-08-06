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

#include "unit-test.hpp"

#include <fmt/printf.h>

#include <baryonyx/core>

void
test_cleaning_affected_variables()
{
    auto ctx = baryonyx::make_context();
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/prepro.lp");

    Ensures(pb.vars.names.size() == static_cast<size_t>(23));

    auto pb_pp = baryonyx::preprocess(ctx, pb);

    Ensures(pb_pp.vars.names.size() == static_cast<size_t>(2));
    Ensures(pb_pp.affected_vars.names.size() == static_cast<size_t>(21));
}

void
test_affect_variable()
{
    auto ctx = baryonyx::make_context();
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/prepro.lp");
    auto pb_pp = baryonyx::preprocess(ctx, pb);

    Ensures(pb_pp.vars.names.size() == static_cast<size_t>(2));
    Ensures(pb_pp.affected_vars.names.size() == static_cast<size_t>(21));

    {
        // + d + e <= 1
        // d = 1 implies e = 0

        auto pb_af = baryonyx::affect(ctx, pb_pp, 0, true);
        Ensures(pb_af.vars.names.size() == static_cast<size_t>(0));
        Ensures(pb_af.affected_vars.names.size() == static_cast<size_t>(23));
    }

    {
        // + d + e <= 1
        // d = 0 implies e = {0, 1}

        auto pb_af = baryonyx::affect(ctx, pb_pp, 0, false);
        Ensures(pb_af.vars.names.size() == static_cast<size_t>(1));
        Ensures(pb_af.affected_vars.names.size() == static_cast<size_t>(22));
    }

    {
        // + d + e <= 1
        // e = 1 implies d = 0

        auto pb_af = baryonyx::affect(ctx, pb_pp, 1, true);
        Ensures(pb_af.vars.names.size() == static_cast<size_t>(0));
        Ensures(pb_af.affected_vars.names.size() == static_cast<size_t>(23));
    }

    {
        // + d + e <= 1
        // e = 0 implies d = {0, 1}

        auto pb_af = baryonyx::affect(ctx, pb_pp, 1, false);
        Ensures(pb_af.vars.names.size() == static_cast<size_t>(1));
        Ensures(pb_af.affected_vars.names.size() == static_cast<size_t>(22));
    }
}

void
test_split()
{
    auto ctx = baryonyx::make_context();
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/prepro.lp");
    auto pb_pp = baryonyx::preprocess(ctx, pb);

    Ensures(pb_pp.vars.names.size() == static_cast<size_t>(2));
    Ensures(pb_pp.affected_vars.names.size() == static_cast<size_t>(21));

    baryonyx::problem p0, p1;
    std::tie(p0, p1) = baryonyx::split(ctx, pb_pp, 0);

    {
        // + d + e <= 1
        // d = 1 implies e = 0

        Ensures(p0.vars.names.size() == static_cast<size_t>(0));
        Ensures(p0.affected_vars.names.size() == static_cast<size_t>(23));
    }

    {
        // + d + e <= 1
        // d = 0 implies e = {0, 1}

        Ensures(p1.vars.names.size() == static_cast<size_t>(1));
        Ensures(p1.affected_vars.names.size() == static_cast<size_t>(22));
    }
}

int
main(int /*argc*/, char* /* argv */ [])
{
    test_cleaning_affected_variables();
    test_affect_variable();
    test_split();

    return unit_test::report_errors();
}
