/* Copyright (C) 2016-2019 INRA
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

#include "problem.hpp"
#include "unit-test.hpp"
#include "utils.hpp"

#include <baryonyx/core>

int
get_variable(const std::vector<std::string_view>& variable_names,
             std::string_view name)
{
    for (int i = 0, e = baryonyx::length(variable_names); i != e; ++i)
        if (variable_names[i] == name)
            return i;

    return -1;
}

void
test_bound_affectation()
{
    auto ctx = baryonyx::make_context(stdout, 7);

    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/bound.lp");

    Ensures(pb);
    Ensures(pb.vars.names.size() == static_cast<size_t>(6));

    // The preprocessing can not remove any constraint or variable from the
    // raw_problem.

    auto pb_pp = baryonyx::preprocess(ctx, pb);
    Ensures(pb_pp.vars.names.size() == static_cast<size_t>(6));
    Ensures(pb_pp.affected_vars.names.size() == static_cast<size_t>(0));
    Ensures(pb_pp.equal_constraints.size() == static_cast<size_t>(0));
    Ensures(pb_pp.less_constraints.size() == static_cast<size_t>(2));
    Ensures(pb_pp.greater_constraints.size() == static_cast<size_t>(0));

    // We add an "affectation bound" for the variable 'f' i.e. in bound section
    // of the lp file format:
    //
    // bounds:
    // f = 0

    auto id = get_variable(pb.vars.names, "f");
    Ensures(id >= 0);
    Ensures(id < baryonyx::length(pb.vars.names));

    Ensures(pb.vars.names.size() == pb.vars.values.size());

    pb.vars.values[id].min = 0;
    pb.vars.values[id].max = 0;

    auto pb_pp2 = baryonyx::preprocess(ctx, pb);
    Ensures(pb_pp2.vars.names.size() == static_cast<size_t>(5));
    Ensures(pb_pp2.affected_vars.names.size() == static_cast<size_t>(1));
    Ensures(pb_pp2.equal_constraints.size() == static_cast<size_t>(0));
    Ensures(pb_pp2.less_constraints.size() == static_cast<size_t>(2));
    Ensures(pb_pp2.greater_constraints.size() == static_cast<size_t>(0));
}

void
test_cleaning_affected_variables()
{
    auto ctx = baryonyx::make_context();
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/prepro.lp");

    Ensures(pb);
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

    Ensures(pb);
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
        // d = 0 implies e = {0, 1} but minimize so => 0.

        auto pb_af = baryonyx::affect(ctx, pb_pp, 0, false);
        Ensures(pb_af.vars.names.size() == static_cast<size_t>(0));
        Ensures(pb_af.affected_vars.names.size() == static_cast<size_t>(23));
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
        // e = 0 implies d = {0, 1} but minimize so => 0.

        auto pb_af = baryonyx::affect(ctx, pb_pp, 1, false);
        Ensures(pb_af.vars.names.size() == static_cast<size_t>(0));
        Ensures(pb_af.affected_vars.names.size() == static_cast<size_t>(23));
    }
}

void
test_split()
{
    auto ctx = baryonyx::make_context();
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/prepro.lp");
    Ensures(pb);
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

        Ensures(p1.vars.names.size() == static_cast<size_t>(0));
        Ensures(p1.affected_vars.names.size() == static_cast<size_t>(23));
    }
}

int
main(int /*argc*/, char* /* argv */ [])
{
    unit_test::checks("bound_affectation", test_bound_affectation);
    unit_test::checks("cleaning_affected_variables",
                      test_cleaning_affected_variables);
    unit_test::checks("affect_variable", test_affect_variable);
    unit_test::checks("split", test_split);

    return unit_test::report_errors();
}
