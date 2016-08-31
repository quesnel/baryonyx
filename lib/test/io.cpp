/* Copyright (C) 2016 INRA
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

#include <lpcore>
#include <lpcore-compare>
#include <iostream>
#include <fstream>
#include <sstream>

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

void test_examples_1()
{
    const char *example_1 = "maximize\n"
        "x1 + 2x2 + 3x3\n"
        "st\n"
        "time:  -x1 + x2 + x3 <= 20\n"
        "labor:  x1 - 3x2 + x3 <= 30\n"
        "bounds\n"
        "x1 <= 40\n"
        "end\n";

    std::istringstream iss(example_1);

    auto pb = lp::make_problem(iss);

    assert(pb.type == lp::objective_function_type::maximize);
    assert(pb.vars.names.size() == 3);
    assert(pb.vars.values.size() == 3);

    assert(pb.vars.names[0] == "x1");
    assert(pb.vars.names[1] == "x2");
    assert(pb.vars.names[2] == "x3");

    assert(pb.vars.values[0].min == 0);
    assert(pb.vars.values[1].min == 0);
    assert(pb.vars.values[2].min == 0);

    assert(pb.vars.values[0].max == 40);
    assert(pb.vars.values[1].max == std::numeric_limits<int>::max());
    assert(pb.vars.values[2].max == std::numeric_limits<int>::max());
}

void test_examples_2()
{
    std::ifstream ifs();

    auto pb = lp::make_problem(EXAMPLES_DIR "/assignment_problem_4.lp");

    assert(pb.vars.names.size() == 16);
    assert(pb.vars.values.size() == 16);

    std::stringstream ss;
    ss << pb;

    std::printf("%s\n", ss.str().c_str());

    auto pb2 = lp::make_problem(ss);

    assert(pb == pb2);
}

void test_examples_3()
{
    auto pb = lp::make_problem(EXAMPLES_DIR "/geom-30a-3-ext_1000_support.lp");

    assert(pb.type == lp::objective_function_type::minimize);
    assert(pb.vars.names.size() == 819);
    assert(pb.vars.values.size() == 819);

    lp::index nb {0};
    for (auto& elem : pb.vars.values)
        if (elem.type == lp::variable_type::binary)
            ++nb;

    assert(nb == 90);
}

void test_examples_4()
{
    // auto pb = lp::make_problem(EXAMPLES_DIR "/verger.lp");

    // assert(pb.vars.names.size() == 16);
    // assert(pb.vars.values.size() == 16);
}

int main(int /* argc */, char */* argv */[])
{
    test_examples_1();
    test_examples_2();
    test_examples_3();
    test_examples_4();

    return 0;
}
