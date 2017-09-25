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

#include "unit-test.hpp"

#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>

#include <baryonyx/core-compare>
#include <baryonyx/core-out>
#include <baryonyx/core>

void
test_examples_1(std::shared_ptr<baryonyx::context> ctx)
{
    const char* example_1 = "maximize\n"
                            "obj: x1 + 2x2 + 3x3 - 100\n"
                            "st\n"
                            "time:  -x1 + x2 + x3 <= 20\n"
                            "labor:  x1 - 3x2 + x3 <= 30\n"
                            "test: x1 - 3x2 + x3 <= -5\n"
                            "bounds\n"
                            "x1 <= 40\n"
                            "end\n";

    std::istringstream iss(example_1);

    auto pb = baryonyx::make_problem(ctx, iss);
    std::cout << __func__ << '\n' << baryonyx::resume(pb) << '\n';

    Ensures(pb.type == baryonyx::objective_function_type::maximize);
    Ensures(pb.objective.elements.size() == 3);
    Ensures(pb.objective.elements[0].factor == 1);
    Ensures(pb.objective.elements[0].variable_index == 0);
    Ensures(pb.objective.elements[1].factor == 2);
    Ensures(pb.objective.elements[1].variable_index == 1);
    Ensures(pb.objective.elements[2].factor == 3);
    Ensures(pb.objective.elements[2].variable_index == 2);
    Ensures(pb.objective.constant == -100);

    Ensures(pb.vars.names.size() == 3);
    Ensures(pb.vars.values.size() == 3);
    Ensures(pb.less_constraints.size() == 3);

    Ensures(pb.less_constraints[0].elements.size() == 3);
    Ensures(pb.less_constraints[0].elements[0].factor == -1);
    Ensures(pb.less_constraints[0].elements[0].variable_index == 0);
    Ensures(pb.less_constraints[0].elements[1].factor == 1);
    Ensures(pb.less_constraints[0].elements[1].variable_index == 1);
    Ensures(pb.less_constraints[0].elements[2].factor == 1);
    Ensures(pb.less_constraints[0].elements[2].variable_index == 2);
    Ensures(pb.less_constraints[0].value == 20);

    Ensures(pb.less_constraints[1].elements.size() == 3);
    Ensures(pb.less_constraints[1].elements[0].factor == 1);
    Ensures(pb.less_constraints[1].elements[0].variable_index == 0);
    Ensures(pb.less_constraints[1].elements[1].factor == -3);
    Ensures(pb.less_constraints[1].elements[1].variable_index == 1);
    Ensures(pb.less_constraints[1].elements[2].factor == 1);
    Ensures(pb.less_constraints[1].elements[2].variable_index == 2);
    Ensures(pb.less_constraints[1].value == 30);

    Ensures(pb.less_constraints[2].elements.size() == 3);
    Ensures(pb.less_constraints[1].elements[0].factor == 1);
    Ensures(pb.less_constraints[2].elements[0].variable_index == 0);
    Ensures(pb.less_constraints[2].elements[1].factor == -3);
    Ensures(pb.less_constraints[2].elements[1].variable_index == 1);
    Ensures(pb.less_constraints[2].elements[2].factor == 1);
    Ensures(pb.less_constraints[2].elements[2].variable_index == 2);
    Ensures(pb.less_constraints[2].value == -5);

    Ensures(pb.vars.names[0] == "x1");
    Ensures(pb.vars.names[1] == "x2");
    Ensures(pb.vars.names[2] == "x3");

    Ensures(pb.vars.values[0].min == 0);
    Ensures(pb.vars.values[1].min == 0);
    Ensures(pb.vars.values[2].min == 0);

    Ensures(pb.vars.values[0].max == 40);
    Ensures(pb.vars.values[1].max == std::numeric_limits<int>::max());
    Ensures(pb.vars.values[2].max == std::numeric_limits<int>::max());
}

void
test_examples_2(std::shared_ptr<baryonyx::context> ctx)
{
    std::ifstream ifs;

    std::vector<std::vector<int>> values(3);

    values[0] = { 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };
    values[1] = { 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0 };
    values[2] = { 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0 };

    for (int i = 1; i != 4; ++i) {
        std::string filepath{ EXAMPLES_DIR "/assignment_problem_" };
        filepath += std::to_string(i);
        filepath += ".lp";

        auto pb = baryonyx::make_problem(ctx, filepath);
        std::cout << __func__ << '\n' << baryonyx::resume(pb) << '\n';

        Ensures(pb.vars.names.size() == 16);
        Ensures(pb.vars.values.size() == 16);

        std::stringstream ss;
        ss << pb;

        auto pb2 = baryonyx::make_problem(ctx, ss);
        Ensures(pb == pb2);
    }
}

void
test_examples_3(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb = baryonyx::make_problem(
      ctx, EXAMPLES_DIR "/geom-30a-3-ext_1000_support.lp");
    std::cout << __func__ << '\n' << baryonyx::resume(pb) << '\n';

    Ensures(pb.type == baryonyx::objective_function_type::minimize);
    Ensures(pb.vars.names.size() == 819);
    Ensures(pb.vars.values.size() == 819);

    baryonyx::index nb{ 0 };
    for (auto& elem : pb.vars.values)
        if (elem.type == baryonyx::variable_type::binary)
            ++nb;

    Ensures(nb == 90);
}

void
test_examples_4(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/general.lp");
    std::cout << __func__ << '\n' << baryonyx::resume(pb) << '\n';

    Ensures(pb.type == baryonyx::objective_function_type::minimize);
    Ensures(pb.vars.names.size() == 3);
    Ensures(pb.vars.values.size() == 3);

    baryonyx::index nb{ 0 };
    for (auto& elem : pb.vars.values)
        if (elem.type == baryonyx::variable_type::general)
            ++nb;

    Ensures(nb == 3);
}

void
test_examples_sudoku(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/sudoku.lp");

    Ensures(pb.vars.names.size() == 81);
    Ensures(pb.vars.values.size() == 81);

    for (auto& vv : pb.vars.values) {
        Ensures(vv.min == 1);
        Ensures(vv.max == 9);
        Ensures(vv.type == baryonyx::variable_type::general);
    }

    std::cout << __func__ << '\n' << baryonyx::resume(pb) << '\n';
}

void
test_examples_8_queens_puzzle(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/8_queens_puzzle.lp");

    Ensures(pb.vars.names.size() == 64);
    Ensures(pb.vars.values.size() == 64);

    for (auto& vv : pb.vars.values) {
        Ensures(vv.min == 0);
        Ensures(vv.max == 1);
        Ensures(vv.type == baryonyx::variable_type::binary);
    }

    std::cout << __func__ << '\n' << baryonyx::resume(pb) << '\n';
}

void
test_examples_vm(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/vm.lp");

    std::cout << __func__ << '\n' << baryonyx::resume(pb) << '\n';
}

void
test_verger_5_5(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/verger_5_5.lp");

    std::cout << __func__ << '\n' << baryonyx::resume(pb) << '\n';
}

int
main(int /* argc */, char* /* argv */ [])
{
    auto ctx = std::make_shared<baryonyx::context>();
    ctx->set_standard_stream_logger();

    test_examples_1(ctx);
    test_examples_2(ctx);
    test_examples_3(ctx);
    test_examples_4(ctx);
    test_examples_sudoku(ctx);
    test_examples_8_queens_puzzle(ctx);
    test_examples_vm(ctx);
    test_verger_5_5(ctx);

    return unit_test::report_errors();
}
