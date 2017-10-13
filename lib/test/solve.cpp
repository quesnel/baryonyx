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
#include <random>
#include <sstream>

#include <baryonyx/core-compare>
#include <baryonyx/core-out>
#include <baryonyx/core>

void
test_preprocessor(std::shared_ptr<baryonyx::context> ctx)
{
    std::stringstream ss;

    {
        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/prepro.lp");
        ctx->set_parameter("norm", "infinity");
        auto result = baryonyx::solve(ctx, pb);

        Ensures(result.affected_vars.names.size() == 20);
        Ensures(result.affected_vars.values[0] == 0);
        Ensures(result.affected_vars.values[1] == 0);
        Ensures(result.affected_vars.values[2] == 1);

        Ensures(result.variable_value.size() == 2);
        Ensures(result.variable_name.size() == 2);

        Ensures(result.status == baryonyx::result_status::success);
        Ensures(baryonyx::is_valid_solution(pb, result.variable_value) ==
                true);

        ss << result;
        if (not ss.good())
            Ensures(ss.good());
    }

    {
        ss.seekg(0, std::ios::beg);

        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/prepro.lp");
        auto re = baryonyx::make_result(ctx, ss);

        Ensures(is_valid_solution(pb, re));
    }
}

void
test_preprocessor_2(std::shared_ptr<baryonyx::context> ctx)
{
    std::stringstream ss;
    double r;

    {
        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/capmo1_direct.lp");
        ctx->set_parameter("preprocessing", "equal,less,greater");
        auto result = baryonyx::solve(ctx, pb);

        r = result.value;
        Ensures(result.value == 6212977);
        Ensures(baryonyx::is_valid_solution(pb, result.variable_value) ==
                true);

        ss << result;
        if (not ss.good())
            Ensures(ss.good());
    }

    {
        ss.seekg(0, std::ios::beg);

        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/capmo1_direct.lp");

        baryonyx::result result = baryonyx::make_result(ctx, ss);

        result.status = baryonyx::result_status::success;

        Ensures(baryonyx::is_valid_solution(pb, result));
        Ensures(baryonyx::compute_solution(pb, result) == r);
    }
}

void
test_assignment_problem(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb =
      baryonyx::make_problem(ctx, EXAMPLES_DIR "/assignment_problem_1.lp");
    ctx->set_parameter("limit", 50);

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

void
test_assignment_problem_random_coast(std::shared_ptr<baryonyx::context> ctx)
{
    ctx->set_parameter("limit", 1000000);
    ctx->set_parameter("theta", 0.5);
    ctx->set_parameter("delta", 0.2);
    ctx->set_parameter("kappa-step", 10e-4);
    ctx->set_parameter("kappa-max", 10.0);
    ctx->set_parameter("alpha", 0.0);
    ctx->set_parameter("w", 20);

    for (int i{ 0 }, e{ 10 }; i != e; ++i) {
        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/assignment_problem_1.lp");

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(1, 100);

        for (auto& elem : pb.objective.elements)
            elem.factor = dis(gen);

        auto result = baryonyx::solve(ctx, pb);

        Ensures(result.status == baryonyx::result_status::success);
        Ensures(baryonyx::is_valid_solution(pb, result.variable_value) ==
                true);
    }
}

void
test_negative_coeff(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff.lp");

    ctx->set_parameter("limit", 50);

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

void
test_negative_coeff2(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff2.lp");

    ctx->set_parameter("limit", 2);

    auto result = baryonyx::solve(ctx, pb);

    for (auto elen : result.variable_value)
        std::cout << elen << ' ';
    std::cout << '\n';

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(result.variable_value[0] == 1);
    Ensures(result.variable_value[1] == 0);
    Ensures(result.variable_value[2] == 0);
    Ensures(result.variable_value[3] == 1);

    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

void
test_negative_coeff3(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff3.lp");

    ctx->set_parameter("limit", 50);

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

void
test_negative_coeff4(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff4.lp");

    ctx->set_parameter("limit", 50);

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

void
test_flat30_7(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/flat30-7.lp");

    ctx->set_parameter("limit", -1);
    ctx->set_parameter("delta", 0.001);
    ctx->set_parameter("kappa-min", 0.3);
    ctx->set_parameter("kappa-step", 1e-4);
    ctx->set_parameter("kappa-max", 10.0);
    ctx->set_parameter("w", 60);

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

void
test_uf50_0448(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/uf50-0448.lp");

    std::map<std::string, baryonyx::parameter> params;
    ctx->set_parameter("limit", -1);
    ctx->set_parameter("theta", 0.5);
    ctx->set_parameter("delta", 1.0);
    ctx->set_parameter("kappa-min", 0.1);
    ctx->set_parameter("kappa-step", 1e-17);
    ctx->set_parameter("kappa-max", 1.0);
    ctx->set_parameter("alpha", 2.0);
    ctx->set_parameter("w", 60);
    ctx->set_parameter("constraint-order", std::string("random-sorting"));

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

void
test_aim_50_1_6_yes1_2(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb =
      baryonyx::make_problem(ctx, EXAMPLES_DIR "/aim-50-1_6-yes1-2.lp");

    ctx->set_parameter("limit", -1);
    ctx->set_parameter("theta", 0.6);
    ctx->set_parameter("delta", 0.01);
    ctx->set_parameter("kappa-step", 2 * 10e-4);
    ctx->set_parameter("kappa-max", 100.0);
    ctx->set_parameter("alpha", 1.0);
    ctx->set_parameter("w", 20);

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

void
test_bibd1n(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/bibd1n.lp");

    ctx->set_parameter("limit", -1);
    ctx->set_parameter("theta", 0.5);
    ctx->set_parameter("delta", 0.000001);
    ctx->set_parameter("kappa-min", 0.0005);
    ctx->set_parameter("kappa-step", 1e-17);
    ctx->set_parameter("kappa-max", 1.0);
    ctx->set_parameter("alpha", 2.0);
    ctx->set_parameter("w", 60);
    ctx->set_parameter("serialize", 1);
    ctx->set_parameter("constraint-order", std::string("reversing"));
    ctx->set_parameter("time-limit", 0.0);
    ctx->set_parameter("preprocessing", "variables-number");

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

void
test_8_queens_puzzle_fixed_cost(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/8_queens_puzzle.lp");

    ctx->set_parameter("limit", -1);
    ctx->set_parameter("theta", 0.5);
    ctx->set_parameter("delta", 0.02);
    ctx->set_parameter("kappa-step", 0.01);
    ctx->set_parameter("kappa-max", 60.0);
    ctx->set_parameter("alpha", 1.0);
    ctx->set_parameter("w", 40);

    std::vector<int> cost{ 25, 89, 12, 22, 84, 3,  61, 14, 93, 97, 68, 5,  51,
                           72, 96, 80, 13, 38, 81, 48, 70, 50, 66, 68, 30, 97,
                           79, 4,  41, 44, 47, 62, 60, 11, 18, 44, 57, 24, 7,
                           11, 66, 87, 9,  17, 27, 60, 95, 45, 94, 47, 60, 87,
                           79, 53, 81, 52, 91, 53, 57, 8,  63, 78, 1,  8 };

    std::size_t i{ 0 };

    for (auto& elem : pb.objective.elements)
        elem.factor = cost[i++];

    auto result = baryonyx::solve(ctx, pb);

    for (int i = 0; i != 8; ++i) {
        for (int j = 0; j != 8; ++j) {
            std::cout << result.variable_value[j * 8 + i] << ' ';
        }
        std::cout << '\n';
    }

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

void
test_8_queens_puzzle_random_cost(std::shared_ptr<baryonyx::context> ctx)
{
    std::map<std::string, baryonyx::parameter> params;
    ctx->set_parameter("limit", -1);
    ctx->set_parameter("theta", 0.5);
    ctx->set_parameter("delta", 0.02);
    ctx->set_parameter("kappa-step", 0.01);
    ctx->set_parameter("kappa-max", 60.0);
    ctx->set_parameter("alpha", 1.0);
    ctx->set_parameter("w", 40);
    ctx->set_parameter("constraint-order", std::string("infeasibility-decr"));
    ctx->set_parameter("preprocessing", std::string("variables-weight"));

    for (int i{ 0 }, e{ 10 }; i != e; ++i) {
        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/8_queens_puzzle.lp");

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(1, 100);

        for (auto& elem : pb.objective.elements)
            elem.factor = dis(gen);

        auto result = baryonyx::solve(ctx, pb);

        Ensures(result.status == baryonyx::result_status::success);
        Ensures(baryonyx::is_valid_solution(pb, result.variable_value) ==
                true);
    }
}

void
test_qap(std::shared_ptr<baryonyx::context> ctx)
{
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/small4.lp");

    std::map<std::string, baryonyx::parameter> params;
    ctx->set_parameter("limit", -1);
    ctx->set_parameter("theta", 0.5);
    ctx->set_parameter("delta", 0.01);
    ctx->set_parameter("kappa-step", 10e-4);
    ctx->set_parameter("kappa-max", 10.0);
    ctx->set_parameter("alpha", 0.0);
    ctx->set_parameter("w", 20);

    auto result = baryonyx::solve(ctx, pb);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

int
main(int /* argc */, char* /* argv */ [])
{
    auto ctx = std::make_shared<baryonyx::context>();
    ctx->set_standard_stream_logger();

    test_preprocessor(ctx);
    test_preprocessor_2(ctx);
    test_assignment_problem(ctx);
    test_assignment_problem_random_coast(ctx);
    test_negative_coeff(ctx);
    test_negative_coeff2(ctx);
    test_negative_coeff3(ctx);
    test_negative_coeff4(ctx);
    test_8_queens_puzzle_fixed_cost(ctx);
    test_8_queens_puzzle_random_cost(ctx);
    test_qap(ctx);
    test_uf50_0448(ctx);
    test_flat30_7(ctx);
    test_aim_50_1_6_yes1_2(ctx);

    // test_bibd1n(ctx);

    return unit_test::report_errors();
}
