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

#include <fstream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>

#include <fmt/printf.h>

#include <baryonyx/core-compare>
#include <baryonyx/core-out>
#include <baryonyx/core>

template<typename T>
inline bool
is_essentially_equal(const T v1, const T v2, const T epsilon)
{
    static_assert(std::is_floating_point<T>::value,
                  "is_essentially_equal required a float/double "
                  "as template arguement");

    return fabs((v1) - (v2)) <=
           ((fabs(v1) > fabs(v2) ? fabs(v2) : fabs(v1)) * (epsilon));
}

static void
test_preprocessor()
{
    auto ctx = baryonyx::make_context();

    std::stringstream ss;

    {
        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/prepro.lp");
        context_set_parameter(ctx, "norm", "infinity");
        auto result = baryonyx::solve(ctx, pb);

        Ensures(result.affected_vars.names.size() == 21);
        Ensures(result.affected_vars.values[0] == 0);
        Ensures(result.affected_vars.values[1] == 0);
        Ensures(result.affected_vars.values[2] == 1);

        Ensures(result.variable_value.size() == 2);
        Ensures(result.variable_name.size() == 2);

        Ensures(result.status == baryonyx::result_status::success);
        Ensures(baryonyx::is_valid_solution(pb, result.variable_value) ==
                true);

        Ensures(is_essentially_equal(result.value, 1000.45, 0.01));

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

static void
test_preprocessor_2()
{
    auto ctx = baryonyx::make_context();

    std::stringstream ss;
    double r;

    {
        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/capmo1_direct.lp");
        context_set_parameter(ctx, "preprocessing", "equal,less,greater");
        auto result = baryonyx::solve(ctx, pb);

        fmt::print("result.value {:F} == 6212977\n", result.value);

        r = result.value;
        Ensures(result.value > 6000000);
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

static void
test_real_cost()
{
    auto ctx = baryonyx::make_context();

    const char* str_pb = "minimize\n"
                         "- 0.1 a - 0.5 b - 0.9 c - 1e-7 d\n"
                         "Subject to:\n"
                         "-a -b -c <= -1\n"
                         "-a -b -c >= -3\n"
                         "-a -c >= -2\n"
                         "-a -c <= -1\n"
                         "a + c >= 1\n"
                         "+ b + c +d >= 2\n"
                         "Binaries\n"
                         "a b c d\n"
                         "End\n";

    std::istringstream iss(str_pb);

    auto pb = baryonyx::make_problem(ctx, iss);
    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);

    Ensures(is_essentially_equal(result.value, -0.6, 0.01));

    if (result)
        Ensures(baryonyx::is_valid_solution(pb, result.variable_value) ==
                true);
}

static void
test_assignment_problem()
{
    auto ctx = baryonyx::make_context();

    auto pb =
      baryonyx::make_problem(ctx, EXAMPLES_DIR "/assignment_problem_1.lp");
    context_set_parameter(ctx, "limit", 50);

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
}

static void
test_assignment_problem_random_coast()
{
    auto ctx = baryonyx::make_context();

    context_set_parameter(ctx, "limit", 1000000);
    context_set_parameter(ctx, "theta", 0.5);
    context_set_parameter(ctx, "delta", 0.2);
    context_set_parameter(ctx, "kappa-step", 10e-4);
    context_set_parameter(ctx, "kappa-max", 10.0);
    context_set_parameter(ctx, "alpha", 0.0);
    context_set_parameter(ctx, "w", 20);

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

static void
test_negative_coeff()
{
    auto ctx = baryonyx::make_context();

    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff.lp");

    context_set_parameter(ctx, "limit", 50);

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

static void
test_negative_coeff2()
{
    auto ctx = baryonyx::make_context();

    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff2.lp");

    context_set_parameter(ctx, "limit", 2);

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(result.affected_vars.names.size() == 4);
}

static void
test_negative_coeff3()
{
    auto ctx = baryonyx::make_context();

    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff3.lp");

    context_set_parameter(ctx, "limit", 50);

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

static void
test_negative_coeff4()
{
    auto ctx = baryonyx::make_context();

    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff4.lp");

    context_set_parameter(ctx, "limit", 50);

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

static void
test_negative_coeff5()
{
    auto ctx = baryonyx::make_context();

    const char* str_pb = "minimize\n"
                         "a b c d\n"
                         "Subject to:\n"
                         "-a -b -c <= -1\n"
                         "-a -b -c >= -3\n"
                         "-a -c >= -2\n"
                         "-a -c <= -1\n"
                         "a + c >= 1\n"
                         "+ b + c +d >= 2\n"
                         "Binaries\n"
                         "a b c d\n"
                         "End\n";

    std::istringstream iss(str_pb);

    auto pb = baryonyx::make_problem(ctx, iss);
    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);

    if (result)
        Ensures(baryonyx::is_valid_solution(pb, result.variable_value) ==
                true);
}

static void
test_8_queens_puzzle_fixed_cost()
{
    auto ctx = baryonyx::make_context();

    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/8_queens_puzzle.lp");

    context_set_parameter(ctx, "limit", -1);
    context_set_parameter(ctx, "theta", 0.5);
    context_set_parameter(ctx, "delta", 0.02);
    context_set_parameter(ctx, "kappa-step", 0.01);
    context_set_parameter(ctx, "kappa-max", 60.0);
    context_set_parameter(ctx, "alpha", 1.0);
    context_set_parameter(ctx, "w", 40);

    {
        std::vector<int> cost{ 25, 89, 12, 22, 84, 3,  61, 14, 93, 97, 68,
                               5,  51, 72, 96, 80, 13, 38, 81, 48, 70, 50,
                               66, 68, 30, 97, 79, 4,  41, 44, 47, 62, 60,
                               11, 18, 44, 57, 24, 7,  11, 66, 87, 9,  17,
                               27, 60, 95, 45, 94, 47, 60, 87, 79, 53, 81,
                               52, 91, 53, 57, 8,  63, 78, 1,  8 };

        int i{ 0 };
        for (auto& elem : pb.objective.elements)
            elem.factor = cost[i++];
    }

    auto result = baryonyx::solve(ctx, pb);

    for (int i = 0; i != 8; ++i) {
        for (int j = 0; j != 8; ++j)
            fmt::print("{} ", result.variable_value[j * 8 + i]);

        fmt::print("\n");
    }

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

static void
test_8_queens_puzzle_random_cost()
{
    auto ctx = baryonyx::make_context();

    std::map<std::string, baryonyx::parameter> params;
    context_set_parameter(ctx, "limit", -1);
    context_set_parameter(ctx, "theta", 0.5);
    context_set_parameter(ctx, "delta", 0.02);
    context_set_parameter(ctx, "kappa-step", 0.01);
    context_set_parameter(ctx, "kappa-max", 60.0);
    context_set_parameter(ctx, "alpha", 1.0);
    context_set_parameter(ctx, "w", 40);
    context_set_parameter(
      ctx, "constraint-order", std::string("infeasibility-decr"));
    context_set_parameter(
      ctx, "preprocessing", std::string("variables-weight"));

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

static void
test_qap()
{
    auto ctx = baryonyx::make_context();

    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/small4.lp");

    std::map<std::string, baryonyx::parameter> params;
    context_set_parameter(ctx, "limit", -1);
    context_set_parameter(ctx, "theta", 0.5);
    context_set_parameter(ctx, "delta", 0.01);
    context_set_parameter(ctx, "kappa-step", 10e-4);
    context_set_parameter(ctx, "kappa-max", 10.0);
    context_set_parameter(ctx, "alpha", 0.0);
    context_set_parameter(ctx, "w", 20);

    auto result = baryonyx::solve(ctx, pb);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

static void
test_flat30_7()
{
    auto ctx = baryonyx::make_context();

    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/flat30-7.lp");

    context_set_parameter(ctx, "limit", -1);
    // context_set_parameter(ctx, "delta", 0.001);
    context_set_parameter(ctx, "kappa-min", 0.3);
    context_set_parameter(ctx, "kappa-step", 1e-10);
    context_set_parameter(ctx, "kappa-max", 1.0);
    // context_set_parameter(ctx, "w", 60);

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

static void
test_uf50_0448()
{
    auto ctx = baryonyx::make_context();

    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/uf50-0448.lp");

    std::map<std::string, baryonyx::parameter> params;
    context_set_parameter(ctx, "limit", -1);
    context_set_parameter(ctx, "theta", 0.5);
    context_set_parameter(ctx, "delta", 1.0);
    context_set_parameter(ctx, "kappa-min", 0.1);
    context_set_parameter(ctx, "kappa-step", 1e-17);
    context_set_parameter(ctx, "kappa-max", 1.0);
    context_set_parameter(ctx, "alpha", 2.0);
    context_set_parameter(ctx, "w", 60);
    context_set_parameter(
      ctx, "constraint-order", std::string("random-sorting"));

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

static void
test_aim_50_1_6_yes1_2()
{
    auto ctx = baryonyx::make_context();

    auto pb =
      baryonyx::make_problem(ctx, EXAMPLES_DIR "/aim-50-1_6-yes1-2.lp");

    context_set_parameter(ctx, "limit", -1);
    context_set_parameter(ctx, "theta", 0.6);
    context_set_parameter(ctx, "delta", 0.01);
    context_set_parameter(ctx, "kappa-step", 2 * 10e-4);
    context_set_parameter(ctx, "kappa-max", 100.0);
    context_set_parameter(ctx, "alpha", 1.0);
    context_set_parameter(ctx, "w", 20);

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}

#if 0
static void
test_Z_coefficient_1()
{
    auto ctx = baryonyx::make_context();

    {
        const char* str_pb = "minimize\n"
                             "15 x1 + 19 x2 + 13 x3 + 12 x4\n"
                             "Subject to:\n"
                             "2 x1 + 1 x2 + 3 x3 + 2 x4 <= 3\n"
                             "Binaries\n"
                             "x1 x2 x3 x4\n"
                             "End\n";

        std::istringstream iss(str_pb);

        auto pb = baryonyx::make_problem(ctx, iss);
        auto result = baryonyx::solve(ctx, pb);

        Ensures(result.status == baryonyx::result_status::success);

        if (result)
            Ensures(baryonyx::is_valid_solution(pb, result.variable_value) ==
                    true);
    }

    {
        const char* str_pb = "minimize\n"
                             "Subject to:\n"
                             "2a + 3b -5c + 7d <= 0\n"
                             "-2b + 2c >= 1\n"
                             "Binaries\n"
                             "a b c d\n"
                             "End\n";

        std::istringstream iss(str_pb);

        auto pb = baryonyx::make_problem(ctx, iss);
        auto result = baryonyx::solve(ctx, pb);

        Ensures(result.status == baryonyx::result_status::success);

        if (result)
            Ensures(baryonyx::is_valid_solution(pb, result.variable_value) ==
                    true);
    }
}

static void
test_bibd1n()
{
    auto ctx = baryonyx::make_context();


    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/bibd1n.lp");

    context_set_parameter(ctx, "limit", -1);
    context_set_parameter(ctx, "theta", 0.5);
    context_set_parameter(ctx, "delta", 1e-7);
    context_set_parameter(ctx, "kappa-min", 0.05);
    context_set_parameter(ctx, "kappa-step", 1e-17);
    context_set_parameter(ctx, "kappa-max", 1.0);
    context_set_parameter(ctx, "alpha", 1.0);
    context_set_parameter(ctx, "w", 60);
    context_set_parameter(ctx, "print-level", 1);
    context_set_parameter(ctx, "constraint-order", std::string("random-sorting"));
    context_set_parameter(ctx, "time-limit", -1.0);
    context_set_parameter(ctx, "preprocessing", std::string("none"));
    context_set_parameter(ctx, "floating-point-type", std::string("float"));

    auto result = baryonyx::solve(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    Ensures(baryonyx::is_valid_solution(pb, result.variable_value) == true);
}
#endif

int
main(int /*argc*/, char* /* argv */ [])
{
    test_preprocessor();
    test_preprocessor_2();
    test_real_cost();
    test_assignment_problem();
    test_assignment_problem_random_coast();
    test_negative_coeff();
    test_negative_coeff2();
    test_negative_coeff3();
    test_negative_coeff4();
    test_negative_coeff5();
    test_8_queens_puzzle_fixed_cost();
    test_8_queens_puzzle_random_cost();
    test_qap();
    test_uf50_0448();
    test_flat30_7();
    test_aim_50_1_6_yes1_2();

#if 0
    test_Z_coefficient_1();
    test_bibd1n();
#endif

    return unit_test::report_errors();
}
