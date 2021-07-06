/* Copyright (C) 2016-2021 INRAE
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
#include "result.hpp"
#include "resume.hpp"

#include <baryonyx/core-compare>
#include <baryonyx/core-out>
#include <baryonyx/core>

#include <boost/ut.hpp>

#include <fmt/format.h>
#include <fmt/ostream.h>

#include <fstream>
#include <map>
#include <numeric>
#include <random>
#include <sstream>

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

int
get_variable(const baryonyx::result& r, std::string variable_name)
{
    for (int i = 0, e = static_cast<int>(r.affected_vars.values.size());
         i != e;
         ++i)
        if (r.affected_vars.names[i] == variable_name)
            return i;

    return -1;
}

bool
get_value(const baryonyx::result& r, int variable)
{
    return r.affected_vars.values[variable];
}

int
main()
{
    using namespace boost::ut;

    "test_preprocessor"_test = [] {
        auto ctx = baryonyx::make_context(4);

        std::stringstream ss;

        {
            auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/prepro.lp");
            baryonyx::solver_parameters params;
            params.cost_norm =
              baryonyx::solver_parameters::cost_norm_type::loo;
            baryonyx::context_set_solver_parameters(ctx, params);

            fmt::print("{}", pb);

            auto result = baryonyx::solve(ctx, pb);
            fmt::print("{}\naffected vars: {}\n",
                       result,
                       result.affected_vars.names.size());

            expect(result.affected_vars.names.size() == 21);

            int var;

            var = get_variable(result, "w");
            expect(var >= 0);
            expect(get_value(result, var) == false);

            var = get_variable(result, "a");
            expect(var >= 0);
            expect(get_value(result, var) == false);

            var = get_variable(result, "t");
            expect(var >= 0);
            expect(get_value(result, var) == false);

            var = get_variable(result, "ZZ");
            expect(var >= 0);
            expect(get_value(result, var) == true);

            var = get_variable(result, "c1");
            expect(var >= 0);
            expect(get_value(result, var) == false);

            var = get_variable(result, "c2");
            expect(var >= 0);
            expect(get_value(result, var) == false);

            var = get_variable(result, "c3");
            expect(var >= 0);
            expect(get_value(result, var) == false);

            var = get_variable(result, "c4");
            expect(var >= 0);
            expect(get_value(result, var) == false);

            var = get_variable(result, "c5");
            expect(var >= 0);
            expect(get_value(result, var) == false);

            var = get_variable(result, "c6");
            expect(var >= 0);
            expect(get_value(result, var) == false);

            var = get_variable(result, "d1");
            expect(var >= 0);
            expect(get_value(result, var) == true);

            var = get_variable(result, "d2");
            expect(var >= 0);
            expect(get_value(result, var) == true);

            var = get_variable(result, "d3");
            expect(var >= 0);
            expect(get_value(result, var) == true);

            var = get_variable(result, "d4");
            expect(var >= 0);
            expect(get_value(result, var) == false);

            var = get_variable(result, "d5");
            expect(var >= 0);
            expect(get_value(result, var) == false);

            var = get_variable(result, "d6");
            expect(var >= 0);
            expect(get_value(result, var) == false);

            var = get_variable(result, "b");
            expect(var >= 0);
            expect(get_value(result, var) == true);

            expect(result.status == baryonyx::result_status::success);
            expect(!result.solutions.empty());
            expect(result.solutions.size() >= 1);
            expect(result.variable_name.size() == 2);

            expect(baryonyx::is_valid_solution(pb, result) == true);

            expect(result.solutions.back().value > 6.0);

            fmt::print(ss, "{}", result);
            if (!ss.good())
                expect(ss.good());
        }

        {
            ss.seekg(0, std::ios::beg);

            auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/prepro.lp");

            baryonyx::result re;

            ss >> re;

            expect(is_valid_solution(pb, re));
        }
    };

    "test_preprocessor_2"_test = [] {
        auto ctx = baryonyx::make_context(4);

        std::stringstream ss;
        double r;

        {
            auto pb =
              baryonyx::make_problem(ctx, EXAMPLES_DIR "/capmo1_direct.lp");

            auto result = baryonyx::solve(ctx, pb);

            expect(!!result);
            expect(baryonyx::is_valid_solution(pb, result) == true);
            expect(baryonyx::compute_solution(pb, result) ==
                   result.solutions.back().value);

            expect(!result.solutions.empty());
            fmt::print("result.value {:F} == 6212977\n",
                       result.solutions.back().value);

            r = result.solutions.back().value;
            expect(result.solutions.back().value < 1156908);

            fmt::print(ss, "{}", result);
            if (!ss.good())
                expect(ss.good());
        }

        {
            ss.seekg(0, std::ios::beg);

            auto pb =
              baryonyx::make_problem(ctx, EXAMPLES_DIR "/capmo1_direct.lp");

            baryonyx::result result;
            ss >> result;

            result.status = baryonyx::result_status::success;

            expect(baryonyx::is_valid_solution(pb, result));
            expect(baryonyx::compute_solution(pb, result) == r);
        }
    };

    "test_quadratic_preprocessor"_test = [] {
        auto ctx = baryonyx::make_context(4);

        const char* str_pb = "minimize\n"
                             "-5a + [ 2a * b + 3c * d ] /2 + 10b\n"
                             "Subject to:\n"
                             "a + b >= 0\n"
                             "b + c + d >= 0\n"

                             "Binaries\n"
                             "a b c d\n"
                             "End\n";

        std::istringstream iss(str_pb);
        auto pb = baryonyx::make_problem(*ctx, iss);
        expect(pb);

        auto result = baryonyx::solve(ctx, pb);
        expect(!!result);
    };

    "test_real_cost"_test = [] {
        auto ctx = baryonyx::make_context(4);

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

        auto pb = baryonyx::make_problem(*ctx, iss);
        auto result = baryonyx::solve(ctx, pb);

        expect(!!result);
        expect(result.status == baryonyx::result_status::success);
        expect(!result.solutions.empty());

        expect(result.solutions.back().value < 0.0);

        if (result)
            expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    "test_assignment_problem"_test = [] {
        auto ctx = baryonyx::make_context(4);

        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/assignment_problem_1.lp");

        baryonyx::solver_parameters params;
        // params.limit = 50;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);
    };

    "test_assignment_problem_with_buffer"_test = [] {
        auto ctx = baryonyx::make_context(4);

        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/assignment_problem_1.lp");

        baryonyx::solver_parameters params;
        params.use_buffer_solver = true;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);
    };

    "test_assignment_problem_random_coast"_test = [] {
        auto ctx = baryonyx::make_context(4);

        baryonyx::solver_parameters params;
        params.limit = 1000000;
        params.theta = 0.5;
        params.delta = 0.2;
        params.kappa_step = 10e-4;
        params.kappa_max = 10.0;
        params.alpha = 0.0;
        params.w = 20;
        baryonyx::context_set_solver_parameters(ctx, params);

        for (int i{ 0 }, e{ 10 }; i != e; ++i) {
            auto pb = baryonyx::make_problem(
              ctx, EXAMPLES_DIR "/assignment_problem_1.lp");

            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> dis(1, 100);

            for (auto& elem : pb.objective.elements)
                elem.factor = dis(gen);

            auto result = baryonyx::solve(ctx, pb);

            expect(result.status == baryonyx::result_status::success);
            expect(baryonyx::is_valid_solution(pb, result) == true);
        }
    };

    "test_assignment_problem_random_coast_with_buffer"_test = [] {
        auto ctx = baryonyx::make_context(4);

        baryonyx::solver_parameters params;
        params.limit = 1000000;
        params.theta = 0.5;
        params.delta = 0.2;
        params.kappa_step = 10e-4;
        params.kappa_max = 10.0;
        params.alpha = 0.0;
        params.w = 20;
        params.use_buffer_solver = true;
        baryonyx::context_set_solver_parameters(ctx, params);

        for (int i{ 0 }, e{ 10 }; i != e; ++i) {
            auto pb = baryonyx::make_problem(
              ctx, EXAMPLES_DIR "/assignment_problem_1.lp");

            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> dis(1, 100);

            for (auto& elem : pb.objective.elements)
                elem.factor = dis(gen);

            auto result = baryonyx::solve(ctx, pb);

            expect(result.status == baryonyx::result_status::success);
            expect(baryonyx::is_valid_solution(pb, result) == true);
        }
    };

    "test_negative_coeff"_test = [] {
        auto ctx = baryonyx::make_context(4);
        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff.lp");

        baryonyx::solver_parameters params;
        params.use_buffer_solver = false;
        params.limit = 50;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);
        expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    "test_negative_coeff_buffer"_test = [] {
        auto ctx = baryonyx::make_context(4);
        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff.lp");

        baryonyx::solver_parameters params;
        params.use_buffer_solver = true;
        params.limit = 50;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);
        expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    "test_negative_coeff2"_test = [] {
        auto ctx = baryonyx::make_context(4);

        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff2.lp");

        baryonyx::solver_parameters params;
        // params.limit = 2;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);
        expect(result.affected_vars.names.size() +
                 result.variable_name.size() ==
               4);
    };

    "test_negative_coeff2_with_buffer"_test = [] {
        auto ctx = baryonyx::make_context(4);

        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff2.lp");

        baryonyx::solver_parameters params;
        params.use_buffer_solver = true;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);
        expect(result.affected_vars.names.size() +
                 result.variable_name.size() ==
               4);
    };

    "test_negative_coeff3"_test = [] {
        auto ctx = baryonyx::make_context(4);

        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff3.lp");
        expect(pb);

        baryonyx::solver_parameters params;
        params.limit = 10000;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);
        expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    "test_negative_coeff3_with_buffer"_test = [] {
        auto ctx = baryonyx::make_context(4);

        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff3.lp");
        expect(pb);

        baryonyx::solver_parameters params;
        params.limit = 10000;
        params.use_buffer_solver = true;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);
        expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    "test_negative_coeff4"_test = [] {
        auto ctx = baryonyx::make_context(4);

        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff4.lp");

        baryonyx::solver_parameters params;
        params.limit = 50;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);
        expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    "test_negative_coeff4_with_buffer"_test = [] {
        auto ctx = baryonyx::make_context(4);

        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/negative-coeff4.lp");

        baryonyx::solver_parameters params;
        params.limit = 50;
        params.use_buffer_solver = true;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);
        expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    "test_negative_coeff5"_test = [] {
        auto ctx = baryonyx::make_context(4);

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

        auto pb = baryonyx::make_problem(*ctx, iss);
        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);

        if (result)
            expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    "test_negative_coeff5_with_buffer"_test = [] {
        auto ctx = baryonyx::make_context(4);

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

        auto pb = baryonyx::make_problem(*ctx, iss);
        baryonyx::solver_parameters params;
        params.use_buffer_solver = true;
        baryonyx::context_set_solver_parameters(ctx, params);
        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);

        if (result)
            expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    "test_8_queens_puzzle_fixed_cost"_test = [] {
        auto ctx = baryonyx::make_context(4);

        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/8_queens_puzzle.lp");

        baryonyx::solver_parameters params;
        params.limit = -1;
        params.theta = 0.5;
        params.delta = 0.02;
        params.kappa_step = 0.01;
        params.kappa_max = 60.0;
        params.alpha = 1.0;
        params.w = 40;
        baryonyx::context_set_solver_parameters(ctx, params);

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
        expect(result.status == baryonyx::result_status::success);
        expect(!result.solutions.empty());

        for (int i = 0; i != 8; ++i) {
            for (int j = 0; j != 8; ++j)
                fmt::print("{} ",
                           result.solutions.back().variables[j * 8 + i]);

            fmt::print("\n");
        }

        expect(result.status == baryonyx::result_status::success);
        expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    "test_8_queens_puzzle_fixed_cost_with_buffer"_test = [] {
        auto ctx = baryonyx::make_context(4);

        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/8_queens_puzzle.lp");

        baryonyx::solver_parameters params;
        params.limit = -1;
        params.theta = 0.5;
        params.delta = 0.02;
        params.kappa_step = 0.01;
        params.kappa_max = 60.0;
        params.alpha = 1.0;
        params.w = 40;
        params.use_buffer_solver = true;
        baryonyx::context_set_solver_parameters(ctx, params);

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
        expect(result.status == baryonyx::result_status::success);
        expect(!result.solutions.empty());

        for (int i = 0; i != 8; ++i) {
            for (int j = 0; j != 8; ++j)
                fmt::print("{} ",
                           result.solutions.back().variables[j * 8 + i]);

            fmt::print("\n");
        }

        expect(result.status == baryonyx::result_status::success);
        expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    "test_8_queens_puzzle_random_cost"_test = [] {
        auto ctx = baryonyx::make_context(4);

        baryonyx::solver_parameters params;
        params.limit = -1;
        params.theta = 0.5;
        params.delta = 0.02;
        params.kappa_step = 0.01;
        params.kappa_max = 60.0;
        params.alpha = 1.0;
        params.w = 40;
        params.order =
          baryonyx::solver_parameters::constraint_order::infeasibility_decr;

        baryonyx::context_set_solver_parameters(ctx, params);

        for (int i{ 0 }, e{ 10 }; i != e; ++i) {
            auto pb =
              baryonyx::make_problem(ctx, EXAMPLES_DIR "/8_queens_puzzle.lp");

            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> dis(1, 100);

            for (auto& elem : pb.objective.elements)
                elem.factor = dis(gen);

            auto result = baryonyx::solve(ctx, pb);

            expect(result.status == baryonyx::result_status::success);
            expect(baryonyx::is_valid_solution(pb, result) == true);
        }
    };

    "test_8_queens_puzzle_random_cost_with_buffer"_test = [] {
        auto ctx = baryonyx::make_context(4);

        baryonyx::solver_parameters params;
        params.limit = -1;
        params.theta = 0.5;
        params.delta = 0.02;
        params.kappa_step = 0.01;
        params.kappa_max = 60.0;
        params.alpha = 1.0;
        params.w = 40;
        params.use_buffer_solver = true;
        params.order =
          baryonyx::solver_parameters::constraint_order::infeasibility_decr;

        baryonyx::context_set_solver_parameters(ctx, params);

        for (int i{ 0 }, e{ 10 }; i != e; ++i) {
            auto pb =
              baryonyx::make_problem(ctx, EXAMPLES_DIR "/8_queens_puzzle.lp");

            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<> dis(1, 100);

            for (auto& elem : pb.objective.elements)
                elem.factor = dis(gen);

            auto result = baryonyx::solve(ctx, pb);

            expect(result.status == baryonyx::result_status::success);
            expect(baryonyx::is_valid_solution(pb, result) == true);
        }
    };

    "test_qap"_test = [] {
        auto ctx = baryonyx::make_context(4);

        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/small4.lp");

        baryonyx::solver_parameters params;
        params.limit = -1;
        params.theta = 0.5;
        params.delta = 0.01;
        params.kappa_step = 10e-4;
        params.kappa_max = 10.0;
        params.alpha = 0.0;
        params.w = 20;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        if (result)
            expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    "test_qap_with_buffer"_test = [] {
        auto ctx = baryonyx::make_context(4);

        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/small4.lp");

        baryonyx::solver_parameters params;
        params.limit = -1;
        params.theta = 0.5;
        params.delta = 0.01;
        params.kappa_step = 10e-4;
        params.kappa_max = 10.0;
        params.alpha = 0.0;
        params.w = 20;
        params.use_buffer_solver = true;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        if (result)
            expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    "test_flat30_7"_test = [] {
        auto ctx = baryonyx::make_context(4);
        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/flat30-7.lp");

        baryonyx::solver_parameters params;
        params.limit = -1;
        params.delta = 0.001;
        params.kappa_min = 0.3;
        params.kappa_step = 1e-10;
        params.kappa_max = 1.0;
        params.order =
          baryonyx::solver_parameters::constraint_order::reversing;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);
        expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    //"test_flat30_7_with_buffer"_test = [] {
    //    auto ctx = baryonyx::make_context(4);
    //    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/flat30-7.lp");

    //    baryonyx::solver_parameters params;
    //    params.limit = -1;
    //    params.delta = 0.001;
    //    params.kappa_min = 0.;
    //    params.kappa_step = 1e-10;
    //    params.kappa_max = 1.0;
    //    params.order =
    //      baryonyx::solver_parameters::constraint_order::random_sorting;
    //    params.use_buffer_solver = true;
    //    baryonyx::context_set_solver_parameters(ctx, params);

    //    auto result = baryonyx::solve(ctx, pb);

    //    expect(result.status == baryonyx::result_status::success);
    //    expect(baryonyx::is_valid_solution(pb, result) == true);
    //};

    "test_uf50_0448"_test = [] {
        auto ctx = baryonyx::make_context(4);
        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/uf50-0448.lp");

        baryonyx::solver_parameters params;
        params.limit = -1;
        params.theta = 0.5;
        params.delta = 1.0;
        params.kappa_min = 0.1;
        params.kappa_step = 1e-17;
        params.kappa_max = 1.0;
        params.alpha = 2.0;
        params.w = 60;
        params.order =
          baryonyx::solver_parameters::constraint_order::random_sorting;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);
        expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    //"test_uf50_0448_with_buffer"_test = [] {
    //    auto ctx = baryonyx::make_context(4);
    //    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/uf50-0448.lp");

    //    baryonyx::solver_parameters params;
    //    params.limit = -1;
    //    params.theta = 0.5;
    //    params.delta = 1.0;
    //    params.kappa_min = 0.1;
    //    params.kappa_step = 1e-17;
    //    params.kappa_max = 1.0;
    //    params.alpha = 2.0;
    //    params.w = 60;
    //    params.use_buffer_solver = true;
    //    params.order =
    //      baryonyx::solver_parameters::constraint_order::random_sorting;
    //    baryonyx::context_set_solver_parameters(ctx, params);

    //    auto result = baryonyx::solve(ctx, pb);

    //    expect(result.status == baryonyx::result_status::success);
    //    expect(baryonyx::is_valid_solution(pb, result) == true);
    //};

    "test_aim_50_1_6_yes1_2"_test = [] {
        auto ctx = baryonyx::make_context(4);

        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/aim-50-1_6-yes1-2.lp");

        baryonyx::solver_parameters params;
        params.limit = -1;
        params.theta = 0.6;
        params.delta = 0.01;
        params.kappa_step = 2 * 10e-4;
        params.kappa_max = 100.0;
        params.alpha = 1.0;
        params.w = 20;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);
        expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    "test_aim_50_1_6_yes1_2_with_buffer"_test = [] {
        auto ctx = baryonyx::make_context(4);

        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/aim-50-1_6-yes1-2.lp");

        baryonyx::solver_parameters params;
        params.limit = -1;
        params.theta = 0.6;
        params.delta = 0.01;
        params.kappa_step = 2 * 10e-4;
        params.kappa_max = 100.0;
        params.alpha = 1.0;
        params.w = 20;
        params.use_buffer_solver = true;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::solve(ctx, pb);

        expect(result.status == baryonyx::result_status::success);
        expect(baryonyx::is_valid_solution(pb, result) == true);
    };

    "test_Z_coefficient_1"_test = [] {
        auto ctx = baryonyx::make_context(4);

        {
            const char* str_pb = "minimize\n"
                                 "15 x1 + 19 x2 + 13 x3 + 12 x4\n"
                                 "Subject to:\n"
                                 "2 x1 + 1 x2 + 3 x3 + 2 x4 <= 3\n"
                                 "Binaries\n"
                                 "x1 x2 x3 x4\n"
                                 "End\n";

            std::istringstream iss(str_pb);

            auto pb = baryonyx::make_problem(*ctx, iss);
            auto result = baryonyx::solve(ctx, pb);

            expect(result.status == baryonyx::result_status::success);

            if (result)
                expect(baryonyx::is_valid_solution(pb, result) == true);
        }

        {
            const char* str_pb = "minimize\n"
                                 "Subject to:\n"
                                 "2a + 3b -5c + 7d <= 0\n"
                                 "-2b + 2c >= 1\n"
                                 "7a + 7d <= 7\n"
                                 "Binaries\n"
                                 "a b c d\n"
                                 "End\n";

            std::istringstream iss(str_pb);

            auto pb = baryonyx::make_problem(*ctx, iss);
            auto result = baryonyx::solve(ctx, pb);

            expect(result.status == baryonyx::result_status::success);

            if (result)
                expect(baryonyx::is_valid_solution(pb, result) == true);
        }
    };

    //"test_bibd1n"_test = [] {
    //    auto ctx = baryonyx::make_context(4);

    //    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/bibd1n.lp");

    //    baryonyx::solver_parameters params;
    //    params.limit = -1;
    //    params.theta = 0.6;
    //    params.delta = 1e-7;
    //    params.kappa_step = 1e-17;
    //    params.kappa_max = 0.05;
    //    params.kappa_max = 1.0;
    //    params.alpha = 1.0;
    //    params.w = 60;
    //    params.use_buffer_solver = true;
    //    baryonyx::context_set_solver_parameters(ctx, params);

    //    auto result = baryonyx::solve(ctx, pb);

    //    expect(result.status == baryonyx::result_status::success);
    //    expect(baryonyx::is_valid_solution(pb, result) == true);
    //};
}
