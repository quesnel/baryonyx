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
#include <lpcore-compare>
#include <lpcore-out>
#include <lpcore>
#include <map>
#include <numeric>
#include <random>
#include <sstream>

void
test_assignment_problem()
{
    auto pb = lp::make_problem(EXAMPLES_DIR "/assignment_problem_1.lp");

    std::map<std::string, lp::parameter> params;
    params["kappa"] = 0.0001;
    params["theta"] = 0.5;
    params["delta"] = 0.5;
    params["limit"] = 50l;

    auto result = lp::solve(pb, params);

    Ensures(result.solution_found == true);
}

void
test_assignment_problem_random_coast()
{
    auto pb = lp::make_problem(EXAMPLES_DIR "/assignment_problem_1.lp");

    std::map<std::string, lp::parameter> params;
    params["kappa"] = 0.0001;
    params["theta"] = 0.5;
    params["delta"] = 0.5;
    params["limit"] = 50l;

    for (int i{ 0 }, e{ 10 }; i != e; ++i) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(1, 100);

        for (auto& elem : pb.objective.elements)
            elem.factor = dis(gen);

        auto result = lp::solve(pb, params);

        Ensures(result.solution_found == true);
    }
}

void
test_inequality()
{
    auto pb = lp::make_problem(EXAMPLES_DIR "/inequality.lp");

    std::map<std::string, lp::parameter> params;
    params["limit"] = 50l;

    auto result = lp::solve(pb, params);
    std::cout << result << '\n';

    Ensures(result.solution_found == true);
}

void
test_inequality_1()
{
    auto pb = lp::make_problem(EXAMPLES_DIR "/inequality0.lp");

    std::map<std::string, lp::parameter> params;
    params["limit"] = 50l;

    auto result = lp::solve(pb, params);
    std::cout << result << '\n';

    Ensures(result.solution_found == true);
    Ensures(result.variable_value[0] == 1);
    Ensures(result.variable_value[1] == 0);
    Ensures(result.variable_value[2] == 0);
    Ensures(result.variable_value[3] == 1);
}

void
test_8_queens_puzzle_fixed_cost()
{
    auto pb = lp::make_problem(EXAMPLES_DIR "/8_queens_puzzle.lp");

    std::map<std::string, lp::parameter> params;
    params["kappa"] = 0.1;
    params["theta"] = 0.5;
    params["delta"] = 0.5;
    params["limit"] = 100l;

    std::vector<int> cost{ 25, 89, 12, 22, 84, 3,  61, 14, 93, 97, 68, 5,  51,
                           72, 96, 80, 13, 38, 81, 48, 70, 50, 66, 68, 30, 97,
                           79, 4,  41, 44, 47, 62, 60, 11, 18, 44, 57, 24, 7,
                           11, 66, 87, 9,  17, 27, 60, 95, 45, 94, 47, 60, 87,
                           79, 53, 81, 52, 91, 53, 57, 8,  63, 78, 1,  8 };

    std::size_t i{ 0 };

    for (auto& elem : pb.objective.elements)
        elem.factor = cost[i++];

    auto result = lp::solve(pb, params);

    std::cout << result << '\n';

    for (int i = 0; i != 8; ++i) {
        for (int j = 0; j != 8; ++j) {
            std::cout << result.variable_value[j * 8 + i] << ' ';
        }
        std::cout << '\n';
    }

    Ensures(result.solution_found == true);
}

void
test_8_queens_puzzle_random_cost()
{
    auto pb = lp::make_problem(EXAMPLES_DIR "/8_queens_puzzle.lp");

    std::map<std::string, lp::parameter> params;
    params["limit"] = 1'000l;

    for (int i{ 0 }, e{ 10 }; i != e; ++i) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(1, 100);

        for (auto& elem : pb.objective.elements)
            elem.factor = dis(gen);

        auto result = lp::solve(pb, params);

        Ensures(result.solution_found == true);
    }
}

void
test_qap()
{
    auto pb = lp::make_problem(EXAMPLES_DIR "/small4.lp");

    std::map<std::string, lp::parameter> params;
    params["limit"] = 10'000l;

    auto result = lp::solve(pb, params);

    std::cout << result << '\n';
}

void
test_verger_5_5()
{
    auto pb = lp::make_problem(EXAMPLES_DIR "/verger_5_5.lp");

    std::map<std::string, lp::parameter> params;
    params["limit"] = 1'000l;
    params["w"] = 3l;

    auto result = lp::solve(pb, params);

    std::cout << result << '\n';
}

int
main(int /* argc */, char* /* argv */ [])
{
    // test_assignment_problem();
    // test_assignment_problem_random_coast();
    // test_inequality();
    // test_inequality_1();
    // test_8_queens_puzzle_fixed_cost();
    // test_8_queens_puzzle_random_cost();
    test_qap();
    // test_verger_5_5();

    return unit_test::report_errors();
}
