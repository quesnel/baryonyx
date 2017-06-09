/* Copyright (C) 2017 INRA
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
#include <lptest>
#include <map>
#include <numeric>
#include <random>
#include <sstream>

void
test_qap()
{
    lp::result result;

    {
        auto pb = lp::make_problem(EXAMPLES_DIR "/small4.lp");

        std::map<std::string, lp::parameter> params;
        params["limit"] = 1000000l;
        params["theta"] = 0.5;
        params["delta"] = 0.1;
        params["kappa-step"] = 1e-3;
        params["kappa-max"] = 10.0;
        params["alpha"] = 1.0;
        params["w"] = 20l;

        params["time-limit"] = 20.0;
        params["pushing-k-factor"] = 0.9;
        params["pushes-limit"] = 50l;
        params["pushing-objective-amplifier"] = 10l;
        params["pushing-iteration-limit"] = 50l;

        params["thread"] = 2l;

        result = lp::optimize(pb, params);

        Ensures(result.status == lp::result_status::success);
        if (result.status == lp::result_status::success)
            Ensures(result.value == 790.0);
    }

    {
        if (result.status == lp::result_status::success) {
            auto pb = lp::make_problem(EXAMPLES_DIR "/small4.lp");

            Ensures(lp::is_valid_solution(pb, result.variable_value) == true);
            Ensures(lp::compute_solution(pb, result.variable_value) == 790.0);
        }
    }
}

void
test_n_queens_problem()
{
    std::vector<bool> valid_solutions(30, false);
    std::vector<double> solutions(30, 0.0);
    std::vector<double> cplex_solutions(30, 0.0);

    { /* Tries to read the cplex solution files produced by CPLEX 12.7.0.0
         and the `script.sh' `n-queens-problem.commands' files. If an
         error occured, the test fails and returns. */

        std::ifstream ifs{ EXAMPLES_DIR "/n-queens/solutions.txt" };

        Ensures(ifs.is_open());
        if (not ifs.is_open())
            return;

        for (auto& elem : cplex_solutions)
            ifs >> elem;

        Ensures(ifs.good());
        if (not ifs.good())
            return;
    }

    std::map<std::string, lp::parameter> params;
    params["limit"] = 100000l;
    params["theta"] = 0.5;
    params["delta"] = 0.000001;
    params["kappa-min"] = 0.30;
    params["kappa-step"] = 0.1;
    params["kappa-max"] = 100.0;
    params["alpha"] = 1.0;
    params["w"] = 60l;
    // params["constraint-order"] = std::string("infeasibility-incr");
    // params["constraint-order"] = std::string("adaptative");
    // params["constraint-order"] = std::string("infeasibility-decr");
    // params["constraint-order"] = std::string("none");
    // params["constraint-order"] = std::string("reversing");
    params["constraint-order"] = std::string("random-sorting");
    params["time-limit"] = 20.0;
    params["pushing-k-factor"] = 0.9;
    params["pushes-limit"] = 50l;
    params["pushing-objective-amplifier"] = 10l;
    params["pushing-iteration-limit"] = 10l;

    for (std::size_t i{ 0 }; i != valid_solutions.size(); ++i) {
        std::string filepath{ EXAMPLES_DIR "/n-queens/n-queens-problem-" };
        filepath += std::to_string(i);
        filepath += ".lp";

        std::ifstream ifs(filepath);
        Ensures(ifs.is_open());
        if (not ifs.is_open())
            return;

        auto pb = lp::make_problem(ifs);
        auto result = lp::optimize(pb, params);

        valid_solutions[i] = (result.remaining_constraints == 0);
        if (valid_solutions[i])
            solutions[i] = result.value;
    }

    auto all_found = std::accumulate(
      valid_solutions.cbegin(),
      valid_solutions.cend(),
      std::size_t{ 0 },
      [](int cumul, const auto& elem) { return elem ? cumul + 1 : cumul; });

    double mean_distance{ 0 };

    for (std::size_t i{ 0 }, e{ solutions.size() }; i != e; ++i) {
        double distance =
          ((cplex_solutions[i] - solutions[i]) / cplex_solutions[i]) * 100.0;
        printf("%zu: %s %f %f %f\n",
               i,
               (valid_solutions[i] ? "true" : "false"),
               solutions[i],
               cplex_solutions[i],
               distance);

        mean_distance += distance;
    }

    printf("Optimum means: %f\n", mean_distance / solutions.size());

    Ensures(all_found == valid_solutions.size());
}

int
main(int /* argc */, char* /* argv */ [])
{
    test_qap();
    test_n_queens_problem();

    return unit_test::report_errors();
}
