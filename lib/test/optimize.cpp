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
test_qap(std::shared_ptr<lp::context> ctx)
{
    lp::result result;

    {
        auto pb = lp::make_problem(ctx, EXAMPLES_DIR "/small4.lp");

        ctx->set_parameter("limit", 10'000'000l);
        ctx->set_parameter("theta", 0.5);
        ctx->set_parameter("delta", 0.2);
        ctx->set_parameter("kappa-step", 10e-4);
        ctx->set_parameter("kappa-max", 10.0);
        ctx->set_parameter("alpha", 0.0);
        ctx->set_parameter("w", 20l);
        ctx->set_parameter("time-limit", 40.0);
        ctx->set_parameter("pushing-k-factor", 0.9);
        ctx->set_parameter("pushes-limit", 50l);
        ctx->set_parameter("pushing-objective-amplifier", 10l);
        ctx->set_parameter("pushing-iteration-limit", 50l);
        ctx->set_parameter("thread", 2l);

        result = lp::optimize(ctx, pb);

        Ensures(result.status == lp::result_status::success);
        if (result.status == lp::result_status::success)
            Ensures(result.value == 790.0);
    }

    {
        if (result.status == lp::result_status::success) {
            auto pb = lp::make_problem(ctx, EXAMPLES_DIR "/small4.lp");

            Ensures(lp::is_valid_solution(pb, result.variable_value) == true);
            Ensures(lp::compute_solution(pb, result.variable_value) == 790.0);
        }
    }
}

void
test_n_queens_problem(std::shared_ptr<lp::context> ctx)
{
    std::vector<bool> valid_solutions(30, false);
    std::vector<double> solutions(30, 0.0);
    std::vector<double> cplex_solutions(30, 0.0);

    {
        // Tries to read the cplex solution files produced by CPLEX 12.7.0.0
        // and the `script.sh' `n-queens-problem.commands' files. If an error
        // occured, the test fails and returns.*/

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

    ctx->set_parameter("limit", 100000l);
    ctx->set_parameter("theta", 0.5);
    ctx->set_parameter("delta", 1.0);
    ctx->set_parameter("kappa-min", 0.30);
    ctx->set_parameter("kappa-step", 1e-2);
    ctx->set_parameter("kappa-max", 100.0);
    ctx->set_parameter("alpha", 1.0);
    ctx->set_parameter("w", 60l);
    ctx->set_parameter("constraint-order", std::string("random-sorting"));
    ctx->set_parameter("time-limit", 20.0);
    ctx->set_parameter("pushing-k-factor", 0.9);
    ctx->set_parameter("pushes-limit", 50l);
    ctx->set_parameter("pushing-objective-amplifier", 10l);
    ctx->set_parameter("pushing-iteration-limit", 10l);

    for (std::size_t i{ 0 }; i != valid_solutions.size(); ++i) {
        std::string filepath{ EXAMPLES_DIR "/n-queens/n-queens-problem-" };
        filepath += std::to_string(i);
        filepath += ".lp";

        auto pb = lp::make_problem(ctx, filepath);
        auto result = lp::optimize(ctx, pb);

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
    auto ctx = std::make_shared<lp::context>();
    ctx->set_standard_stream_logger();

    test_qap(ctx);
    test_n_queens_problem(ctx);

    return unit_test::report_errors();
}
