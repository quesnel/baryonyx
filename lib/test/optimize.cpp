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

#include <baryonyx/core-compare>
#include <baryonyx/core-out>
#include <baryonyx/core-test>
#include <baryonyx/core>

#include <fmt/printf.h>

void
test_qap(const baryonyx::context_ptr& ctx)
{
    auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/small4.lp");

    context_set_parameter(ctx, "limit", -1);
    context_set_parameter(ctx, "theta", 0.5);
    context_set_parameter(ctx, "delta", 0.2);
    context_set_parameter(ctx, "kappa-step", 10e-4);
    context_set_parameter(ctx, "kappa-max", 10.0);
    context_set_parameter(ctx, "alpha", 0.0);
    context_set_parameter(ctx, "w", 20);
    context_set_parameter(ctx, "time-limit", 40.0);
    context_set_parameter(ctx, "pushing-k-factor", 0.9);
    context_set_parameter(ctx, "pushes-limit", 50);
    context_set_parameter(ctx, "pushing-objective-amplifier", 10);
    context_set_parameter(ctx, "pushing-iteration-limit", 50);
    context_set_parameter(ctx, "thread", 2);

    auto result = baryonyx::optimize(ctx, pb);

    Ensures(result.status == baryonyx::result_status::success);
    if (result.status == baryonyx::result_status::success)
        Ensures(result.solutions.back().value == 790.0);

    if (result.status == baryonyx::result_status::success) {
        pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/small4.lp");

        Ensures(baryonyx::is_valid_solution(
                  pb, result.solutions.back().variables) == true);
        Ensures(baryonyx::compute_solution(
                  pb, result.solutions.back().variables) == 790.0);
    }
}

void
test_n_queens_problem(const baryonyx::context_ptr& ctx)
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

    context_set_parameter(ctx, "limit", 100000);
    context_set_parameter(ctx, "theta", 0.5);
    context_set_parameter(ctx, "delta", 1.0);
    context_set_parameter(ctx, "kappa-min", 0.30);
    context_set_parameter(ctx, "kappa-step", 1e-2);
    context_set_parameter(ctx, "kappa-max", 100.0);
    context_set_parameter(ctx, "alpha", 1.0);
    context_set_parameter(ctx, "w", 60);
    context_set_parameter(
      ctx, "constraint-order", std::string("random-sorting"));
    context_set_parameter(ctx, "time-limit", 20.0);
    context_set_parameter(ctx, "pushing-k-factor", 0.9);
    context_set_parameter(ctx, "pushes-limit", 50);
    context_set_parameter(ctx, "pushing-objective-amplifier", 10);
    context_set_parameter(ctx, "pushing-iteration-limit", 10);

    for (std::size_t i{ 0 }; i != valid_solutions.size(); ++i) {
        std::string filepath{ EXAMPLES_DIR "/n-queens/n-queens-problem-" };
        filepath += std::to_string(i);
        filepath += ".lp";

        auto pb = baryonyx::make_problem(ctx, filepath);
        auto result = baryonyx::optimize(ctx, pb);

        valid_solutions[i] = (result.remaining_constraints == 0);
        if (valid_solutions[i])
            solutions[i] = result.solutions.back().value;
    }

    auto all_found =
      std::accumulate(valid_solutions.cbegin(),
                      valid_solutions.cend(),
                      static_cast<std::size_t>(0),
                      [](std::size_t cumul, const auto& elem) -> std::size_t {
                          return elem ? cumul + 1 : cumul;
                      });

    double mean_distance{ 0 };

    for (std::size_t i{ 0 }, e{ solutions.size() }; i != e; ++i) {
        double distance =
          ((cplex_solutions[i] - solutions[i]) / cplex_solutions[i]) * 100.0;

        fmt::print("{}: {} {} {} {}\n",
                   i,
                   valid_solutions[i],
                   solutions[i],
                   cplex_solutions[i],
                   distance);

        mean_distance += distance;
    }

    fmt::print("Optimum means: {}\n",
               mean_distance / static_cast<double>(solutions.size()));

    Ensures(all_found == valid_solutions.size());
}

int
main(int /* argc */, char* /* argv */ [])
{
    auto ctx = baryonyx::make_context();

    test_qap(ctx);
    test_n_queens_problem(ctx);

    return unit_test::report_errors();
}
