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

#include <baryonyx/core-compare>
#include <baryonyx/core-test>
#include <baryonyx/core>

#include <boost/ut.hpp>

#include <fmt/core.h>

#include <fstream>

int
main()
{
    using namespace boost::ut;

    "test_bibd1n"_test = [] {
        auto ctx = baryonyx::make_context(4);
        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/bibd1n.lp");

        baryonyx::solver_parameters params;
        params.delta = 1e-2;
        params.time_limit = 1.0;
        params.limit = 5000;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::optimize(ctx, pb);
        expect(result.status != baryonyx::result_status::internal_error);
        expect((result.status == baryonyx::result_status::success) ==
               (result.remaining_constraints == 0));

        if (result.status == baryonyx::result_status::success) {
            pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/bibd1n.lp");
            expect(baryonyx::is_valid_solution(pb, result));
        }
    };

    "test_optimize_returns_its_solution"_test = [] {
        for (const int seed : { 1, 2, 3 }) {
            auto ctx = baryonyx::make_context(0);
            auto pb = baryonyx::make_problem(
              ctx, EXAMPLES_DIR "/assignment_problem_1.lp");

            baryonyx::solver_parameters params;
            params.time_limit = 0.5;
            params.seed = seed;
            params.thread = 4;
            baryonyx::context_set_solver_parameters(ctx, params);

            auto result = baryonyx::optimize(ctx, pb);
            expect((result.status == baryonyx::result_status::success) >>
                   fatal);
            expect(baryonyx::is_valid_solution(pb, result));
            expect(baryonyx::compute_solution(pb, result) ==
                   result.solutions.back().value);
        }
    };

    "test_qap"_test = [] {
        auto ctx = baryonyx::make_context(4);
        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/small4.lp");

        baryonyx::solver_parameters params;
        params.time_limit = 1.0;
        params.limit = 5000;
        params.theta = 0.5;
        params.delta = 0.2;
        params.alpha = 0.0;
        params.w = 20;
        params.pushing_k_factor = 0.9;
        params.pushes_limit = 50;
        params.pushing_objective_amplifier = 10;
        params.pushing_iteration_limit = 50;
        params.thread = 2;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::optimize(ctx, pb);
        expect((result.status == baryonyx::result_status::success) >> fatal);

        pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/small4.lp");
        expect(baryonyx::is_valid_solution(pb, result));
        expect(baryonyx::compute_solution(pb, result) == 790.0);
        expect(result.solutions.back().value == 790.0);
    };

    "test_n_queens_problem"_test = [] {
        auto ctx = baryonyx::make_context(4);
        std::vector<double> cplex_solutions(30, 0.0);

        {
            // Tries to read the cplex solution files produced by
            // CPLEX 12.7.0.0 and the `script.sh'
            // files. If an error occured, the test fails and returns.*/

            std::ifstream ifs{ EXAMPLES_DIR "/n-queens/solutions.txt" };

            expect(ifs.is_open());
            if (!ifs.is_open())
                return;

            for (auto& elem : cplex_solutions)
                ifs >> elem;

            expect(ifs.good());
            if (!ifs.good())
                return;
        }

        baryonyx::solver_parameters params;
        params.time_limit = 0.5;
        params.limit = 5000;
        params.theta = 0.5;
        params.delta = 1.0;
        params.alpha = 1.0;
        params.w = 60;
        params.pushing_k_factor = 0.9;
        params.pushes_limit = 50;
        params.pushing_objective_amplifier = 10;
        params.pushing_iteration_limit = 10;
        params.thread = 2;
        params.order =
          baryonyx::solver_parameters::constraint_order::random_sorting;
        baryonyx::context_set_solver_parameters(ctx, params);

        double total_gap = 0.0;

        for (std::size_t i{ 0 }; i != cplex_solutions.size(); ++i) {
            std::string filepath{ EXAMPLES_DIR "/n-queens/n-queens-problem-" };
            filepath += std::to_string(i);
            filepath += ".lp";

            auto pb = baryonyx::make_problem(ctx, filepath);
            auto result = baryonyx::optimize(ctx, pb);
            expect((result.status == baryonyx::result_status::success) >>
                   fatal)
              << filepath;

            pb = baryonyx::make_problem(ctx, filepath);
            const auto value = result.solutions.back().value;
            expect(baryonyx::is_valid_solution(pb, result)) << filepath;
            expect(baryonyx::compute_solution(pb, result) == value)
              << filepath;
            expect(value >= cplex_solutions[i]) << filepath;

            total_gap += (value - cplex_solutions[i]) / cplex_solutions[i];
        }

        const auto mean_gap =
          100.0 * total_gap / static_cast<double>(cplex_solutions.size());
        fmt::print("n-queens mean gap to CPLEX: {:.2f}%\n", mean_gap);
        expect(mean_gap <= 1.0);
    };
}
