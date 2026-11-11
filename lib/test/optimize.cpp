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

#include <baryonyx/core-compare>
#include <baryonyx/core-test>
#include <baryonyx/core>

#include <boost/ut.hpp>

#include <fmt/core.h>

#include <chrono>
#include <fstream>
#include <thread>
#include <sstream>
#include <string>

int
main()
{
    using namespace boost::ut;
    using optimizer_type = baryonyx::solver_parameters::optimizer_type;

    constexpr optimizer_type optimizers[] = { optimizer_type::genetic,
                                              optimizer_type::restart };
    const auto name = [](optimizer_type o) {
        return o == optimizer_type::genetic ? "genetic" : "restart";
    };

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

    "test_optimize_returns_its_solution"_test = [&] {
        for (const auto optimizer : optimizers) {
            for (const int seed : { 1, 2, 3 }) {
                auto ctx = baryonyx::make_context(0);
                auto pb = baryonyx::make_problem(
                  ctx, EXAMPLES_DIR "/assignment_problem_1.lp");

                baryonyx::solver_parameters params;
                params.time_limit = 0.5;
                params.seed = seed;
                params.thread = 4;
                params.optimizer = optimizer;
                baryonyx::context_set_solver_parameters(ctx, params);

                auto result = baryonyx::optimize(ctx, pb);
                expect((result.status ==
                        baryonyx::result_status::success) >>
                       fatal)
                  << name(optimizer);
                expect(baryonyx::is_valid_solution(pb, result));
                expect(baryonyx::compute_solution(pb, result) ==
                       result.solutions.back().value);
                expect(result.solutions.back().value == 15.0)
                  << name(optimizer);
            }
        }
    };

    "test_optimize_maximization"_test = [&] {
        const char* str_pb = "maximize\n"
                             "-3 x1 - 8 x2 - 6 x3 - 4 x4 - 4 x5 - 2 x6 - 3 x7"
                             " - 2 x8 - 2 x9 - 7 x10 - 10 x11 - 8 x12 - 6 x13"
                             " - 8 x14 - 10 x15 - 5 x16\n"
                             "subject to:\n"
                             "x1 + x2 + x3 + x4 = 1\n"
                             "x5 + x6 + x7 + x8 = 1\n"
                             "x9 + x10 + x11 + x12 = 1\n"
                             "x13 + x14 + x15 + x16 = 1\n"
                             "x1 + x5 + x9 + x13 = 1\n"
                             "x2 + x6 + x10 + x14 = 1\n"
                             "x3 + x7 + x11 + x15 = 1\n"
                             "x4 + x8 + x12 + x16 = 1\n"
                             "binary\n"
                             "x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14"
                             " x15 x16\n"
                             "end\n";

        for (const auto optimizer : optimizers) {
            auto ctx = baryonyx::make_context(0);
            std::istringstream iss(str_pb);
            auto pb = baryonyx::make_problem(*ctx, iss);
            const auto copy = pb;

            baryonyx::solver_parameters params;
            params.time_limit = 0.5;
            params.seed = 1;
            params.thread = 2;
            params.optimizer = optimizer;
            baryonyx::context_set_solver_parameters(ctx, params);

            auto result = baryonyx::optimize(ctx, pb);
            expect((result.status == baryonyx::result_status::success) >>
                   fatal);
            expect(baryonyx::is_valid_solution(copy, result));
            expect(result.solutions.back().value == -15.0) << name(optimizer);
        }
    };

    "test_tuned_and_lns_optimizers"_test = [] {
        // Both run within the time limit and return their best result, with
        // the sign of the objective of the problem.
        using optimizer_type = baryonyx::solver_parameters::optimizer_type;

        struct test_case
        {
            optimizer_type optimizer;
            const char* sense;
            double expected;
        };

        constexpr test_case cases[] = {
            { optimizer_type::tuned, "minimize", 15.0 },
            { optimizer_type::tuned, "maximize", -15.0 },
            { optimizer_type::lns, "minimize", 15.0 },
            { optimizer_type::lns, "maximize", -15.0 },
        };

        for (const auto& [optimizer, sense, expected] : cases) {
            auto ctx = baryonyx::make_context(0);
            auto pb = baryonyx::make_problem(
              ctx, EXAMPLES_DIR "/assignment_problem_1.lp");

            if (std::string_view(sense) == "maximize") {
                pb.type = baryonyx::objective_function_type::maximize;
                for (auto& elem : pb.objective.elements)
                    elem.factor = -elem.factor;
            }

            const auto copy = pb;

            baryonyx::solver_parameters params;
            params.time_limit = 1.0;
            params.seed = 1;
            params.thread = 2;
            params.optimizer = optimizer;
            baryonyx::context_set_solver_parameters(ctx, params);

            const auto start = std::chrono::steady_clock::now();
            auto result = baryonyx::optimize(ctx, pb);
            const auto elapsed = std::chrono::duration<double>(
                                   std::chrono::steady_clock::now() - start)
                                   .count();

            expect((result.status == baryonyx::result_status::success) >>
                   fatal)
              << sense;
            expect(baryonyx::is_valid_solution(copy, result));
            expect(result.solutions.back().value == expected) << sense;
            expect(elapsed < 2.0) << sense;
        }
    };

    "test_request_stop"_test = [] {
        // Without time limit the optimization only ends on request_stop(),
        // with its best solution; the next optimization is not stopped.
        auto ctx = baryonyx::make_context(0);
        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/assignment_problem_1.lp");

        baryonyx::solver_parameters params;
        params.time_limit = 0.0;
        params.seed = 1;
        params.thread = 2;
        baryonyx::context_set_solver_parameters(ctx, params);

        const auto elapsed_since = [](auto start) {
            return std::chrono::duration<double>(
                     std::chrono::steady_clock::now() - start)
              .count();
        };

        auto start = std::chrono::steady_clock::now();
        std::thread stopper([] {
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
            baryonyx::request_stop();
        });
        auto result = baryonyx::optimize(ctx, pb);
        stopper.join();

        expect((result.status == baryonyx::result_status::success) >> fatal);
        expect(baryonyx::is_valid_solution(pb, result));
        expect(elapsed_since(start) < 2.0);

        params.time_limit = 0.5;
        baryonyx::context_set_solver_parameters(ctx, params);

        start = std::chrono::steady_clock::now();
        result = baryonyx::optimize(ctx, pb);
        expect(elapsed_since(start) > 0.4);
        expect(result.solutions.back().value == 15.0);
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
        params.thread = 2;
        baryonyx::context_set_solver_parameters(ctx, params);

        auto result = baryonyx::optimize(ctx, pb);
        expect((result.status == baryonyx::result_status::success) >> fatal);

        pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/small4.lp");
        expect(baryonyx::is_valid_solution(pb, result));
        expect(baryonyx::compute_solution(pb, result) == 790.0);
        expect(result.solutions.back().value == 790.0);
    };

    "test_n_queens_problem"_test = [&] {
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

        for (const auto optimizer : optimizers) {
            baryonyx::solver_parameters params;
            params.optimizer = optimizer;
            params.time_limit = 0.5;
            params.limit = 5000;
            params.theta = 0.5;
            params.delta = 1.0;
            params.alpha = 1.0;
            params.w = 60;
            params.thread = 2;
            params.order =
              baryonyx::solver_parameters::constraint_order::random_sorting;
            baryonyx::context_set_solver_parameters(ctx, params);

            double total_gap = 0.0;

            for (std::size_t i{ 0 }; i != cplex_solutions.size(); ++i) {
                const auto filepath = fmt::format(
                  EXAMPLES_DIR "/n-queens/n-queens-problem-{}.lp", i);

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
            fmt::print("n-queens mean gap to CPLEX ({}): {:.2f}%\n",
                       name(optimizer),
                       mean_gap);
            expect(mean_gap <= 1.0) << name(optimizer);
        }
    };
}
