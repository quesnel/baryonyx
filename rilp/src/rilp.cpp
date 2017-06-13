/*
 * Copyright (C) 2017 INRA
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in he Software without restriction, including
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

#include <Rcpp.h>
#include <lpcore>

using namespace Rcpp;

/*
 * Convert the lp::result_status enumeration to an integer reusabe in R
 * code.
 */
constexpr int status_get(lp::result_status status) noexcept;

//' Tries to solve the 01 linear programming problem.
//'
//' @return A list:
//' - Integer: 0: internal_error, 1: success, 2: time limit reached, 3:
//'   kappa max reched and 4: loop limit reached.
//' - Real: the solution found.
//' - Real: the duration to found a solution.
//' - Integer: the constraints remaining if a limit is reached.
//'
//' @useDynLib rilp
//' @importFrom Rcpp sourceCpp
//'
//' @export
// [[Rcpp::export]]
List
solve_01lp_problem(std::string file_path,
                   long int limit = 1000,
                   double theta = 0.5,
                   double delta = 1e-4,
                   std::string constraint_order = "none",
                   double kappa_min = 0.1,
                   double kappa_step = 1e-4,
                   double kappa_max = 1.0,
                   double alpha = 1.0,
                   long int w = 60,
                   double time_limit = 10.0,
                   double pushing_k_factor = 0.9,
                   long int pushes_limit = 50,
                   double pushing_objective_amplifier = 10,
                   long int pushing_iteration_limit = 10,
                   long int thread = 1) noexcept
{
    try {
        auto pb = lp::make_problem(file_path);

        std::map<std::string, lp::parameter> params;
        params["limit"] = limit, params["theta"] = theta;
        params["delta"] = delta;
        params["kappa-min"] = kappa_min;
        params["kappa-step"] = kappa_step;
        params["kappa-max"] = kappa_max;
        params["alpha"] = alpha;
        params["constraint-order"] = constraint_order;
        params["w"] = w;

        params["time-limit"] = time_limit;
        params["pushing-k-factor"] = pushing_k_factor;
        params["pushes-limit"] = pushes_limit;
        params["pushing-objective-amplifier"] = pushing_objective_amplifier;
        params["pushing-iteration-limit"] = pushing_iteration_limit;

        params["thread"] = thread;

        auto result = lp::solve(pb, params);

        return List::create(status_get(result.status),
                            result.value,
                            result.duration,
                            result.remaining_constraints);

    } catch (const std::bad_alloc& e) {
        Rcpp::Rcout << "lp: " << e.what() << '\n';
    } catch (const std::exception& e) {
        Rcpp::Rcout << "lp: " << e.what() << '\n';
    } catch (...) {
        Rcpp::Rcout << "lp: unknown error\n";
    }

    return List::create();
}

constexpr int
status_get(lp::result_status status) noexcept
{
    switch (status) {
    case lp::result_status::success:
        return 1;
    case lp::result_status::time_limit_reached:
        return 2;
    case lp::result_status::kappa_max_reached:
        return 3;
        break;
    case lp::result_status::limit_reached:
        return 4;
        break;
    default:
        return 0;
        break;
    }
}
