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

//' Tries to solve the 01 linear programming problem.
//'
//' @param constraint_order: 0-none, 1-reversing, 2-random-sorting,
//' 3-infeasibility-decr, 4- infeasibility-incr
//'
//' @return A list of two value: The remaining constraint and the value of
//'   the solution (Integer, Real):
//'   - Integer: The number of remaining constraints. 0 means a solution
//'      found, greater than 0 means remaining constraints and value is
//'      0. -1 means an error occured during solving. Use the @c `lp`
//'      program to check out the error.
//'   - Real: The value of the solution found (if the remaining
//'     constraint equal 0).
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
                   long int constraint_order = 0,
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
                   long int seed = -1,
                   long int thread = 1) noexcept
{
    try {
        auto pb = lp::make_problem(file_path);

        std::map<std::string, lp::parameter> params;
        params["limit"] = limit;
        params["theta"] = theta;
        params["delta"] = delta;

        params["kappa-min"] = kappa_min;
        params["kappa-step"] = kappa_step;
        params["kappa-max"] = kappa_max;
        params["alpha"] = alpha;
        params["w"] = w;

        if (seed > 0)
            params["seed"] = seed;

        switch (constraint_order) {
        case 1:
            params["constraint-order"] = std::string("reversing");
            break;
        case 2:
            params["constraint-order"] = std::string("random-sorting");
            break;
        case 3:
            params["constraint-order"] = std::string("infeasibility-decr");
            break;
        case 4:
            params["constraint-order"] = std::string("infeasibility-incr");
            break;
        case 0:
        default:
            params["constraint-order"] = std::string("none");
            break;
        }

        params["time-limit"] = time_limit;
        params["pushing-k-factor"] = pushing_k_factor;
        params["pushes-limit"] = pushes_limit;
        params["pushing-objective-amplifier"] = pushing_objective_amplifier;
        params["pushing-iteration-limit"] = pushing_iteration_limit;
        params["thread"] = thread;

        auto result = lp::solve(pb, params);

        if (result.status == lp::result_status::success)
            return List::create(result.remaining_constraints, result.value);

        return List::create(result.remaining_constraints, NA_REAL);
    } catch (const std::bad_alloc& e) {
        Rcpp::Rcout << "lp memory error: " << e.what() << '\n';
    } catch (const std::exception& e) {
        Rcpp::Rcout << "lp error: " << e.what() << '\n';
    } catch (...) {
        Rcpp::Rcout << "lp error: unknown error\n";
    }

    return List::create(NA_INTEGER, NA_REAL);
}
