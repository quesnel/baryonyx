/* Copyright (C) 2017 INRA
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

#define R_USE_C99_IN_CXX
#include <Rcpp.h>

#include <lpcore>

using namespace Rcpp;

class Rcontext : public lp::context::logger
{
public:
    Rcontext() = default;
    ~Rcontext() noexcept override = default;

    void write(int priority,
               const char* file,
               int line,
               const char* fn,
               const char* format,
               va_list args) noexcept override
    {
        Rprintf("lp: %d at %d in function %s form file %s: ",
                priority,
                line,
                fn,
                file);

        Rvprintf(format, args);
    }

    void write(lp::context::message_type m,
               const char* format,
               va_list args) noexcept override
    {
        switch (m) {
        case lp::context::message_type::emerg:
            Rprintf("lp: system is unusable\n");
            break;
        case lp::context::message_type::alert:
            Rprintf("lp: action must be taken immediately\n");
            break;
        case lp::context::message_type::crit:
            Rprintf("lp: critical conditions\n");
            break;
        case lp::context::message_type::err:
            Rprintf("lp: error conditions\n");
            break;
        case lp::context::message_type::warning:
            Rprintf("lp: warning conditions\n");
            break;
        case lp::context::message_type::notice:
        case lp::context::message_type::info:
        case lp::context::message_type::debug:
            break;
        }

        Rvprintf(format, args);
    }
};

//' Tries to solve the 01 linear programming problem.
//'
//' @param constraint_order: 0-none, 1-reversing, 2-random-sorting,
//' 3-infeasibility-decr, 4- infeasibility-incr
//'
//' @return A list of three values: The remaining constraint and the value
//' of the solution (Integer, Real, Real):
//'   - Integer: The number of remaining constraints. 0 means a solution
//'      found, greater than 0 means remaining constraints and value is
//'      0. -1 means an error occurred during solving. Use the @c `lp`
//'      program to check out the error.
//'   - Real: The value of the solution found (if the remaining
//'     constraint equal 0).
//'   - Real: The value [objective_function_min, objective_function_max +
//'     constraints number].
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
                   long int thread = 1,
                   bool verbose = true) noexcept
{
    try {
        auto ctx = std::make_shared<lp::context>();
        ctx->set_logger(std::make_unique<Rcontext>());
        ctx->set_log_priority(!verbose ? 6 : 3);

        auto pb = lp::make_problem(ctx, file_path);
        auto mm = lp::compute_min_max_objective_function(pb);

        if (pb.type == lp::objective_function_type::maximize) {
            std::swap(std::get<0>(mm), std::get<1>(mm));
            std::get<0>(mm) = -std::get<0>(mm);
            std::get<1>(mm) = -std::get<1>(mm);
        }

        std::map<std::string, lp::parameter> params;
        params["limit"] = limit;
        params["theta"] = theta;
        params["delta"] = delta;

        params["kappa-min"] = kappa_min;
        params["kappa-step"] = kappa_step;
        params["kappa-max"] = kappa_max;
        params["alpha"] = alpha;
        params["w"] = w;

        if (seed > 0) {
            Rprintf("solver uses a PRNG with %ld as seed\n", seed);
            params["seed"] = seed;
        } else {
            Rprintf("solver uses a PRNG with a random seed.\n");
        }

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

        auto result = lp::solve(ctx, pb, params);
        if (result.status == lp::result_status::success) {
            if (pb.type == lp::objective_function_type::maximize)
                return List::create(
                  result.remaining_constraints, result.value, -result.value);
            else
                return List::create(
                  result.remaining_constraints, result.value, result.value);
        }

        return List::create(result.remaining_constraints,
                            result.value,
                            std::get<1>(mm) + result.remaining_constraints);

    } catch (const std::bad_alloc& e) {
        Rprintf("lp memory error: %s\n", e.what());
    } catch (const std::exception& e) {
        Rprintf("lp error: %s\n", e.what());
    } catch (...) {
        Rprintf("lp error: unknown error\n");
    }

    return List::create(
      NA_INTEGER, NA_REAL, std::numeric_limits<double>::infinity());
}

//' Tries to optimize the 01 linear programming problem.
//'
//' @param constraint_order: 0-none, 1-reversing, 2-random-sorting,
//' 3-infeasibility-decr, 4- infeasibility-incr
//'
//' @return an uniq scalar.
//'   - If a solution is found:
//'     - if the problem is a minimization: the value of the solution found.
//'     - if the problem is a maximization: the inverse of the solution
//'       found.
//'   - If no solution is found, we use the limits of the objective
//'     function (minimal and maximal value possible.
//'     - if the problem is a minimization: the maximal value possible +
//'       the remaining constraints.
//'     - if the problem is a maximization: the inverse of the minimal
//'       value possible + the remaining constraints.
//'   - If a error occurred (not enough memory, problem error etc.):
//'     - if the problem is a minimization: the maximal value possible +
//'       the number of constraints .
//'     - if the problem is a maximization: the inverse of the minimal
//'       value possible + the number of constraints.
//'
//' @useDynLib rilp
//' @importFrom Rcpp sourceCpp
//'
//' @export
// [[Rcpp::export]]
double
optimize_01lp_problem(std::string file_path,
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
                      long int thread = 1,
                      bool verbose = true) noexcept
{
    double ret{ 0 };

    try {
        auto ctx = std::make_shared<lp::context>();
        ctx->set_logger(std::make_unique<Rcontext>());
        ctx->set_log_priority(verbose ? 6 : 3);

        auto pb = lp::make_problem(ctx, file_path);
        auto mm = lp::compute_min_max_objective_function(pb);

        if (pb.type == lp::objective_function_type::maximize)
            ret = -std::get<0>(mm) + size(pb);
        else
            ret = std::get<1>(mm) + size(pb);

        if (pb.type == lp::objective_function_type::maximize) {
            std::swap(std::get<0>(mm), std::get<1>(mm));
            std::get<0>(mm) = -std::get<0>(mm);
            std::get<1>(mm) = -std::get<1>(mm);
        }

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

        auto result = lp::optimize(ctx, pb, params);

        if (result.status == lp::result_status::success) {
            if (pb.type == lp::objective_function_type::maximize)
                ret = -result.value;
            else
                ret = result.value;
        } else {
            if (pb.type == lp::objective_function_type::maximize)
                ret = -std::get<0>(mm) + result.remaining_constraints;
            else
                ret = std::get<1>(mm) + result.remaining_constraints;
        }
    } catch (const std::bad_alloc& e) {
        Rprintf("lp memory error: %s\n", e.what());
    } catch (const std::exception& e) {
        Rprintf("lp error: %s\n", e.what());
    } catch (...) {
        Rprintf("lp error: unknown error\n");
    }

    Rprintf("optimizer returns: %f (kappa: %f %f %f delta: %f theta: %f)\n",
            ret,
            kappa_min,
            kappa_step,
            kappa_max,
            delta,
            theta);

    return ret;
}
