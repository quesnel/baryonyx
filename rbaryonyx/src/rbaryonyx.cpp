
/* Copyright (C) 2016-2018 INRA
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

#include <baryonyx/core>

using namespace Rcpp;

static void
r_write(int level, std::string msg)
{
    switch (level) {
    case 0:
        Rprintf("baryonyx: System is unusable: ");
        break;
    case 1:
        Rprintf("baryonyx: Action must be taken immediately: ");
        break;
    case 2:
        Rprintf("baryonyx: critical conditions: ");
        break;
    case 3:
        Rprintf("baryonyx: error conditions: ");
        break;
    case 4:
        Rprintf("baryonyx: warning conditions: ");
        break;
    default:
        break;
    }

    Rprintf("%s", msg.c_str());
};

static void
assign_parameters(const baryonyx::context_ptr& ctx,
                  int limit,
                  double theta,
                  double delta,
                  int constraint_order,
                  double kappa_min,
                  double kappa_step,
                  double kappa_max,
                  double alpha,
                  int w,
                  double time_limit,
                  int seed,
                  int thread,
                  int norm,
                  double pushing_k_factor,
                  double pushing_objective_amplifier,
                  int pushes_limit,
                  int pushing_iteration_limit,
                  int init_policy,
                  double init_random,
                  int float_type)
{
    baryonyx::context_set_parameter(ctx, "limit", limit);
    baryonyx::context_set_parameter(ctx, "theta", theta);
    baryonyx::context_set_parameter(ctx, "delta", delta);

    std::string order("none");
    switch (constraint_order) {
    case 1:
        order = "reversing";
        break;
    case 2:
        order = "random-sorting";
        break;
    case 3:
        order = "infeasibility-decr";
        break;
    case 4:
        order = "infeasibility-incr";
        break;
    default:
        break;
    }
    baryonyx::context_set_parameter(ctx, "constraint-order", order);

    baryonyx::context_set_parameter(ctx, "kappa-min", kappa_min);
    baryonyx::context_set_parameter(ctx, "kappa-step", kappa_step);
    baryonyx::context_set_parameter(ctx, "kappa-max", kappa_max);
    baryonyx::context_set_parameter(ctx, "alpha", alpha);
    baryonyx::context_set_parameter(ctx, "w", w);
    baryonyx::context_set_parameter(ctx, "time-limit", time_limit);

    if (seed > 0)
        baryonyx::context_set_parameter(ctx, "seed", seed);

    baryonyx::context_set_parameter(ctx, "thread", thread);

    baryonyx::context_set_parameter(
      ctx,
      "norm",
      (norm == 0)
        ? "none"
        : (norm == 1) ? "rng"
                      : (norm == 2) ? "l1" : (norm == 3) ? "l2" : "inf");

    baryonyx::context_set_parameter(ctx, "pushing-k-factor", pushing_k_factor);
    baryonyx::context_set_parameter(
      ctx, "pushing-objective-amplifier", pushing_objective_amplifier);
    baryonyx::context_set_parameter(ctx, "pushes-limit", pushes_limit);
    baryonyx::context_set_parameter(
      ctx, "pushing-iteration-limit", pushing_iteration_limit);

    std::string policy("bastert");
    switch (init_policy) {
    case 1:
        policy = "random";
        break;
    case 2:
        policy = "best";
        break;
    default:
        break;
    }
    baryonyx::context_set_parameter(ctx, "init-policy", policy);
    baryonyx::context_set_parameter(ctx, "init-random", init_random);

    baryonyx::context_set_parameter(
      ctx,
      "floating-point-type",
      (float_type == 1) ? "double"
                        : (float_type == 2) ? "longdouble" : "float");
}

// static double
// linearize_result(const baryonyx::result& r,
//                  double obj_min,
//                  double obj_max,
//                  baryonyx::objective_function_type type)
// {
//     if (r.status == baryonyx::result_status::success) {
//         if (type == baryonyx::objective_function_type::maximize)
//             return -r.value;
//         else
//             return r.value;
//     }

//     if (type == baryonyx::objective_function_type::maximize)
//         return -obj_min + r.remaining_constraints;

//     return obj_max + r.remaining_constraints;
// }

static List
convert_result(const baryonyx::result& res,
               bool minimize,
               const std::vector<double>& objective)
{
    bool solution;
    bool error;
    double value = 0;
    double duration = res.duration;
    int variables = res.variables;
    int constraints = res.constraints;
    int remaining = res.remaining_constraints;

    switch (res.status) {
    case baryonyx::result_status::success:
        Rprintf("Solution found: %f\n", res.solutions.back().value);
        value = res.solutions.back().value;
        solution = true;
        error = false;
        break;
    case baryonyx::result_status::internal_error:
        Rprintf("Internal error\n");
        solution = false;
        error = true;
        break;
    default:
        Rprintf("No solution found. Constraints remaining: %d\n",
                res.remaining_constraints);
        solution = false;
        error = false;
        break;
    }

    return List::create(_["solution_found"] = solution,
                        _["error_found"] = error,
                        _["value"] = value,
                        _["duration"] = duration,
                        _["variables"] = variables,
                        _["constraints"] = constraints,
                        _["remaining_constraints"] = remaining,
                        _["minimize"] = minimize,
                        _["objective_function"] = objective);
}

//' Tries to solve the 01 linear programming problem.
//'
//' @param constraint_order between each run over R. DEfault is to use the
//'    order defined in model (lp file).
//'    - 0: none
//'    - 1: reversing
//'    - 2: random-sorting
//'    - 3: infeasibility-decr
//'    - 4: infeasibility-incr
//'
//' @param norm the normalization used to reduces cost dependencies. Default
//'    is to use the infinity norm. Other values are:
//'    - 0: none (use the default cost from the model (lp file)).
//'    - 1: rng
//'    - 2: l1-norm
//'    - 3: l2-norm
//'    - 4: infinity-norm
//'
//' @param policy_type the type of (re)initilization used into the solver.
//'     Default is to use bastert and 0.5 for policy_random.
//'    - 0: bastert
//'    - 1: random
//'    - 2: best
//'
//' @param float_type the type of real used into the solver. Default is to
//'     use the C/C++ double representation.
//'    - 0: float
//'    - 1: double
//'    - 2: longdouble
//'
//' @return a named list:
//'   \item{solution_found}{Boolean, TRUE is a solution is found.}
//'   \item{error_found}{Boolean, TRUE if an error occurred in Baryonyx,
//'     other parameter may be useless and undefined.}
//'   \item{value}{Real, the value of the solution if solution_found is
//'     TRUE.}
//'   \item{duration}{Real, time in second taken to found or not found a
//'     solution.}
//'   \item{variables}{Integer, number of variables in the problem.}
//'   \item{constraints}{Integer, number of constraints in the problem.}
//'   \item{remaining_constraints}{Integer, constraints remaining, 0 if
//'     solution_found is TRUE, between 1 and constraints if
//'      solution_found is FALSE.}
//'   \item{minimize}{Boolean, TRUE is problem is a minimization, FALSE if
//'     the problem is a maximization.}
//'   \item{objective_function}{Vector, bound of the objective function.}
//'
//' @useDynLib rbaryonyx
//' @importFrom Rcpp sourceCpp
//'
//' @export
// [[Rcpp::export]]
List
solve_01lp_problem(std::string file_path,
                   int limit = 1000,
                   double theta = 0.5,
                   double delta = -1,
                   int constraint_order = 0,
                   double kappa_min = 0.0,
                   double kappa_step = 1.0e-3,
                   double kappa_max = 0.6,
                   double alpha = 1.0,
                   int w = 20,
                   double time_limit = 10.0,
                   int seed = -1,
                   int thread = 1,
                   int norm = 4,
                   double pushing_k_factor = 0.9,
                   double pushing_objective_amplifier = 5.0,
                   int pushes_limit = 100,
                   int pushing_iteration_limit = 50,
                   int policy_type = 0,
                   double policy_random = 0.5,
                   int float_type = 1,
                   bool verbose = true) noexcept

{
    baryonyx::result ret;
    std::vector<double> objective(2);
    bool minimize = true;

    try {
        auto ctx = baryonyx::make_context(&r_write, verbose ? 6 : 4);
        auto pb = baryonyx::make_problem(ctx, file_path);

        std::tie(objective[0], objective[1]) =
          baryonyx::compute_min_max_objective_function(pb);

        minimize = (pb.type == baryonyx::objective_function_type::minimize);

        assign_parameters(ctx,
                          limit,
                          theta,
                          delta,
                          constraint_order,
                          kappa_min,
                          kappa_step,
                          kappa_max,
                          alpha,
                          w,
                          time_limit,
                          seed,
                          thread,
                          norm,
                          pushing_k_factor,
                          pushing_objective_amplifier,
                          pushes_limit,
                          pushing_iteration_limit,
                          policy_type,
                          policy_random,
                          float_type);

        ret = baryonyx::solve(ctx, pb);
    } catch (const std::bad_alloc& e) {
        Rprintf("Baryonyx memory error: %s\n", e.what());
        ret.status = baryonyx::result_status::internal_error;
    } catch (const std::exception& e) {
        Rprintf("Baryonyx error: %s\n", e.what());
        ret.status = baryonyx::result_status::internal_error;
    } catch (...) {
        Rprintf("Baryonyx error: unknown error\n");
        ret.status = baryonyx::result_status::internal_error;
    }

    return convert_result(ret, minimize, objective);
}

//' Tries to optimize the 01 linear programming problem.
//'
//' @param constraint_order between each run over R. DEfault is to use the
//'    order defined in model (lp file).
//'    - 0: none
//'    - 1: reversing
//'    - 2: random-sorting
//'    - 3: infeasibility-decr
//'    - 4: infeasibility-incr
//'
//' @param norm the normalization used to reduces cost dependencies. Default
//'    is to use the infinity norm. Other values are:
//'    - 0: none (use the default cost from the model (lp file)).
//'    - 1: rng
//'    - 2: l1-norm
//'    - 3: l2-norm
//'    - 4: infinity-norm
//'
//' @param policy_type the type of (re)initilization used into the solver.
//'     Default is to use bastert and 0.5 for policy_random.
//'    - 0: bastert
//'    - 1: random
//'    - 2: best
//'
//' @param float_type the type of real used into the solver. Default is to
//'     use the C/C++ double representation.
//'    - 0: float
//'    - 1: double
//'    - 2: longdouble
//'
//' @return a named list:
//'   \item{solution_found}{Boolean, TRUE is a solution is found.}
//'   \item{error_found}{Boolean, TRUE if an error occurred in Baryonyx,
//'     other parameter may be useless and undefined.}
//'   \item{value}{Real, the value of the solution if solution_found is
//'     TRUE.}
//'   \item{duration}{Real, time in second taken to found or not found a
//'     solution.}
//'   \item{variables}{Integer, number of variables in the problem.}
//'   \item{constraints}{Integer, number of constraints in the problem.}
//'   \item{remaining_constraints}{Integer, constraints remaining, 0 if
//'     solution_found is TRUE, between 1 and constraints if
//'      solution_found is FALSE.}
//'   \item{minimize}{Boolean, TRUE is problem is a minimization, FALSE if
//'     the problem is a maximization.}
//'   \item{objective_function}{Vector, bound of the objective function.}
//'
//' @useDynLib rbaryonyx
//' @importFrom Rcpp sourceCpp
//'
//' @export
// [[Rcpp::export]]
List
optimize_01lp_problem(std::string file_path,
                      int limit = 1000,
                      double theta = 0.5,
                      double delta = -1,
                      int constraint_order = 0,
                      double kappa_min = 0.0,
                      double kappa_step = 1.0e-3,
                      double kappa_max = 0.6,
                      double alpha = 1.0,
                      int w = 20,
                      double time_limit = 10.0,
                      int seed = -1,
                      int thread = 1,
                      int norm = 4,
                      double pushing_k_factor = 0.9,
                      double pushing_objective_amplifier = 5.0,
                      int pushes_limit = 100,
                      int pushing_iteration_limit = 50,
                      int policy_type = 0,
                      double policy_random = 0.5,
                      int float_type = 1,
                      bool verbose = true) noexcept
{
    baryonyx::result ret;
    std::vector<double> objective(2);
    bool minimize = true;

    try {
        auto ctx = baryonyx::make_context(&r_write, verbose ? 6 : 4);
        auto pb = baryonyx::make_problem(ctx, file_path);

        std::tie(objective[0], objective[1]) =
          baryonyx::compute_min_max_objective_function(pb);

        minimize = (pb.type == baryonyx::objective_function_type::minimize);

        assign_parameters(ctx,
                          limit,
                          theta,
                          delta,
                          constraint_order,
                          kappa_min,
                          kappa_step,
                          kappa_max,
                          alpha,
                          w,
                          time_limit,
                          seed,
                          thread,
                          norm,
                          pushing_k_factor,
                          pushing_objective_amplifier,
                          pushes_limit,
                          pushing_iteration_limit,
                          policy_type,
                          policy_random,
                          float_type);

        ret = baryonyx::optimize(ctx, pb);
    } catch (const std::bad_alloc& e) {
        Rprintf("Baryonyx memory error: %s\n", e.what());
        ret.status = baryonyx::result_status::internal_error;
    } catch (const std::exception& e) {
        Rprintf("Baryonyx error: %s\n", e.what());
        ret.status = baryonyx::result_status::internal_error;
    } catch (...) {
        Rprintf("Baryonyx error: unknown error\n");
        ret.status = baryonyx::result_status::internal_error;
    }

    return convert_result(ret, minimize, objective);
}
