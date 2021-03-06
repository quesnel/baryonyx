/* Copyright (C) 2016-2019 INRA
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
r_solver_started_cb(const baryonyx::solver_parameters& /*params*/)
{
    R_FlushConsole();
    R_ProcessEvents();
}

static void
r_solver_updated_cb(int remaining_constraints,
                    double value,
                    int loop,
                    double duration,
                    long int reinit_number)
{
    if (remaining_constraints > 0) {
        Rprintf("  - Constraints remaining: %d (loop: %d t: %fs reinit:%ld)\n",
                remaining_constraints,
                loop,
                duration,
                reinit_number);
    } else {
        if (loop >= 0)
            Rprintf("  - Solution found: %f (loop: %d t: %fs reinit:%ld)\n",
                    value,
                    loop,
                    duration,
                    reinit_number);
        else
            Rprintf(
              "  - Solution found via push: %f (loop: %d t: %fs reinit:%ld)\n",
              value,
              -loop,
              duration,
              reinit_number);
    }

    R_FlushConsole();
    R_ProcessEvents();
}

static void
r_solver_finished_cb(const baryonyx::result& r)
{
    switch (r.status) {
    case baryonyx::result_status::success:
        if (r.loop >= 0)
            Rprintf("Best solution found: %f in %d loop and %fs\n",
                    r.solutions.back().value,
                    r.loop,
                    r.duration);
        else
            Rprintf("Best solution found via push: %f in %d loop and %fs\n",
                    r.solutions.back().value,
                    -r.loop,
                    r.duration);
        break;
    case baryonyx::result_status::internal_error:
        Rprintf("No solution. Internal error\n");
        break;
    case baryonyx::result_status::uninitialized:
        Rprintf("No solution. Uninitialized error\n");
        break;
    case baryonyx::result_status::kappa_max_reached:
        Rprintf(
          "No solution. Constraint remaining: %d. Kappa reached in %fs.\n",
          r.remaining_constraints,
          r.duration);
        break;
    case baryonyx::result_status::time_limit_reached:
        Rprintf("No solution. Constraint remaining: %d. Time limit reached "
                "at %fs.\n",
                r.remaining_constraints,
                r.duration);
        break;
    case baryonyx::result_status::limit_reached:
        Rprintf("No solution. Constraint remaining: %d. Loop limit reached "
                "in %fs.\n",
                r.remaining_constraints,
                r.duration);
        break;
    }

    R_FlushConsole();
    R_ProcessEvents();
}

static baryonyx::solver_parameters::pre_constraint_order
get_pre_constraint_order(int prep)
{
    namespace bx = baryonyx;

    switch (prep) {
    case 0:
        return bx::solver_parameters::pre_constraint_order::none;
    case 2:
        return bx::solver_parameters::pre_constraint_order::less_greater_equal;
    case 3:
        return bx::solver_parameters::pre_constraint_order::less_equal_greater;
    case 4:
        return bx::solver_parameters::pre_constraint_order::greater_less_equal;
    case 5:
        return bx::solver_parameters::pre_constraint_order::greater_equal_less;
    case 6:
        return bx::solver_parameters::pre_constraint_order::equal_less_greater;
    case 7:
        return bx::solver_parameters::pre_constraint_order::equal_greater_less;
    case 8:
        return bx::solver_parameters::pre_constraint_order::p1;
    case 9:
        return bx::solver_parameters::pre_constraint_order::p2;
    case 10:
        return bx::solver_parameters::pre_constraint_order::p3;
    case 11:
        return bx::solver_parameters::pre_constraint_order::p4;
    default:
        return bx::solver_parameters::pre_constraint_order::memory;
    }
}

static baryonyx::solver_parameters::constraint_order
get_constrait_order(int order)
{
    namespace bx = baryonyx;

    switch (order) {
    case 1:
        return bx::solver_parameters::constraint_order::reversing;
    case 2:
        return bx::solver_parameters::constraint_order::random_sorting;
    case 3:
        return bx::solver_parameters::constraint_order::infeasibility_decr;
    case 4:
        return bx::solver_parameters::constraint_order::infeasibility_incr;
    case 5:
        return bx::solver_parameters::constraint_order::lagrangian_decr;
    case 6:
        return bx::solver_parameters::constraint_order::lagrangian_incr;
    case 7:
        return bx::solver_parameters::constraint_order::pi_sign_change;
    case 8:
        return bx::solver_parameters::constraint_order::cycle;
    default:
        return bx::solver_parameters::constraint_order::none;
    }
}

static baryonyx::solver_parameters::cost_norm_type
get_cost_norm(int norm)
{
    namespace bx = baryonyx;

    switch (norm) {
    case 0:
        return bx::solver_parameters::cost_norm_type::none;
    case 1:
        return bx::solver_parameters::cost_norm_type::random;
    case 2:
        return bx::solver_parameters::cost_norm_type::l1;
    case 3:
        return bx::solver_parameters::cost_norm_type::l2;
    default:
        return bx::solver_parameters::cost_norm_type::loo;
    }
}

static baryonyx::solver_parameters::init_policy_type
get_init_policy(int type)
{
    namespace bx = baryonyx;

    switch (type) {
    case 0:
        return bx::solver_parameters::init_policy_type::bastert;
    case 1:
        return bx::solver_parameters::init_policy_type::pessimistic_solve;
    case 2:
        return bx::solver_parameters::init_policy_type::optimistic_solve;
    default:
        return bx::solver_parameters::init_policy_type::bastert;
    }
}

static baryonyx::solver_parameters::floating_point_type
get_floating_point_type(int type)
{
    namespace bx = baryonyx;

    switch (type) {
    case 0:
        return bx::solver_parameters::floating_point_type::float_type;
    case 2:
        return bx::solver_parameters::floating_point_type::longdouble_type;
    default:
        return bx::solver_parameters::floating_point_type::double_type;
    }
}

static baryonyx::solver_parameters::storage_type
get_storage_type(int type)
{
    namespace bx = baryonyx;

    switch (type) {
    case 1:
        return bx::solver_parameters::storage_type::bound;
    case 2:
        return bx::solver_parameters::storage_type::five;
    default:
        return bx::solver_parameters::storage_type::one;
    }
}

static List
convert_result(const baryonyx::result& res, bool minimize)
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

    std::vector<double> solutions(res.solutions.size());
    std::transform(
      res.solutions.cbegin(),
      res.solutions.cend(),
      solutions.begin(),
      [](const baryonyx::solution& elem) -> double { return elem.value; });

    return List::create(_["solution_found"] = solution,
                        _["error_found"] = error,
                        _["value"] = value,
                        _["duration"] = duration,
                        _["variables"] = variables,
                        _["constraints"] = constraints,
                        _["remaining_constraints"] = remaining,
                        _["minimize"] = minimize,
                        _["solutions"] = solutions);
}

//' Tries to solve the 01 linear programming problem.
//'
//' @param pre_constraint_order order of the raw problem.
//'     0: none
//'     1: memory
//'     2: less_greater_equal
//'     3: less_equal_greater
//'     4: greater_less_equal
//'     5: greater_equal_less
//'     6: equal_less_greater
//'     7: equal_greater_less
//'     8: p1
//'     9: p2
//'     10: p3
//'     11: p4
//'
//' @param constraint_order order between each run over R. DEfault is to use
//'    the order defined in model (lp file).
//'    - 0: none
//'    - 1: reversing
//'    - 2: random-sorting
//'    - 3: infeasibility-decr
//'    - 4: infeasibility-incr
//'    - 5: lagrangian_decr
//'    - 6: lagrangian_incr
//'    - 7: pi_sign_change
//'    - 8: cycle
//'
//' @param norm the normalization used to reduces cost dependencies. Default
//'    is to use the infinity norm. Other values are:
//'    - 0: none (use the default cost from the model (lp file)).
//'    - 1: rng
//'    - 2: l1-norm
//'    - 3: l2-norm
//'    - 4: infinity-norm
//'
//' @param init_policy the type of initilization used into the solver.
//'    - 0: bastert,
//'    - 1: pessimistic_solve,
//'    - 2: optimistic_solve,
//'
//' @param init_policy_random the type percentage of random in the init_policy.
//'
//' @param float_type the type of real used into the solver. Default is to
//'     use the C/C++ double representation.
//'    - 0: float
//'    - 1: double
//'    - 2: longdouble
//'
//' @param storage_type the type of solution storage. Default is to store
//'     the best solution.
//'    - 0: stores only the best solution found.
//'    - 1: stores the best and the bad solution found.
//'    - 2: stores the five best solutions found.
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
//'   \item{solutions}{Vector, all solution found (may be empty).}
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
                   int pre_constraint_order = 1,
                   int constraint_order = 0,
                   double kappa_min = 0.0,
                   double kappa_step = 1.0e-3,
                   double kappa_max = 0.6,
                   double alpha = 1.0,
                   double w = 0.05,
                   double time_limit = 10.0,
                   int seed = -1,
                   int thread = 1,
                   int norm = 4,
                   double pushing_k_factor = 0.9,
                   double pushing_objective_amplifier = 5.0,
                   int pushes_limit = 100,
                   int pushing_iteration_limit = 50,
                   int init_policy = 0,
                   double init_policy_random = 0.5,
                   int float_type = 1,
                   int storage_type = 1,
                   bool verbose = true) noexcept

{
    baryonyx::result ret;
    bool minimize = true;

    try {
        auto ctx = baryonyx::make_context(verbose ? 6 : 4);
        baryonyx::context_register(
          ctx, r_solver_started_cb, r_solver_updated_cb, r_solver_finished_cb);

        auto pb = baryonyx::make_problem(ctx, file_path);

        minimize = (pb.type == baryonyx::objective_function_type::minimize);

        baryonyx::solver_parameters params;

        params.time_limit = time_limit;
        params.theta = theta;
        params.delta = delta;
        params.kappa_min = kappa_min;
        params.kappa_step = kappa_step;
        params.kappa_max = kappa_max;
        params.alpha = alpha;
        params.pushing_k_factor = pushing_k_factor;
        params.pushing_objective_amplifier = pushing_objective_amplifier;
        params.init_policy_random = init_policy_random;
        params.seed = seed;
        params.thread = thread;
        params.limit = limit;
        params.w = w;
        params.pushes_limit = pushes_limit;
        params.pushing_iteration_limit = pushing_iteration_limit;
        params.pre_order = get_pre_constraint_order(pre_constraint_order);
        params.order = get_constrait_order(constraint_order);
        params.cost_norm = get_cost_norm(norm);
        params.init_policy = get_init_policy(init_policy);
        params.float_type = get_floating_point_type(float_type);
        params.storage = get_storage_type(storage_type);

        baryonyx::context_set_solver_parameters(ctx, params);

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

    return convert_result(ret, minimize);
}

//' Tries to optimize the 01 linear programming problem.
//'
//' @param pre_constraint_order order of the raw problem.
//'     0: none
//'     1: memory
//'     2: less_greater_equal
//'     3: less_equal_greater
//'     4: greater_less_equal
//'     5: greater_equal_less
//'     6: equal_less_greater
//'     7: equal_greater_less
//'     8: p1
//'     9: p2
//'     10: p3
//'     11: p4
//'
//' @param constraint_order between each run over R. DEfault is to use the
//'    order defined in model (lp file).
//'    - 0: none
//'    - 1: reversing
//'    - 2: random-sorting
//'    - 3: infeasibility-decr
//'    - 4: infeasibility-incr
//'    - 5: lagrangian_decr
//'    - 6: lagrangian_incr
//'    - 7: pi_sign_change
//'    - 8: cycle
//'
//' @param norm the normalization used to reduces cost dependencies. Default
//'    is to use the infinity norm. Other values are:
//'    - 0: none (use the default cost from the model (lp file)).
//'    - 1: rng
//'    - 2: l1-norm
//'    - 3: l2-norm
//'    - 4: infinity-norm
//'
//' @param float_type the type of real used into the solver. Default is to
//'     use the C/C++ double representation.
//'    - 0: float
//'    - 1: double
//'    - 2: longdouble
//'
//' @param storage_type the type of solution storage. Default is to store
//'     the best solution.
//'    - 0: stores only the best solution found.
//'    - 1: stores the best and the bad solution found.
//'    - 2: stores the five best solutions found.
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
//'   \item{solutions}{Vector, all solution found (may be empty).}
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
                      int pre_constraint_order = 1,
                      int constraint_order = 0,
                      double kappa_min = 0.0,
                      double kappa_step = 1.0e-3,
                      double kappa_max = 0.6,
                      double alpha = 1.0,
                      double w = 0.05,
                      double time_limit = 10.0,
                      int seed = -1,
                      int thread = 1,
                      int norm = 4,
                      double pushing_k_factor = 0.9,
                      double pushing_objective_amplifier = 5.0,
                      int pushes_limit = 100,
                      int pushing_iteration_limit = 50,
                      double init_crossover_bastert_insertion = 0.01,
                      double init_crossover_solution_selection_mean = 0.0,
                      double init_crossover_solution_selection_stddev = 0.3,
                      double init_mutation_variable_mean = 0.0001,
                      double init_mutation_variable_stddev = 0.001,
                      double init_mutation_value_mean = 0.5,
                      double init_mutation_value_stddev = 0.2,
                      double init_kappa_improve_start = 0,
                      double init_kappa_improve_increase = 0.02,
                      double init_kappa_improve_stop = 0.2,
                      int float_type = 1,
                      int storage_type = 1,
                      bool verbose = true) noexcept
{
    baryonyx::result ret;
    bool minimize = true;

    try {
        auto ctx = baryonyx::make_context(verbose ? 6 : 4);
        baryonyx::context_register(
          ctx, r_solver_started_cb, r_solver_updated_cb, r_solver_finished_cb);

        auto pb = baryonyx::make_problem(ctx, file_path);

        minimize = (pb.type == baryonyx::objective_function_type::minimize);

        baryonyx::solver_parameters params;

        params.time_limit = time_limit;
        params.theta = theta;
        params.delta = delta;
        params.kappa_min = kappa_min;
        params.kappa_step = kappa_step;
        params.kappa_max = kappa_max;
        params.alpha = alpha;
        params.pushing_k_factor = pushing_k_factor;
        params.pushing_objective_amplifier = pushing_objective_amplifier;

        params.init_crossover_bastert_insertion =
          init_crossover_bastert_insertion;
        params.init_crossover_solution_selection_mean =
          init_crossover_solution_selection_mean;
        params.init_crossover_solution_selection_stddev =
          init_crossover_solution_selection_stddev;
        params.init_mutation_variable_mean = init_mutation_variable_mean;
        params.init_mutation_variable_stddev = init_mutation_variable_stddev;
        params.init_mutation_value_mean = init_mutation_value_mean;
        params.init_mutation_value_stddev = init_mutation_value_stddev;
        params.init_kappa_improve_start = init_kappa_improve_start;
        params.init_kappa_improve_increase = init_kappa_improve_increase;
        params.init_kappa_improve_stop = init_kappa_improve_stop;

        params.seed = seed;
        params.thread = thread;
        params.limit = limit;
        params.w = w;
        params.pushes_limit = pushes_limit;
        params.pushing_iteration_limit = pushing_iteration_limit;
        params.pre_order = get_pre_constraint_order(pre_constraint_order);
        params.order = get_constrait_order(constraint_order);
        params.cost_norm = get_cost_norm(norm);
        params.float_type = get_floating_point_type(float_type);
        params.storage = get_storage_type(storage_type);

        baryonyx::context_set_solver_parameters(ctx, params);

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

    return convert_result(ret, minimize);
}
