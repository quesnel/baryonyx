
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

#include <baryonyx/core>

using namespace Rcpp;

namespace {

class Rcontext : public baryonyx::context::logger
{
public:
    Rcontext() = default;
    ~Rcontext() noexcept override = default;

    void write(int priority,
               const char* file,
               int line,
               const char* fn,
               const char* format,
               va_list args) noexcept override final
    {
        Rprintf("baryonyx: %d at %d in function %s form file %s: ",
                priority,
                line,
                fn,
                file);

        Rvprintf(format, args);
    }

    void write(baryonyx::context::message_type m,
               const char* format,
               va_list args) noexcept override final
    {
        switch (m) {
        case baryonyx::context::message_type::emerg:
            Rprintf("baryonyx: System is unusable: ");
            break;
        case baryonyx::context::message_type::alert:
            Rprintf("baryonyx: Action must be taken immediately: ");
            break;
        case baryonyx::context::message_type::crit:
            Rprintf("baryonyx: critical conditions: ");
            break;
        case baryonyx::context::message_type::err:
            Rprintf("baryonyx: error conditions: ");
            break;
        case baryonyx::context::message_type::warning:
            Rprintf("baryonyx: warning conditions: ");
            break;
        case baryonyx::context::message_type::notice:
        case baryonyx::context::message_type::info:
        case baryonyx::context::message_type::debug:
            break;
        }

        Rvprintf(format, args);
    }
};

} // anonymous namespace

static void
assign_parameters(std::shared_ptr<baryonyx::context> ctx,
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
                  int float_type)
{
    ctx->set_parameter("limit", limit);
    ctx->set_parameter("theta", theta);
    ctx->set_parameter("delta", delta);

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
    ctx->set_parameter("constraint-order", order);

    ctx->set_parameter("kappa-min", kappa_min);
    ctx->set_parameter("kappa-step", kappa_step);
    ctx->set_parameter("kappa-max", kappa_max);
    ctx->set_parameter("alpha", alpha);
    ctx->set_parameter("w", w);
    ctx->set_parameter("time-limit", time_limit);

    if (seed > 0)
        ctx->set_parameter("seed", seed);

    ctx->set_parameter("thread", thread);

    ctx->set_parameter(
      "norm",
      (norm == 0)
        ? "none"
        : (norm == 1) ? "rng"
                      : (norm == 2) ? "l1" : (norm == 3) ? "l2" : "inf");

    ctx->set_parameter("pushing-k-factor", pushing_k_factor);
    ctx->set_parameter("pushing-objective-amplifier",
                       pushing_objective_amplifier);
    ctx->set_parameter("pushes-limit", pushes_limit);
    ctx->set_parameter("pushing-iteration-limit", pushing_iteration_limit);

    ctx->set_parameter("floating-point-type",
                       (float_type == 1)
                         ? "double"
                         : (float_type == 2) ? "longdouble" : "float");
}

static double
linearize_result(const baryonyx::result& r,
                 double obj_min,
                 double obj_max,
                 baryonyx::objective_function_type type)
{
    if (r.status == baryonyx::result_status::success) {
        if (type == baryonyx::objective_function_type::maximize)
            return -r.value;
        else
            return r.value;
    }

    if (type == baryonyx::objective_function_type::maximize)
        return -obj_min + r.remaining_constraints;

    return obj_max + r.remaining_constraints;
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
//' @param norm: the normalization used to reduces cost dependencies. Default
//'    is to use the infinity norm. Other values are:
//'    - 0: none (use the default cost from the model (lp file)).
//'    - 1: rng
//'    - 2: l1-norm
//'    - 3: l2-norm
//'    - 4: infinity-norm
//'
//' @param float_type: the type of real used into the solver. Default is to
//'     use the C/C++ double representation.
//'    - 0: float
//'    - 1: double
//'    - 2: longdouble
//'
//' @return an unique scalar.
//'   - If a solution is found:
//'     - if the problem is a minimization: the value of the solution
//'       found.
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
//' @useDynLib rbaryonyx
//' @importFrom Rcpp sourceCpp
//'
//' @export
// [[Rcpp::export]]
double
solve_01lp_problem(std::string file_path,
                   int limit = 1000,
                   double theta = 0.5,
                   double delta = 1e-4,
                   int constraint_order = 0,
                   double kappa_min = 0.1,
                   double kappa_step = 1e-4,
                   double kappa_max = 1.0,
                   double alpha = 1.0,
                   int w = 500,
                   double time_limit = 10.0,
                   int seed = -1,
                   int thread = 1,
                   int norm = 4,
                   double pushing_k_factor = 0.9,
                   double pushing_objective_amplifier = 5.0,
                   int pushes_limit = 10,
                   int pushing_iteration_limit = 20,
                   int float_type = 1,
                   bool verbose = true) noexcept
{
    double ret{ 0 };

    try {
        auto ctx = std::make_shared<baryonyx::context>();
        ctx->set_logger(std::make_unique<Rcontext>());
        ctx->set_log_priority(!verbose ? 6 : 3);

        auto pb = baryonyx::make_problem(ctx, file_path);
        auto mm = baryonyx::compute_min_max_objective_function(pb);
        auto type = pb.type;

        if (type == baryonyx::objective_function_type::maximize)
            ret = -std::get<0>(mm) + size(pb);
        else
            ret = std::get<1>(mm) + size(pb);

        if (type == baryonyx::objective_function_type::maximize) {
            std::swap(std::get<0>(mm), std::get<1>(mm));
            std::get<0>(mm) = -std::get<0>(mm);
            std::get<1>(mm) = -std::get<1>(mm);
        }

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
                          float_type);

        auto result = baryonyx::solve(ctx, pb);
        ret = linearize_result(result, std::get<0>(mm), std::get<1>(mm), type);
    } catch (const std::bad_alloc& e) {
        Rprintf("lp memory error: %s\n", e.what());
    } catch (const std::exception& e) {
        Rprintf("lp error: %s\n", e.what());
    } catch (...) {
        Rprintf("lp error: unknown error\n");
    }

    Rprintf("solver returns: %f\n", ret);

    return ret;
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
//' @param norm: the normalization used to reduces cost dependencies.
// Default
//'    is to use the infinity norm. Other values are:
//'    - 0: none (use the default cost from the model (lp file)).
//'    - 1: rng
//'    - 2: l1-norm
//'    - 3: l2-norm
//'    - 4: infinity-norm
//'
//' @param float_type: the type of real used into the solver. Default is to
//'     use the C/C++ double representation.
//'    - 0: float
//'    - 1: double
//'    - 2: longdouble
//'
//' @return an unique scalar.
//'   - If a solution is found:
//'     - if the problem is a minimization: the value of the solution
//'       found.
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
//' @useDynLib rbaryonyx
//' @importFrom Rcpp sourceCpp
//'
//' @export
// [[Rcpp::export]]
double
optimize_01lp_problem(std::string file_path,
                      int limit = 1000,
                      double theta = 0.5,
                      double delta = 1e-4,
                      int constraint_order = 0,
                      double kappa_min = 0.1,
                      double kappa_step = 1e-4,
                      double kappa_max = 1.0,
                      double alpha = 1.0,
                      int w = 500,
                      double time_limit = 10.0,
                      int seed = -1,
                      int thread = 1,
                      int norm = 4,
                      double pushing_k_factor = 0.9,
                      double pushing_objective_amplifier = 5.0,
                      int pushes_limit = 10,
                      int pushing_iteration_limit = 20,
                      int float_type = 1,
                      bool verbose = true) noexcept

{
    double ret{ 0 };

    try {
        auto ctx = std::make_shared<baryonyx::context>();
        ctx->set_logger(std::make_unique<Rcontext>());
        ctx->set_log_priority(verbose ? 6 : 3);

        auto pb = baryonyx::make_problem(ctx, file_path);
        auto mm = baryonyx::compute_min_max_objective_function(pb);
        auto type = pb.type;

        if (type == baryonyx::objective_function_type::maximize)
            ret = -std::get<0>(mm) + size(pb);
        else
            ret = std::get<1>(mm) + size(pb);

        if (type == baryonyx::objective_function_type::maximize) {
            std::swap(std::get<0>(mm), std::get<1>(mm));
            std::get<0>(mm) = -std::get<0>(mm);
            std::get<1>(mm) = -std::get<1>(mm);
        }

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
                          float_type);

        auto result = baryonyx::optimize(ctx, pb);
        ret = linearize_result(result, std::get<0>(mm), std::get<1>(mm), type);
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
