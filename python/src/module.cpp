/* Copyright (C) 2016-2026 INRAE
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

#include <baryonyx/core>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <atomic>
#include <chrono>
#include <exception>
#include <optional>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

namespace py = pybind11;
namespace bx = baryonyx;

namespace {

template<typename Enum>
struct enum_name
{
    std::string_view name;
    Enum value;
};

using order = bx::solver_parameters::constraint_order;
using norm = bx::solver_parameters::cost_norm_type;
using optimizer = bx::solver_parameters::optimizer_type;

// The names of the command line parameters.
constexpr std::array<enum_name<order>, 9> order_names{ {
  { "none", order::none },
  { "reversing", order::reversing },
  { "random-sorting", order::random_sorting },
  { "infeasibility-decr", order::infeasibility_decr },
  { "infeasibility-incr", order::infeasibility_incr },
  { "lagrangian-decr", order::lagrangian_decr },
  { "lagrangian-incr", order::lagrangian_incr },
  { "pi-sign-change", order::pi_sign_change },
  { "cycle", order::cycle },
} };

constexpr std::array<enum_name<norm>, 5> norm_names{ {
  { "none", norm::none },
  { "random", norm::random },
  { "l1", norm::l1 },
  { "l2", norm::l2 },
  { "loo", norm::loo },
} };

constexpr std::array<enum_name<optimizer>, 4> optimizer_names{ {
  { "genetic", optimizer::genetic },
  { "restart", optimizer::restart },
  { "tuned", optimizer::tuned },
  { "lns", optimizer::lns },
} };

template<typename Enum, std::size_t N>
Enum
parse_enum(const char* parameter,
           const std::array<enum_name<Enum>, N>& names,
           std::string_view value)
{
    for (const auto& elem : names)
        if (elem.name == value)
            return elem.value;

    std::string allowed;
    for (const auto& elem : names) {
        if (!allowed.empty())
            allowed += ", ";
        allowed += elem.name;
    }

    throw py::value_error(std::string(parameter) + ": unknown value '" +
                          std::string(value) + "', expected one of " +
                          allowed);
}

std::string_view
status_name(const bx::result_status status) noexcept
{
    switch (status) {
    case bx::result_status::success:
        return "success";
    case bx::result_status::internal_error:
        return "internal-error";
    case bx::result_status::uninitialized:
        return "uninitialized";
    case bx::result_status::time_limit_reached:
        return "time-limit-reached";
    case bx::result_status::kappa_max_reached:
        return "kappa-max-reached";
    case bx::result_status::limit_reached:
        return "limit-reached";
    case bx::result_status::empty_context:
        return "empty-context";
    }

    return "unknown";
}

std::string_view
file_error_name(const bx::file_format_error_tag tag) noexcept
{
    switch (tag) {
    case bx::file_format_error_tag::success:
        return "success";
    case bx::file_format_error_tag::file_not_found:
        return "file not found";
    case bx::file_format_error_tag::bad_end_of_file:
        return "bad end of file";
    case bx::file_format_error_tag::bad_general:
        return "bad general section";
    case bx::file_format_error_tag::bad_binary:
        return "bad binary section";
    case bx::file_format_error_tag::bad_objective_function_type:
        return "bad objective function type";
    case bx::file_format_error_tag::bad_objective:
        return "bad objective";
    case bx::file_format_error_tag::bad_objective_quadratic:
        return "bad quadratic objective";
    case bx::file_format_error_tag::bad_bound:
        return "bad bound";
    case bx::file_format_error_tag::bad_end:
        return "bad end";
    case bx::file_format_error_tag::bad_constraint:
        return "bad constraint";
    case bx::file_format_error_tag::too_many_variables:
        return "too many variables";
    case bx::file_format_error_tag::too_many_constraints:
        return "too many constraints";
    case bx::file_format_error_tag::bad_name:
        return "bad name";
    case bx::file_format_error_tag::bad_value:
        return "bad value";
    case bx::file_format_error_tag::empty_context:
        return "empty context";
    }

    return "unknown error";
}

std::string_view
solver_error_name(const bx::solver_error_tag tag) noexcept
{
    switch (tag) {
    case bx::solver_error_tag::no_solver_available:
        return "no solver available";
    case bx::solver_error_tag::unrealisable_constraint:
        return "unrealisable constraint";
    case bx::solver_error_tag::not_enough_memory:
        return "not enough memory";
    case bx::solver_error_tag::not_binary_variable:
        return "not binary variable";
    }

    return "unknown error";
}

std::string_view
definition_error_name(const bx::problem_definition_error_tag tag) noexcept
{
    switch (tag) {
    case bx::problem_definition_error_tag::empty_variables:
        return "empty variables";
    case bx::problem_definition_error_tag::empty_objective_function:
        return "empty objective function";
    case bx::problem_definition_error_tag::variable_not_used:
        return "variable not used";
    case bx::problem_definition_error_tag::bad_bound:
        return "bad bound";
    case bx::problem_definition_error_tag::multiple_constraint:
        return "multiple constraint";
    }

    return "unknown error";
}

/** The Python view of @c bx::result, owning its strings. */
struct result
{
    std::string status;
    bool minimize = true;
    std::optional<double> value;
    double duration = 0.0;
    long int loop = 0;
    int variables = 0;
    int constraints = 0;
    int remaining_constraints = 0;
    std::vector<double> values;
    std::optional<py::dict> solution;

    bool success() const noexcept
    {
        return status == "success";
    }
};

result
convert(const bx::result& r, const bool minimize)
{
    result ret;
    ret.status = status_name(r.status);
    ret.minimize = minimize;
    ret.duration = r.duration;
    ret.loop = r.loop;
    ret.variables = r.variables;
    ret.constraints = r.constraints;
    ret.remaining_constraints = r.remaining_constraints;

    for (const auto& elem : r.solutions)
        ret.values.emplace_back(elem.value);

    if (!r || r.solutions.empty())
        return ret;

    ret.value = r.solutions.back().value;

    // The preprocessor removes the variables it fixes from the solved
    // problem: the full assignment joins both parts.
    const auto& x = r.solutions.back().variables;
    if (x.size() != r.variable_name.size())
        return ret;

    py::dict solution;
    for (std::size_t i = 0, e = r.affected_vars.names.size(); i != e; ++i)
        solution[py::str(r.affected_vars.names[i].data(),
                         r.affected_vars.names[i].size())] =
          static_cast<int>(r.affected_vars.values[i]);
    for (std::size_t i = 0, e = x.size(); i != e; ++i)
        solution[py::str(r.variable_name[i].data(),
                         r.variable_name[i].size())] = static_cast<int>(x[i]);

    ret.solution = std::move(solution);
    return ret;
}

int
verbose_level(const py::object& verbose)
{
    // As the command line: 6 by default, 3 with --quiet.
    if (py::isinstance<py::bool_>(verbose))
        return verbose.cast<bool>() ? 6 : 3;

    return verbose.cast<int>();
}

using progress_cb = std::function<void(int, double, long int, double)>;

enum class mode
{
    solve,
    optimize
};

/**
 * Runs the solver in a worker thread: the calling thread waits with the GIL
 * released and checks the Python signals, so that Ctrl-C stops the run with
 * its best result even when no callback is called for a long time, as in
 * the search phase of lns.
 */
result
run(const mode mode,
    const std::string& file_path,
    const bx::solver_parameters& params,
    const py::object& verbose,
    const std::optional<progress_cb>& progress)
{
    auto ctx = bx::make_context(verbose_level(verbose));
    bx::context_set_solver_parameters(ctx, params);

    const auto pb = bx::make_problem(ctx, file_path);
    if (pb.status == bx::file_format_error_tag::file_not_found) {
        PyErr_SetString(PyExc_FileNotFoundError,
                        (file_path + ": cannot open the file").c_str());
        throw py::error_already_set();
    }
    if (!pb)
        throw py::value_error(file_path + ": " +
                              std::string(file_error_name(pb.status)));

    std::exception_ptr callback_error;
    if (progress) {
        bx::context_register(
          ctx,
          bx::solver_started_cb(),
          [&progress, &callback_error](int remaining,
                                       double value,
                                       long int loop,
                                       double duration,
                                       long int /*reinit*/) {
              py::gil_scoped_acquire gil;
              if (callback_error)
                  return;

              try {
                  (*progress)(remaining, value, loop, duration);
              } catch (...) {
                  callback_error = std::current_exception();
                  bx::request_stop();
              }
          });
    }

    const bool minimize = pb.type == bx::objective_function_type::minimize;
    bx::result r;
    std::exception_ptr solver_error;
    std::atomic_bool done{ false };

    std::optional<py::error_already_set> interrupted;
    {
        py::gil_scoped_release release;

        std::thread worker([&]() {
            try {
                r = mode == mode::solve ? bx::solve(ctx, pb)
                                        : bx::optimize(ctx, pb);
            } catch (...) {
                solver_error = std::current_exception();
            }
            done.store(true);
        });

        while (!done.load()) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));

            py::gil_scoped_acquire gil;
            if (!interrupted && PyErr_CheckSignals() != 0)
                interrupted.emplace();

            // solve and optimize clear the request when they start: repeat
            // it until the worker returns.
            if (interrupted || callback_error)
                bx::request_stop();
        }

        worker.join();
    }

    if (interrupted)
        throw std::move(*interrupted);
    if (callback_error)
        std::rethrow_exception(callback_error);
    if (solver_error)
        std::rethrow_exception(solver_error);

    return convert(r, minimize);
}

result
run_with(const mode mode,
         const std::string& file_path,
         const long int limit,
         const double time_limit,
         const int thread,
         const unsigned int seed,
         const double theta,
         const double delta,
         const double kappa_min,
         const double kappa_step,
         const double alpha,
         const std::string& constraint_order,
         const std::string& norm,
         const std::string& optimizer,
         const bool preprocessing,
         const py::object& verbose,
         const std::optional<progress_cb>& progress)
{
    bx::solver_parameters params;
    params.limit = limit;
    params.time_limit = time_limit;
    params.thread = thread;
    params.seed = seed;
    params.theta = theta;
    params.delta = delta;
    params.kappa_min = kappa_min;
    params.kappa_step = kappa_step;
    params.alpha = alpha;
    params.order =
      parse_enum("constraint_order", order_names, constraint_order);
    params.cost_norm = parse_enum("norm", norm_names, norm);
    params.optimizer = parse_enum("optimizer", optimizer_names, optimizer);
    params.preprocessor =
      preprocessing ? bx::solver_parameters::preprocessor_options::all
                    : bx::solver_parameters::preprocessor_options::none;

    return run(mode, file_path, params, verbose, progress);
}

} // namespace

PYBIND11_MODULE(_baryonyx, m)
{
    m.doc() = "Baryonyx: 0-1 integer linear programming with the Wedelin "
              "heuristic";

    py::class_<result>(m, "Result")
      .def_readonly("status", &result::status)
      .def_property_readonly("success", &result::success)
      .def_readonly("minimize", &result::minimize)
      .def_readonly("value", &result::value)
      .def_readonly("duration", &result::duration)
      .def_readonly("loop", &result::loop)
      .def_readonly("variables", &result::variables)
      .def_readonly("constraints", &result::constraints)
      .def_readonly("remaining_constraints", &result::remaining_constraints)
      .def_readonly("values", &result::values)
      .def_readonly("solution", &result::solution)
      .def("__bool__", &result::success)
      .def("__repr__", [](const result& r) {
          return "<baryonyx.Result status=" + r.status + " value=" +
                 (r.value ? std::to_string(*r.value) : std::string("None")) +
                 " remaining_constraints=" +
                 std::to_string(r.remaining_constraints) + ">";
      });

    const bx::solver_parameters defaults;
    const auto add = [&](const char* name, mode mode, const char* doc) {
        m.def(
          name,
          [mode](const std::string& file_path,
                 long int limit,
                 double time_limit,
                 int thread,
                 unsigned int seed,
                 double theta,
                 double delta,
                 double kappa_min,
                 double kappa_step,
                 double alpha,
                 const std::string& constraint_order,
                 const std::string& norm,
                 const std::string& optimizer,
                 bool preprocessing,
                 const py::object& verbose,
                 const std::optional<progress_cb>& progress) {
              return run_with(mode,
                              file_path,
                              limit,
                              time_limit,
                              thread,
                              seed,
                              theta,
                              delta,
                              kappa_min,
                              kappa_step,
                              alpha,
                              constraint_order,
                              norm,
                              optimizer,
                              preprocessing,
                              verbose,
                              progress);
          },
          doc,
          py::arg("file_path"),
          py::kw_only(),
          py::arg("limit") = defaults.limit,
          py::arg("time_limit") = defaults.time_limit,
          py::arg("thread") = defaults.thread,
          py::arg("seed") = defaults.seed,
          py::arg("theta") = defaults.theta,
          py::arg("delta") = defaults.delta,
          py::arg("kappa_min") = defaults.kappa_min,
          py::arg("kappa_step") = defaults.kappa_step,
          py::arg("alpha") = defaults.alpha,
          py::arg("constraint_order") = "none",
          py::arg("norm") = "l2",
          py::arg("optimizer") = "lns",
          py::arg("preprocessing") = true,
          py::arg("verbose") = py::bool_(false),
          py::arg("progress") = py::none());
    };

    add("solve",
        mode::solve,
        "Searches one solution of the 0-1 linear problem of the LP file.");
    add("optimize",
        mode::optimize,
        "Searches the best solution of the 0-1 linear problem of the LP "
        "file within the time limit.");

    py::register_exception_translator([](std::exception_ptr p) {
        try {
            if (p)
                std::rethrow_exception(p);
        } catch (const bx::file_format_failure& e) {
            std::string msg(file_error_name(e.failure()));
            if (e.line() >= 0)
                msg += " at line " + std::to_string(e.line()) +
                       ", column " + std::to_string(e.column());
            if (!e.element().empty())
                msg += ": " + e.element();
            PyErr_SetString(PyExc_ValueError, msg.c_str());
        } catch (const bx::problem_definition_failure& e) {
            std::string msg(definition_error_name(e.failure()));
            if (!e.element().empty())
                msg += ": " + e.element();
            PyErr_SetString(PyExc_ValueError, msg.c_str());
        } catch (const bx::solver_failure& e) {
            PyErr_SetString(PyExc_RuntimeError,
                            std::string(solver_error_name(e.failure())).c_str());
        }
    });

    m.def("request_stop",
          &bx::request_stop,
          "Asks the running solve or optimize to return its best result.");
}
