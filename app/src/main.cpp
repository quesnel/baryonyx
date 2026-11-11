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

#include <baryonyx/core-out>
#include <baryonyx/core-utils>

#include "main.hpp"

#include <fmt/chrono.h>
#include <fmt/color.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include <chrono>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>

#include <cassert>
#include <cctype>
#include <csignal>
#include <cstring>
#include <ctime>

#ifndef _WIN32
#include <unistd.h>
#else
#include <windows.h>
#endif

#ifndef _WIN32
pid_t
get_pid() noexcept
{
    return ::getpid();
}
#else
DWORD
get_pid() noexcept
{
    return ::GetCurrentProcessId();
}
#endif

static std::tm
local_time_now() noexcept
{
    const auto now =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    std::tm tm{};
#ifdef _WIN32
    ::localtime_s(&tm, &now);
#else
    ::localtime_r(&now, &tm);
#endif
    return tm;
}

static constexpr std::string_view
to_string(baryonyx::solver_parameters::optimizer_type type) noexcept
{
    switch (type) {
    case baryonyx::solver_parameters::optimizer_type::genetic:
        return "genetic";
    case baryonyx::solver_parameters::optimizer_type::restart:
        return "restart";
    case baryonyx::solver_parameters::optimizer_type::tuned:
        return "tuned";
    case baryonyx::solver_parameters::optimizer_type::lns:
        return "lns";
    }

    return "unknown";
}

static void
solver_started_cb(const baryonyx::solver_parameters& params)
{
    fmt::print("Solver starts\n");

    fmt::print(" * Global parameters:\n"
               "  - limit: {}\n"
               "  - time-limit: {:.10g}s\n"
               "  - auto-tune: {}\n"
               "  - observation: {}\n"
               "  - optimizer: {}\n",
               params.limit,
               params.time_limit,
               params.mode,
               params.observer,
               to_string(params.optimizer));

    fmt::print(" * In The Middle parameters:\n"
               "  - preprocessing: {}\n"
               "  - constraint-order: {}\n"
               "  - theta: {:.10g}\n"
               "  - delta: {:.10g}\n"
               "  - kappa: {:.10g} {:.10g} {:.10g}\n"
               "  - alpha: {:.10g}\n"
               "  - w: {:.10g}\n"
               "  - norm: {}\n",
               params.pre_order,
               params.order,
               params.theta,
               params.delta,
               params.kappa_min,
               params.kappa_step,
               params.kappa_max,
               params.alpha,
               params.w,
               params.cost_norm);

    fmt::print(" * Pushes system parameters:\n"
               "  - pushes-limit: {}\n"
               "  - pushing-objective-amplifier: {:.10g}\n"
               "  - pushing-iteration-limit: {}\n"
               "  - pushing-k-factor: {:.10g}\n",
               params.pushes_limit,
               params.pushing_objective_amplifier,
               params.pushing_iteration_limit,
               params.pushing_k_factor);

    fmt::print(" * Solver initialization parameters:\n"
               "  - init-policy-random: {}\n",
               params.init_policy_random);
}

static void
solver_updated_cb(int remaining_constraints,
                  double value,
                  long int loop,
                  double duration,
                  long int reinit_number)
{
    if (remaining_constraints > 0) {
        fmt::print(
          "  - Constraints remaining: {} (loop: {} t: {}s reinit: {})\n",
          remaining_constraints,
          loop,
          duration,
          reinit_number);
    } else {
        if (loop >= 0)
            fmt::print(
              "  - Solution found: {:f} (loop: {} t: {}s reinit: {})\n",
              value,
              loop,
              duration,
              reinit_number);
        else
            fmt::print("  - Solution found via push: {:f} (loop: {} t: {}s "
                       "reinit: {})\n",
                       value,
                       -loop,
                       duration,
                       reinit_number);
    }
}

static void
solver_finished_cb(const baryonyx::result& r)
{
    fmt::print("Solver finished\n");

    switch (r.status) {
    case baryonyx::result_status::success:
        if (!r.solutions.empty()) {
            if (r.solutions.back().variables.empty()) {
                fmt::print(
                  "Best solution found via preprocessor: {:.10g} in {}s\n",
                  r.solutions.back().value,
                  r.duration);
            } else {
                if (r.loop >= 0)
                    fmt::print(
                      "Best solution found: {:.10g} in {} loop and {}s\n",
                      r.solutions.back().value,
                      r.loop,
                      r.duration);
                else
                    fmt::print("Best solution found via push: {:.10g} in {} "
                               "loop and {}s\n",
                               r.solutions.back().value,
                               -r.loop,
                               r.duration);
            }
        }
        break;
    case baryonyx::result_status::internal_error:
        fmt::print("No solution. Internal error\n");
        break;
    case baryonyx::result_status::uninitialized:
        fmt::print("No solution. Uninitialized error\n");
        break;
    case baryonyx::result_status::kappa_max_reached:
        fmt::print(
          "No solution. Constraint remaining: {}. Kappa reached in {}s.\n",
          r.remaining_constraints,
          r.duration);
        break;
    case baryonyx::result_status::time_limit_reached:
        fmt::print("No solution. Constraint remaining: {}. Time limit reached "
                   "at {}s.\n",
                   r.remaining_constraints,
                   r.duration);
        break;
    case baryonyx::result_status::limit_reached:
        fmt::print("No solution. Constraint remaining: {}. Loop limit reached "
                   "in {}s.\n",
                   r.remaining_constraints,
                   r.duration);
        break;
    case baryonyx::result_status::empty_context:
        fmt::print("Context uninitialized\n");
        break;
    }
}

/**
 * @brief Tries to split the argument between name and value.
 * @details The split consists to return left and right part arround the
 *     characters `=' or `:'.
 *
 * @param param String to split.
 * @return Return left and right part in success, otherwise only the left
 * part if the splitting characters are not found or an empty tuple if all
 * is empty.
 */
constexpr static std::tuple<std::string_view, std::string_view>
split_argument(const std::string_view param)
{
    auto position = param.find_first_of(":=");

    if (position == std::string_view::npos) /* If nothing is found, return
                                               an empty tuple.*/
        return {};

    if (position + 1 >= param.size())
        return std::make_tuple(param, std::string_view{});

    return std::make_tuple(param.substr(0, position),
                           param.substr(position + 1));
}

struct get_param
{
    constexpr get_param(int argc_, const char** argv_)
      : argv(argv_)
      , argc(argc_)
      , i(0)
    {
    }

    constexpr std::optional<std::string_view> operator()(
      int arg_position,
      const std::string_view longp,
      const std::string_view shortp) noexcept
    {
        std::string_view arg(argv[arg_position]);
        i = arg_position;

        if (arg.starts_with(longp) && arg.size() > longp.size() &&
            (arg[longp.size()] == '=' || arg[longp.size()] == ':'))
            return arg.substr(longp.size() + 1);

        if (arg.starts_with(shortp) && arg.size() > shortp.size())
            return arg.substr(shortp.size());

        if (arg_position + 1 < argc) {
            if (arg == longp) {
                i = arg_position + 1;
                return argv[i];
            }

            if (arg == shortp) {
                i = arg_position + 1;
                return argv[i];
            }
        }

        return std::nullopt;
    }

    const char** argv;
    int argc;
    int i;
};

static void
help() noexcept
{
    fmt::print(
      "Baryonyx v{}.{}.{}", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);

    fmt::print(
      "\nUsage: baryonyx [options] file.lp\n"
      "\nGeneral options:\n"
      "  --help|-h                     This help message\n"
      "  --optimize|-O                 Optimize the objective (default: "
      "feasibility search only)\n"
      "  --param|-p name[:|=]value     Set a solver parameter (see below)\n"
      "  --auto|-a branch              Split the problem on variables\n"
      "  --check|-C file.sol           Check a solution file and print its "
      "value\n"
      "  --disable-preprocessing|-np   Disable preprocessing\n"
      "  --verbose|-v int              Verbose level [0, 7] (default: 6)\n"
      "  --quiet|-q                    Remove most verbose messages\n"
      "  --bench|-b name               Benchmark mode, store results in the "
      "name column of the csv files\n"
      "  --debug                       Enable the debug mode\n"
      "\nParameters (-p name=value)\n"
      " * Global parameters\n"
      "  - limit: integer [-1, +oo[ loop number, <= 0 means no limit "
      "(default: 1000)\n"
      "  - time-limit: real in seconds, <= 0 means no limit\n"
      "  - thread: integer [0, +oo[, 0 means automatic\n"
      "  - seed: integer [0, +oo[, 0 means random\n"
      "  - observer-type: none (default) pnm file\n"
      " * In The Middle parameters\n"
      "  - theta: real [0, 1] (default: 0.5)\n"
      "  - delta: real ]0, +oo[ (default: automatic)\n"
      "  - kappa-min: real [0, 1] (default: 0)\n"
      "  - kappa-step: real [0, 1] (default: 0.001)\n"
      "  - alpha: real [0, 2] (default: 1)\n"
      "  - w: real [0, +oo[ (default: 0.05)\n"
      "  - norm: loo (default) l1 l2 random none\n"
      "  - constraint-order: none (default) reversing random-sorting "
      "infeasibility-decr infeasibility-incr lagrangian-decr "
      "lagrangian-incr pi-sign-change cycle\n"
      "  - preprocessing: none (default) memory less-greater-equal "
      "less-equal-greater greater-less-equal greater-equal-less "
      "equal-less-greater equal-greater-less\n"
      "  - buffer-solver: 0 (default) or 1 (experimental parallel -1 0 1 "
      "solver)\n"
      " * Initialization parameters\n"
      "  - init-policy-random: real [0, 1] (default: 0.5) probability to "
      "invert a variable of the initial solution\n"
      " * Optimizer parameters\n"
      "  - optimizer: lns (default) restart genetic tuned\n");
}

constexpr static bool
is_equal(std::string_view name, const char* longf, char shortf = '\0')
{
    if (name == longf)
        return true;

    if (shortf != '\0' && name.size() == 2 && name[1] == shortf)
        return true;

    return false;
}

constexpr static std::optional<double>
assign_0oo(std::string_view value)
{
    auto result = ::to_double(value);

    if (result.has_value() && *result > 0)
        return *result;
    else
        return std::nullopt;
}

constexpr static std::optional<double>
assign_01(std::string_view value)
{
    auto result = ::to_double(value);

    if (result.has_value() && *result >= 0 && *result <= 1)
        return *result;
    else
        return std::nullopt;
}

constexpr static std::optional<int>
assign(std::string_view value, int mindef, int maxdef)
{
    auto result = ::to_int(value);

    if (result.has_value() && *result >= mindef && *result <= maxdef)
        return result;
    else
        return std::nullopt;
}

constexpr static std::optional<long int>
assign(std::string_view value, long int mindef, long int maxdef)
{
    auto result = ::to_long(value);

    if (result.has_value() && *result >= mindef && *result <= maxdef)
        return result;
    else
        return std::nullopt;
}

constexpr static std::optional<double>
assign_d(std::string_view value, double mindef, double maxdef)
{
    auto result = ::to_double(value);

    if (result.has_value() && *result >= mindef && *result <= maxdef)
        return *result;
    else
        return std::nullopt;
}

enum class command_line_status
{
    success,
    unknown,
    unknown_option,
    parameter_missing,
    verbose_error,
    limit_error,
    time_limit_error,
    observer_type_error,
    preprocessing_error,
    constraint_order_error,
    theta_error,
    delta_error,
    kappa_min_error,
    kappa_step_error,
    alpha_error,
    w_error,
    norm_error,
    init_policy_random_error,
    buffer_solver_error,
    optimizer_error,
    thread_error,
    seed_error
};

constexpr const std::string_view command_line_status_string[] = {
    "success",
    "unknown",
    "unknown_option",
    "parameter_missing",
    "verbose",
    "limit",
    "time_limit",
    "observer_type",
    "preprocessing",
    "constraint_order",
    "theta",
    "delta",
    "kappa_min",
    "kappa_step",
    "alpha",
    "w",
    "norm",
    "init_policy_random",
    "buffer_solver",
    "optimizer",
    "thread",
    "seed"
};

static_assert(std::size(command_line_status_string) ==
              static_cast<std::size_t>(command_line_status::seed_error) + 1);

constexpr std::string_view
to_string(command_line_status s) noexcept
{
    const auto x = static_cast<std::underlying_type_t<command_line_status>>(s);

    assert(x >= 0 && std::cmp_less(x, std::size(command_line_status_string)));

    return command_line_status_string[x];
}

constexpr static command_line_status
assign_parameter(baryonyx::solver_parameters& params,
                 std::string_view name,
                 std::string_view value)
{
    if (is_equal(name, "limit", 'l')) {
        if (auto v =
              assign(value, -1L, (std::numeric_limits<long int>::max)());
            !v)
            return command_line_status::limit_error;
        else
            params.limit = *v;

    } else if (is_equal(name, "time-limit")) {
        if (auto v = to_double(value); !v)
            return command_line_status::time_limit_error;
        else
            params.time_limit = *v <= 0.0 ? -1.0 : *v;

    } else if (is_equal(name, "observer-type")) {
        if (value == "none")
            params.observer = baryonyx::solver_parameters::observer_type::none;
        else if (value == "pnm")
            params.observer = baryonyx::solver_parameters::observer_type::pnm;
        else if (value == "file")
            params.observer = baryonyx::solver_parameters::observer_type::file;
        else
            return command_line_status::observer_type_error;

    } else if (is_equal(name, "preprocessing")) {
        if (value == "none")
            params.pre_order =
              baryonyx::solver_parameters::pre_constraint_order::none;
        else if (value == "memory")
            params.pre_order =
              baryonyx::solver_parameters::pre_constraint_order::memory;
        else if (value == "less-greater-equal")
            params.pre_order = baryonyx::solver_parameters::
              pre_constraint_order::less_greater_equal;
        else if (value == "less-equal-greater")
            params.pre_order = baryonyx::solver_parameters::
              pre_constraint_order::less_equal_greater;
        else if (value == "greater-less-equal")
            params.pre_order = baryonyx::solver_parameters::
              pre_constraint_order::greater_less_equal;
        else if (value == "greater-equal-less")
            params.pre_order = baryonyx::solver_parameters::
              pre_constraint_order::greater_equal_less;
        else if (value == "equal-less-greater")
            params.pre_order = baryonyx::solver_parameters::
              pre_constraint_order::equal_less_greater;
        else if (value == "equal-greater-less")
            params.pre_order = baryonyx::solver_parameters::
              pre_constraint_order::equal_greater_less;
        else
            return command_line_status::preprocessing_error;

    } else if (is_equal(name, "constraint-order")) {
        if (value == "none")
            params.order = baryonyx::solver_parameters::constraint_order::none;
        else if (value == "reversing")
            params.order =
              baryonyx::solver_parameters::constraint_order::reversing;
        else if (value == "random-sorting")
            params.order =
              baryonyx::solver_parameters::constraint_order::random_sorting;
        else if (value == "infeasibility-decr")
            params.order = baryonyx::solver_parameters::constraint_order::
              infeasibility_decr;
        else if (value == "infeasibility-incr")
            params.order = baryonyx::solver_parameters::constraint_order::
              infeasibility_incr;
        else if (value == "lagrangian-decr")
            params.order =
              baryonyx::solver_parameters::constraint_order::lagrangian_decr;
        else if (value == "lagrangian-incr")
            params.order =
              baryonyx::solver_parameters::constraint_order::lagrangian_incr;
        else if (value == "pi-sign-change")
            params.order =
              baryonyx::solver_parameters::constraint_order::pi_sign_change;
        else if (value == "cycle")
            params.order =
              baryonyx::solver_parameters::constraint_order::cycle;
        else
            return command_line_status::constraint_order_error;

    } else if (is_equal(name, "theta")) {
        if (auto v = assign_01(value); !v)
            return command_line_status::theta_error;
        else
            params.theta = *v;

    } else if (is_equal(name, "delta")) {
        if (auto v = assign_0oo(value); !v)
            return command_line_status::delta_error;
        else
            params.delta = *v;

    } else if (is_equal(name, "kappa-min")) {
        if (auto v = assign_01(value); !v)
            return command_line_status::kappa_min_error;
        else
            params.kappa_min = *v;

    } else if (is_equal(name, "kappa-step")) {
        if (auto v = assign_01(value); !v)
            return command_line_status::kappa_step_error;
        else
            params.kappa_step = *v;

    } else if (is_equal(name, "alpha")) {
        if (auto v = assign_d(value, 0, 2); !v)
            return command_line_status::alpha_error;
        else
            params.alpha = *v;

    } else if (is_equal(name, "w")) {
        if (auto v =
              assign_d(value, 0.0, (std::numeric_limits<double>::max)());
            !v)
            return command_line_status::w_error;
        else
            params.w = *v;

    } else if (is_equal(name, "norm")) {
        if (value == "none")
            params.cost_norm =
              baryonyx::solver_parameters::cost_norm_type::none;
        else if (value == "random")
            params.cost_norm =
              baryonyx::solver_parameters::cost_norm_type::random;
        else if (value == "l1")
            params.cost_norm = baryonyx::solver_parameters::cost_norm_type::l1;
        else if (value == "l2")
            params.cost_norm = baryonyx::solver_parameters::cost_norm_type::l2;
        else if (value == "loo")
            params.cost_norm =
              baryonyx::solver_parameters::cost_norm_type::loo;
        else
            return command_line_status::norm_error;

    } else if (is_equal(name, "init-policy-random")) {
        if (auto v = assign_01(value); !v)
            return command_line_status::init_policy_random_error;
        else
            params.init_policy_random = *v;

    } else if (is_equal(name, "buffer-solver")) {
        if (auto v = assign(value, 0, 1); !v)
            return command_line_status::buffer_solver_error;
        else
            params.use_buffer_solver = *v == 1;

    } else if (is_equal(name, "optimizer")) {
        if (value == "genetic")
            params.optimizer =
              baryonyx::solver_parameters::optimizer_type::genetic;
        else if (value == "restart")
            params.optimizer =
              baryonyx::solver_parameters::optimizer_type::restart;
        else if (value == "tuned")
            params.optimizer =
              baryonyx::solver_parameters::optimizer_type::tuned;
        else if (value == "lns")
            params.optimizer =
              baryonyx::solver_parameters::optimizer_type::lns;
        else
            return command_line_status::optimizer_error;

    } else if (is_equal(name, "thread")) {
        if (auto v = assign(value, 0, (std::numeric_limits<int>::max)()); !v)
            return command_line_status::thread_error;
        else
            params.thread = *v;

    } else if (is_equal(name, "seed")) {
        if (auto v = assign(value, 0, (std::numeric_limits<int>::max)()); !v)
            return command_line_status::seed_error;
        else
            params.seed = *v;

    } else
        return command_line_status::unknown;

    return command_line_status::success;
}

struct main_parameters
{
    baryonyx::solver_parameters parameters;
    std::vector<std::string> filenames;
    std::string check_filename;
    std::string bench_name;
    int verbose = 6;
    int parse_arg = 0;
    bool check = false;
    bool optimize = false;
    bool quiet = false;
    bool bench = false;
    command_line_status parse_status = command_line_status::success;
};

static main_parameters
parse(int argc, const char* argv[])
{
    main_parameters ret;
    get_param get(argc, argv);

    for (int i = 1; i < argc; ++i) {
        std::string_view arg(argv[i]);

        ret.parse_arg = i;

        if (arg == "--help" || arg == "-h") {
            ::help();
            continue;
        }

        if (arg == "--quiet" || arg == "-q") {
            ret.quiet = true;
            continue;
        }

        if (arg == "--optimize" || arg == "-O") {
            ret.optimize = true;
            continue;
        }

        if (auto opt = get(i, "--bench", "-b"); opt) {
            ret.bench = true;
            ret.bench_name = *opt;
            i = get.i;
            continue;
        }

        if (arg == "--disable-preprocessing" || arg == "-np") {
            ret.parameters.preprocessor =
              baryonyx::solver_parameters::preprocessor_options::none;
            continue;
        }

        if (auto opt = get(i, "--auto", "-a"); opt) {
            if (*opt == "branch")
                ret.parameters.mode =
                  baryonyx::solver_parameters::mode_type::branch;
            else {
                ret.parse_status = command_line_status::unknown_option;
                return ret;
            }

            i = get.i;
            continue;
        }

        if (auto opt = get(i, "--verbose", "-v"); opt) {
            if (auto v = ::assign(*opt, 0, 7); v) {
                ret.verbose = *v;
            } else {
                ret.parse_status = command_line_status::verbose_error;
                return ret;
            }

            i = get.i;
            continue;
        }

        if (auto opt = get(i, "--check", "-C"); opt) {
            ret.check_filename = *opt;
            ret.check = true;
            i = get.i;
            continue;
        }

        if (auto opt = get(i, "--param", "-p"); opt) {
            ret.parse_arg = get.i;
            auto [name, value] = split_argument(*opt);

            if (name.empty() || value.empty()) {
                ret.parse_status = command_line_status::parameter_missing;
                return ret;
            }

            if (auto v = assign_parameter(ret.parameters, name, value);
                v != command_line_status::success) {
                ret.parse_status = v;
                return ret;
            }

            i = get.i;
            continue;
        }

        if (arg == "--debug") {
            ret.parameters.debug = true;
            continue;
        }

        if (arg.size() > 1 && arg.front() == '-') {
            ret.parse_status = command_line_status::unknown_option;
            return ret;
        }

        ret.filenames.emplace_back(argv[i]);
    }

    return ret;
}

static void
resume(const baryonyx::raw_problem& pb) noexcept
{
    int real{ 0 }, general{ 0 }, binary{ 0 };

    for (const auto& var : pb.vars.values) {
        switch (var.type) {
        case baryonyx::variable_type::real:
            ++real;
            break;
        case baryonyx::variable_type::binary:
            ++binary;
            break;
        case baryonyx::variable_type::general:
            ++general;
            break;
        }
    }

    fmt::print("  - objective: {}\n"
               "  - mode: {}\n"
               "  - variables: {}/{}/{} (real/general/binary)\n"
               "  - constraints: {}/{}/{} (equal/greater/less)\n",
               (pb.type == baryonyx::objective_function_type::maximize
                  ? "maximize"
                  : "minimize"),
               (pb.objective.qelements.empty() ? "linear" : "quadratic"),
               real,
               general,
               binary,
               pb.equal_constraints.size(),
               pb.greater_constraints.size(),
               pb.less_constraints.size());
}

static void
resume(const baryonyx::raw_problem& pb, std::ostream& os) noexcept
{
    int real{ 0 }, general{ 0 }, binary{ 0 };

    for (const auto& var : pb.vars.values) {
        switch (var.type) {
        case baryonyx::variable_type::real:
            ++real;
            break;
        case baryonyx::variable_type::binary:
            ++binary;
            break;
        case baryonyx::variable_type::general:
            ++general;
            break;
        }
    }

    fmt::print(os,
               "\\ objective : {}\n"
               "\\ mode: {}\n"
               "\\ variables: {}/{}/{} (real/general/binary)\n"
               "\\ constraints: {}/{}/{} (equal/greater/less)\n",
               (pb.type == baryonyx::objective_function_type::maximize
                  ? "maximize"
                  : "minimize"),
               (pb.objective.qelements.empty() ? "linear" : "quadratic"),
               real,
               general,
               binary,
               pb.equal_constraints.size(),
               pb.greater_constraints.size(),
               pb.less_constraints.size());
}

static void
resume(const baryonyx::result& result, std::ostream& os) noexcept
{
    fmt::print(os,
               "\\ solver................: {}\n"
               "\\ constraints...........: {}\n"
               "\\ variables.............: {}\n"
               "\\ duration..............: {:f}\n"
               "\\ loop..................: {}\n"
               "\\ status................: ",
               result.method,
               result.constraints,
               result.variables,
               result.duration,
               result.loop);

    switch (result.status) {
    case baryonyx::result_status::internal_error:
        fmt::print(os, "internal error reported\n");
        break;
    case baryonyx::result_status::uninitialized:
        fmt::print(os, "uninitialized\n");
        break;
    case baryonyx::result_status::success:
        fmt::print(os, "solution found\n");

        if (result.solutions.empty()) // Baryonyx ensures solutions are not
            break;                    // empty.

        fmt::print(os,
                   "\\ value.................: {:f}\n"
                   "\\ other value...........: ",
                   result.solutions.back().value);

        for (const auto& elem : result.solutions)
            fmt::print(os, "{} ", elem.value);
        fmt::print(os, "\n");
        fmt::print(os, "\\ variables.............: \n");

        for (std::size_t i{ 0 }, e{ result.affected_vars.names.size() };
             i != e;
             ++i)
            fmt::print(os,
                       "{}={}\n",
                       result.affected_vars.names[i],
                       (result.affected_vars.values[i] ? 1 : 0));

        for (std::size_t i{ 0 }, e{ result.solutions.back().variables.size() };
             i != e;
             ++i)
            fmt::print(os,
                       "{}={}\n",
                       result.variable_name[i],
                       (result.solutions.back().variables[i] ? 1 : 0));

        break;
    case baryonyx::result_status::time_limit_reached:
        fmt::print(os,
                   "time limit reached\n"
                   "\\ remaining constraints.: {}\n",
                   result.remaining_constraints);
        break;
    case baryonyx::result_status::kappa_max_reached:
        fmt::print(os,
                   "kappa max reached\n"
                   "\\ remaining constraints.: {}\n",
                   result.remaining_constraints);
        break;
    case baryonyx::result_status::limit_reached:
        fmt::print(os,
                   "limit reached\n"
                   "\\ remaining constraints.: {}\n",
                   result.remaining_constraints);
        break;
    case baryonyx::result_status::empty_context:
        fmt::print(os, "context uninitialized\n");
        break;
    }
}

// The first Ctrl-C (or SIGTERM) stops the computation, which then returns
// and writes its best solution; a second one kills the process as usual.
extern "C" void
stop_on_signal(int signal_number)
{
    baryonyx::request_stop();
    std::signal(signal_number, SIG_DFL);
}

static inline auto
solve_or_optimize(const baryonyx::context_ptr& ctx,
                  const baryonyx::raw_problem& pb,
                  bool with_optimization) -> baryonyx::result
{
    std::signal(SIGINT, stop_on_signal);
    std::signal(SIGTERM, stop_on_signal);

    auto result = with_optimization ? baryonyx::optimize(ctx, pb)
                                    : baryonyx::solve(ctx, pb);

    std::signal(SIGINT, SIG_DFL);
    std::signal(SIGTERM, SIG_DFL);

    return result;
}

static bool
check_solution(const baryonyx::context_ptr& ctx,
               const std::string& problem_filename,
               const std::string& solution_filename)
{
    try {
        const auto pb = baryonyx::make_problem(ctx, problem_filename);
        if (!pb) {
            fmt::print(stderr, "Fail to read file: {}\n", pb.status);
            return false;
        }

        const auto result = baryonyx::make_result(ctx, solution_filename);
        if (result.status == baryonyx::result_status::time_limit_reached) {
            fmt::print("Check {} with {}: no solution in file\n",
                       problem_filename,
                       solution_filename);
            return false;
        }

        if (!result) {
            fmt::print(stderr,
                       "Fail to read solution file: {}\n",
                       solution_filename);
            return false;
        }

        const bool valid = baryonyx::is_valid_solution(pb, result);
        fmt::print("Check {} with {}: {} (value: {:.10g})\n",
                   problem_filename,
                   solution_filename,
                   valid ? "success" : "failure",
                   baryonyx::compute_solution(pb, result));

        return valid;
    } catch (const baryonyx::file_format_failure& e) {
        fmt::print(stderr,
                   "check failure: {} at line {}: `{}'\n",
                   e.failure(),
                   e.line(),
                   e.element());
        return false;
    } catch (const std::exception& e) {
        fmt::print(stderr, "check failure: {}\n", e.what());
        return false;
    }
}

// Message of the exception being handled, for the console and the result
// files.
static std::string
current_failure_message()
{
    try {
        throw;
    } catch (const baryonyx::precondition_failure& e) {
        return fmt::format("internal failure: {}", e.what());
    } catch (const baryonyx::postcondition_failure& e) {
        return fmt::format("internal failure: {}", e.what());
    } catch (const baryonyx::numeric_cast_failure& e) {
        return fmt::format("numeric cast internal failure: {}", e.what());
    } catch (const baryonyx::file_access_failure& e) {
        return fmt::format("file `{}' fail {}: {}",
                           e.file(),
                           e.error(),
                           std::strerror(e.error()));
    } catch (const baryonyx::file_format_failure& e) {
        return fmt::format("file format error at line {} column {} {}",
                           e.line(),
                           e.column(),
                           e.failure());
    } catch (const baryonyx::problem_definition_failure& e) {
        return fmt::format(
          "definition problem error at {}: {}", e.element(), e.failure());
    } catch (const baryonyx::solver_failure& e) {
        return fmt::format("solver error: {}", e.failure());
    } catch (const std::exception& e) {
        return fmt::format("failure: {}", e.what());
    } catch (...) {
        return "unknown failure";
    }
}

int
main(int argc, const char* argv[])
{
    auto params = parse(argc, argv);

    if (params.parse_status == command_line_status::unknown_option) {
        fmt::print(stderr,
                   fmt::emphasis::bold | fmt::fg(fmt::terminal_color::red),
                   "Unknown option `{}`, see --help\n",
                   argv[params.parse_arg]);
        return EXIT_FAILURE;
    }

    if (params.parse_status != command_line_status::success) {
        fmt::print(stderr,
                   fmt::emphasis::bold | fmt::fg(fmt::terminal_color::red),
                   "Parameter `{}` unknown value `{}`\n",
                   to_string(params.parse_status),
                   argv[params.parse_arg]);
        return EXIT_FAILURE;
    }

    if (params.filenames.empty()) {
        fmt::print(stderr,
                   fmt::emphasis::bold | fmt::fg(fmt::terminal_color::red),
                   "Missing lp or benchmark file(s), see --help\n");
        return EXIT_SUCCESS;
    }

    params.verbose = params.quiet ? 3 : params.verbose;
    auto ctx = baryonyx::make_context(params.verbose);

    context_set_solver_parameters(ctx, params.parameters);

    if (params.check) {
        if (params.filenames.size() != 1) {
            fmt::print(stderr, "--check needs exactly one lp file\n");
            return EXIT_FAILURE;
        }

        return check_solution(
                 ctx, params.filenames.front(), params.check_filename)
                 ? EXIT_SUCCESS
                 : EXIT_FAILURE;
    }

    context_register(
      ctx, solver_started_cb, solver_updated_cb, solver_finished_cb);

    if (params.bench) {
        for (auto filename : params.filenames) {
            if (params.bench_name.empty())
                params.bench_name = fmt::format(
                  "bx-{}.{}.{}", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);

            fmt::print("Start benchmark for file {} (name: {})\n",
                       filename,
                       params.bench_name);

            if (!benchmark(ctx, filename, params.bench_name))
                fmt::print(stderr,
                           fmt::emphasis::bold |
                             fmt::fg(fmt::terminal_color::red),
                           "Benchmark for {} failed.\n",
                           filename);
        }
    } else {
        if (params.filenames.size() == 1) {
            auto pb = baryonyx::make_problem(ctx, params.filenames.front());
            if (!pb) {
                fmt::print(stderr, "Fail to read file: {}\n", pb.status);
            } else {
                std::ofstream ofs;

                try {
                    auto filename = fmt::format(
                      "{}-{}.sol", params.filenames.front(), get_pid());
                    fmt::print("  - output file: {}\n", filename);

                    resume(pb);
                    ofs.open(filename);

                    resume(pb, ofs);

                    fmt::print(ofs,
                               "\\ solver starts: {:%Y-%m-%d %X}\n",
                               local_time_now());

                    auto ret = solve_or_optimize(ctx, pb, params.optimize);

                    fmt::print(ofs,
                               "\\ solver finishes: {:%Y-%m-%d %X}\n",
                               local_time_now());

                    if (ret.status == baryonyx::result_status::success) {
                        fmt::print(ofs,
                                   "\\ Solution found: {:f}\n",
                                   ret.solutions.back().value);
                        resume(ret, ofs);
                    } else {
                        fmt::print(ofs,
                                   "\\ Solution not found. Missing "
                                   "constraints: {}\n",
                                   ret.remaining_constraints);
                    }
                } catch (...) {
                    // The solution file tells why it holds no solution.
                    const auto message = current_failure_message();
                    fmt::print(stderr, "{}\n", message);

                    if (ofs.is_open())
                        fmt::print(ofs,
                                   "\\ Solution not found. Error: {}\n",
                                   message);
                }
            }
        } else {
            auto filename = fmt::format("baryonyx-{}.res", get_pid());
            std::ofstream ofs(filename);

            for (auto& elem : params.filenames) {
                auto pb = baryonyx::make_problem(ctx, elem);
                if (!pb) {
                    fmt::print(stderr, "Fail to read file: {}\n", pb.status);
                    continue;
                }

                try {
                    fmt::print(ofs, "{} ", elem);

                    auto ret = solve_or_optimize(ctx, pb, params.optimize);

                    if (ret.status == baryonyx::result_status::success) {
                        fmt::print(ofs,
                                   "{:f} {}s ",
                                   ret.solutions.back().value,
                                   ret.duration);

                        for (const auto& sol : ret.solutions)
                            fmt::print(ofs, "{} ", sol.value);
                        fmt::print(ofs, "\n");

                    } else {
                        fmt::print(ofs, "No solution found.\n");
                    }
                    ofs.flush();
                } catch (...) {
                    const auto message = current_failure_message();
                    fmt::print(stderr, "{}\n", message);
                    fmt::print(ofs, "error: {}\n", message);
                }
            }
        }
    }

    return EXIT_SUCCESS;
}
