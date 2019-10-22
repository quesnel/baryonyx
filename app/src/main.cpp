/* Copyright (C) 2016-2019 INRA
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

#include <cctype>
#include <cstring>

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

static void
solver_started_cb(const baryonyx::solver_parameters& params)
{
    fmt::print("Solver starts\n");

    fmt::print(" * Global parameters:\n"
               "  - limit: {}\n"
               "  - time-limit: {:.10g}s\n"
               "  - floating-point-type: {}\n"
               "  - print-level: {}\n"
               "  - auto-tune: {}\n"
               "  - observation: {}\n",
               params.limit,
               params.time_limit,
               baryonyx::floating_point_type_to_string(params.float_type),
               params.print_level,
               baryonyx::mode_type_to_string(params.mode),
               baryonyx::observer_type_to_string(params.observer));

    if (params.solver == baryonyx::solver_parameters::solver_type::bastert) {
        fmt::print(" * In The Middle parameters:\n"
                   "  - preprocessing: {}\n"
                   "  - constraint-order: {}\n"
                   "  - theta: {:.10g}\n"
                   "  - delta: {:.10g}\n"
                   "  - kappa: {:.10g} {:.10g} {:.10g}\n"
                   "  - alpha: {:.10g}\n"
                   "  - w: {:.10g}\n"
                   "  - norm: {}\n",
                   baryonyx::pre_constraint_order_to_string(params.pre_order),
                   baryonyx::constraint_order_to_string(params.order),
                   params.theta,
                   params.delta,
                   params.kappa_min,
                   params.kappa_step,
                   params.kappa_max,
                   params.alpha,
                   params.w,
                   baryonyx::cost_norm_type_to_string(params.cost_norm));

        fmt::print(" * Pushes system parameters:\n"
                   "  - pushes-limit: {}\n"
                   "  - pushing-objective-amplifier: {:.10g}\n"
                   "  - pushing-iteration-limit: {}\n"
                   "  - pushing-k-factor: {:.10g}\n",
                   params.pushes_limit,
                   params.pushing_objective_amplifier,
                   params.pushing_iteration_limit,
                   params.pushing_k_factor);

        fmt::print(" * Initialization parameters:\n"
                   "  - init-policy: {}\n"
                   "  - init-policy-random: {}\n"
                   "  - init-random: {:.10g}\n",
                   baryonyx::init_policy_type_to_string(params.init_policy),
                   params.init_policy_random,
                   params.init_random);
    } else {
        fmt::print(" * Random solvers:\n"
                   "  - random: bernouilli with p=0.5\n");
    }
}

static void
solver_updated_cb(const baryonyx::result& r)
{
    if (r.status != baryonyx::result_status::success) {
        fmt::print("  - Constraints remaining: {} (loop: {} t: {}s)\n",
                   r.remaining_constraints,
                   r.loop,
                   r.duration);
    } else {
        if (r.loop >= 0)
            fmt::print("  - Solution found: {:f} (loop: {} t: {}s)\n",
                       r.solutions.back().value,
                       r.loop,
                       r.duration);
        else
            fmt::print("  - Solution found via push: {:f} (loop: {} t: {}s)\n",
                       r.solutions.back().value,
                       -r.loop,
                       r.duration);
    }
}

static void
solver_finished_cb(const baryonyx::result& r)
{
    fmt::print("Solver finished\n");

    switch (r.status) {
    case baryonyx::result_status::success:
        if (r.loop >= 0)
            fmt::print("Best solution found: {:.10g} in {} loop and {}s\n",
                       r.solutions.back().value,
                       r.loop,
                       r.duration);
        else
            fmt::print(
              "Best solution found via push: {:.10g} in {} loop and {}s\n",
              r.solutions.back().value,
              -r.loop,
              r.duration);
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
    }
}

/**
 * @brief Tries to split the argument between name and value.
 * @details The split consists to return left and right part arround the
 *     characters `=' or `:'.
 *
 * @param param String to split.
 * @return Return left and right part in success, otherwise only the left part
 *     if the splitting characters are not found or an empty tuple if all is
 *     empty.
 */
constexpr static std::tuple<std::string_view, std::string_view>
split_argument(const std::string_view param)
{
    auto position = param.find_first_of(":=");

    if (position == std::string_view::npos) /* If nothing is found, return an
                                               empty tuple.*/
        return {};

    if (position + 1 >= param.size())
        return std::make_tuple(param, std::string_view{});

    return std::make_tuple(param.substr(0, position),
                           param.substr(position + 1));
}

constexpr static bool
starts_with(std::string_view str, std::string_view to_found) noexcept
{
    return !(str.find(to_found));
}

struct get_param
{
    constexpr get_param(int argc_, const char** argv_)
      : argv(argv_)
      , argc(argc_)
      , i(0)
    {}

    constexpr std::optional<std::string_view> operator()(
      int arg_position,
      const std::string_view longp,
      const std::string_view shortp) noexcept
    {
        std::string_view arg(argv[arg_position]);
        i = arg_position;

        if (starts_with(arg, longp) && arg.size() > longp.size() &&
            (arg[longp.size()] == '=' || arg[longp.size()] == ':'))
            return arg.substr(longp.size() + 1);

        if (starts_with(arg, shortp) && arg.size() > shortp.size())
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
      "\nGeneral options:\n"
      "  --help|-h                     This help message\n"
      "  --param|-p [name][:|=][value] Add a new parameter (name is"
      " [a-z][A-Z]_) value can be a double, an integer otherwise a"
      " string.\n"
      "  --disable-preprocessing|-np   Disable preprocessing\n"
      "  --optimize|-O                 Optimize model (default "
      "feasibility search only)\n"
      "  --auto:[manual,nlopt]|-a      Automatic parameters optimization\n"
      "  --check filename.sol          Check if the solution is correct."
      "  --random                      Use the random solver instead of "
      "bastert and wedelin\n"
      "\n"
      "  --quiet                       Remove any verbose message\n"
      "  --verbose|-v int              Set verbose level\n"
      "  --bench|-b file_name.csv      Select the bench mode and store in "
      "file_name.csv\n\n"
      "Parameter list for in the middle heuristic\n"
      " * Global parameters"
      "  - limit: integer ]-oo, +oo[ in loop number\n"
      "  - time-limit: real [0, +oo[ in seconds\n"
      "  - floating-point-type: float double longdouble\n"
      "  - observer-type: none (default), pnm, file\n"
      "  - storage-type: one bound five\n"
      "  - print-level: [0, 2]\n"
      " * In The Middle parameters\n"
      "  - preprocessing: none memory less-greater-equal (or any "
      "combination), p1, p2, p3, p4\n"
      "  - constraint-order: none reversing random-sorting "
      "infeasibility-decr infeasibility-incr lagrangian-decr lagrangian-incr "
      "pi-sign-change\n"
      "  - theta: real [0, 1]\n"
      "  - delta: real [0, +oo[\n"
      "  - kappa-min: real [0, 1[\n"
      "  - kappa-step: real [0, 1[\n"
      "  - kappa-max: real [0, 1[\n"
      "  - alpha: integer [0, 2]\n"
      "  - w: integer [0, +oo[\n"
      "  - norm: l1 l2 loo none random\n"
      " * Pushes system parameters\n"
      "  - pushes-limit: integer [0, +oo[\n"
      "  - pushing-objective-amplifier: real [0, +oo[\n"
      "  - pushing-iteration-limit: integer [0, +oo[\n"
      "  - pushing-k-factor: real [0, +oo[\n"
      " * Initialization parameters\n"
      "  - init-policy: bastert pessimistic-solve optimistic-solve cycle\n"
      "  - init-random: real [0, 1]\n");
}

constexpr static bool
is_equal(std::string_view name, const char* longf, char shortf = '\0')
{
    if (name.compare(0, std::string::npos, longf) == 0)
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
    parameter_missing,
    verbose_error,
    limit_error,
    time_limit_error,
    floating_point_type_error,
    observer_type_error,
    print_level_error,
    preprocessing_error,
    constraint_order_error,
    storage_type_error,
    theta_error,
    delta_error,
    kappa_min_error,
    kappa_step_error,
    kappa_max_error,
    alpha_error,
    w_error,
    norm_error,
    pushes_limit_error,
    pushing_objective_amplifier_error,
    pushing_iteration_limit_error,
    pushing_k_factor_error,
    init_policy_error,
    init_policy_random_error,
    init_random_error,
    thread_error,
    seed_error
};

constexpr const std::string_view command_line_status_string[] = {
    "success",
    "unknown",
    "parameter_missing",
    "verbose_error"
    "limit",
    "time_limit",
    "verbose_error",
    "floating_point_type",
    "observer_type",
    "print_level",
    "preprocessing",
    "constraint_order",
    "storage_type",
    "theta",
    "delta",
    "kappa_min",
    "kappa_step",
    "kappa_max",
    "alpha",
    "w",
    "norm",
    "pushes_limit",
    "pushing_objective_amplifier",
    "pushing_iteration_limit",
    "pushing_k_factor",
    "init_policy",
    "init_policy_random",
    "init_random",
    "thread",
    "seed"
};

template<typename Integer>
constexpr typename std::make_unsigned<Integer>::type
to_unsigned(Integer value)
{
    assert(value >= 0 && "Signed to unsigned error: negative value");

    return static_cast<typename std::make_unsigned<Integer>::type>(value);
}

template<typename Integer>
constexpr typename std::make_signed<Integer>::type
to_signed(Integer value)
{
    assert(static_cast<std::uintmax_t>(value) <
             static_cast<std::uintmax_t>(
               std::numeric_limits<
                 typename std::make_signed<Integer>::type>::max()) &&
           "Unsigned to signed error: too big unsigned");

    return static_cast<typename std::make_signed<Integer>::type>(value);
}

constexpr const std::string_view
to_string(command_line_status s) noexcept
{
    auto x = static_cast<std::underlying_type<command_line_status>::type>(s);

    assert(x < to_signed(std::size(command_line_status_string)));

    return command_line_status_string[x];
}

constexpr static command_line_status
assign_parameter(baryonyx::solver_parameters& params,
                 std::string_view name,
                 std::string_view value)
{
    if (is_equal(name, "limit", 'l')) {
        if (auto v = assign(value, -1, std::numeric_limits<int>::max()); !v)
            return command_line_status::limit_error;
        else
            params.limit = *v;

    } else if (is_equal(name, "time-limit")) {
        if (auto v = to_double(value); !v)
            return command_line_status::time_limit_error;
        else
            params.time_limit = *v <= 0.0 ? -1.0 : *v;

    } else if (is_equal(name, "floating-point-type")) {
        if (value == "float")
            params.float_type =
              baryonyx::solver_parameters::floating_point_type::float_type;
        else if (value == "double")
            params.float_type =
              baryonyx::solver_parameters::floating_point_type::double_type;
        else if (value == "longdouble")
            params.float_type = baryonyx::solver_parameters::
              floating_point_type::longdouble_type;
        else
            return command_line_status::floating_point_type_error;

    } else if (is_equal(name, "observer-type")) {
        if (value == "none")
            params.observer = baryonyx::solver_parameters::observer_type::none;
        else if (value == "pnm")
            params.observer = baryonyx::solver_parameters::observer_type::pnm;
        else if (value == "file")
            params.observer = baryonyx::solver_parameters::observer_type::file;
        else
            return command_line_status::observer_type_error;

    } else if (is_equal(name, "print-level")) {
        if (auto v = assign(value, 0, 2); !v)
            return command_line_status::print_level_error;
        else
            params.print_level = *v;

    } else if (is_equal(name, "preprocessing")) {
        if (value == "none")
            params.pre_order =
              baryonyx::solver_parameters::pre_constraint_order::none;
        else if (value == "memory")
            params.pre_order =
              baryonyx::solver_parameters::pre_constraint_order::memory;
        else if (value == "p1")
            params.pre_order =
              baryonyx::solver_parameters::pre_constraint_order::p1;
        else if (value == "p2")
            params.pre_order =
              baryonyx::solver_parameters::pre_constraint_order::p2;
        else if (value == "p3")
            params.pre_order =
              baryonyx::solver_parameters::pre_constraint_order::p3;
        else if (value == "p4")
            params.pre_order =
              baryonyx::solver_parameters::pre_constraint_order::p4;
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
        else
            return command_line_status::constraint_order_error;

    } else if (is_equal(name, "storage-type")) {
        if (value == "five")
            params.storage = baryonyx::solver_parameters::storage_type::five;
        else if (value == "bound")
            params.storage = baryonyx::solver_parameters::storage_type::bound;
        else if (value == "one")
            params.storage = baryonyx::solver_parameters::storage_type::one;
        else
            return command_line_status::storage_type_error;

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

    } else if (is_equal(name, "kappa-max")) {
        if (auto v = assign_01(value); !v)
            return command_line_status::kappa_max_error;
        else
            params.kappa_max = *v;

    } else if (is_equal(name, "alpha")) {
        if (auto v = assign_d(value, 0, 2); !v)
            return command_line_status::alpha_error;
        else
            params.alpha = *v;

    } else if (is_equal(name, "w")) {
        if (auto v = assign_d(value, 0.0, std::numeric_limits<double>::max());
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

    } else if (is_equal(name, "pushes-limit")) {
        if (auto v = assign(value, 0, std::numeric_limits<int>::max()); !v)
            return command_line_status::pushes_limit_error;
        else
            params.pushes_limit = *v;

    } else if (is_equal(name, "pushing-objective-amplifier")) {
        if (auto v = assign_0oo(value); !v)
            return command_line_status::pushing_objective_amplifier_error;
        else
            params.pushing_objective_amplifier = *v;

    } else if (is_equal(name, "pushing-iteration-limit")) {
        if (auto v = assign(value, 0, std::numeric_limits<int>::max()); !v)
            return command_line_status::pushing_iteration_limit_error;
        else
            params.pushing_iteration_limit = *v;

    } else if (is_equal(name, "pushing-k-factor")) {
        if (auto v = assign_0oo(value); !v)
            return command_line_status::pushing_k_factor_error;
        else
            params.pushing_k_factor = *v;

    } else if (is_equal(name, "init-policy")) {
        if (value == "bastert")
            params.init_policy =
              baryonyx::solver_parameters::init_policy_type::bastert;
        else if (value == "pessimistic-solve")
            params.init_policy =
              baryonyx::solver_parameters::init_policy_type::pessimistic_solve;
        else if (value == "optimistic-solve")
            params.init_policy =
              baryonyx::solver_parameters::init_policy_type::optimistic_solve;
        else if (value == "cycle")
            params.init_policy =
              baryonyx::solver_parameters::init_policy_type::cycle;
        else
            return command_line_status::init_policy_error;

    } else if (is_equal(name, "init-random")) {
        if (auto v = assign_01(value); !v)
            return command_line_status::init_random_error;
        else
            params.init_random = *v;

    } else if (is_equal(name, "init-policy-random")) {
        if (auto v = assign_01(value); !v)
            return command_line_status::init_policy_random_error;
        else
            params.init_policy_random = *v;

    } else if (is_equal(name, "thread")) {
        if (auto v = assign(value, 0, std::numeric_limits<int>::max()); !v)
            return command_line_status::thread_error;
        else
            params.thread = *v;

    } else if (is_equal(name, "seed")) {
        if (auto v = assign(value, 0, std::numeric_limits<int>::max()); !v)
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
            if (*opt == "manual")
                ret.parameters.mode =
                  baryonyx::solver_parameters::mode_type::manual;
            else if (*opt == "nlopt")
                ret.parameters.mode =
                  baryonyx::solver_parameters::mode_type::nlopt;
            else if (*opt == "branch")
                ret.parameters.mode =
                  baryonyx::solver_parameters::mode_type::branch;
            else if (*opt == "branch-manual")
                ret.parameters.mode =
                  baryonyx::solver_parameters::mode_type::manual |
                  baryonyx::solver_parameters::mode_type::branch;
            else if (*opt == "branch-nlopt")
                ret.parameters.mode =
                  baryonyx::solver_parameters::mode_type::nlopt |
                  baryonyx::solver_parameters::mode_type::branch;
            else
                ret.parameters.mode =
                  baryonyx::solver_parameters::mode_type::none;

            i = get.i;
            continue;
        }

        if (auto opt = get(i, "--verbose", "-v"); opt) {
            if (auto v = ::assign(*opt, 0, 6); v) {
                ret.verbose = *v;
            } else {
                ret.parse_status = command_line_status::verbose_error;
                return ret;
            }
        }

        if (auto opt = get(i, "--check", "-C"); opt) {
            ret.check_filename = *opt;
            i = get.i;
            continue;
        }

        if (auto opt = get(i, "--param", "-p"); opt) {
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

        if (arg == "--random") {
            ret.parameters.solver =
              baryonyx::solver_parameters::solver_type::random;
            continue;
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
    }
}

static inline auto
solve_or_optimize(const baryonyx::context_ptr& ctx,
                  const baryonyx::raw_problem& pb,
                  bool with_optimization) -> baryonyx::result
{
    return with_optimization ? baryonyx::optimize(ctx, pb)
                             : baryonyx::solve(ctx, pb);
}

int
main(int argc, const char* argv[])
{
    auto params = parse(argc, argv);

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
    auto ctx = baryonyx::make_context(stdout, params.verbose);
    context_set_solver_parameters(ctx, params.parameters);

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
                fmt::print(stderr,
                           "Fail to read file: {} ({})\n",
                           file_format_error_format(pb.status),
                           static_cast<int>(pb.status));
            } else {
                try {
                    auto filename = fmt::format(
                      "{}-{}.sol", params.filenames.front(), get_pid());
                    fmt::print("  - output file: {}\n", filename);

                    if (params.check) {
                        auto result =
                          baryonyx::make_result(ctx, params.check_filename);

                        if (result) {
                            auto valid =
                              baryonyx::is_valid_solution(pb, result);
                            fmt::print("Check {} with {}: {}",
                                       params.filenames.front(),
                                       params.check_filename,
                                       (valid ? "success" : "failure"));
                        }
                    }

                    resume(pb);
                    std::ofstream ofs(filename);

                    auto now = std::chrono::system_clock::now();
                    auto in_time_t = std::chrono::system_clock::to_time_t(now);
                    resume(pb, ofs);

                    fmt::print(ofs,
                               "\\ solver starts: \n",
                               std::put_time(std::localtime(&in_time_t),
                                             "%Y-%m-%d %X"));

                    auto ret = solve_or_optimize(ctx, pb, params.optimize);
                    in_time_t = std::chrono::system_clock::to_time_t(now);

                    fmt::print(ofs,
                               "\\ solver finishes: {}\n",
                               std::put_time(std::localtime(&in_time_t),
                                             "%Y-%m-%d %X"));

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
                } catch (const baryonyx::precondition_failure& e) {
                    fmt::print(stderr, "internal failure: {}\n", e.what());
                } catch (const baryonyx::postcondition_failure& e) {
                    fmt::print(stderr, "internal failure: {}\n", e.what());
                } catch (const baryonyx::numeric_cast_failure& e) {
                    fmt::print(
                      stderr, "numeric cast internal failure: {}\n", e.what());
                } catch (const baryonyx::file_access_failure& e) {
                    fmt::print(stderr,
                               "file `{}' fail {}: {}\n",
                               e.file(),
                               e.error(),
                               std::strerror(e.error()));
                } catch (const baryonyx::file_format_failure& e) {
                    fmt::print(stderr,
                               "file format error at line {} column {} "
                               "{}\n",
                               e.line(),
                               e.column(),
                               file_format_error_format(e.failure()));
                } catch (const baryonyx::problem_definition_failure& e) {
                    fmt::print(stderr,
                               "definition problem error at {}: {}\n",
                               e.element(),
                               problem_definition_error_format(e.failure()));
                } catch (const baryonyx::solver_failure& e) {
                    fmt::print(stderr,
                               "solver error: {}\n",
                               solver_error_format(e.failure()));
                } catch (const std::exception& e) {
                    fmt::print(stderr, "failure: {}.\n", e.what());
                }
            }
        } else {
            auto filename = fmt::format("baryonyx-{}.res", get_pid());
            std::ofstream ofs(filename);

            for (auto& elem : params.filenames) {
                auto pb = baryonyx::make_problem(ctx, elem);
                if (!pb) {
                    fmt::print(stderr,
                               "Fail to read file: {} ({})\n",
                               file_format_error_format(pb.status),
                               static_cast<int>(pb.status));
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
                } catch (const baryonyx::precondition_failure& e) {
                    fmt::print(stderr, "internal failure: {}\n", e.what());
                } catch (const baryonyx::postcondition_failure& e) {
                    fmt::print(stderr, "internal failure: {}\n", e.what());
                } catch (const baryonyx::numeric_cast_failure& e) {
                    fmt::print(
                      stderr, "numeric cast internal failure: {}\n", e.what());
                } catch (const baryonyx::file_access_failure& e) {
                    fmt::print(stderr,
                               "file `{}' fail {}: {}\n",
                               e.file(),
                               e.error(),
                               std::strerror(e.error()));
                } catch (const baryonyx::file_format_failure& e) {
                    fmt::print(stderr,
                               "file format error at line {} column {} "
                               "{}\n",
                               e.line(),
                               e.column(),
                               file_format_error_format(e.failure()));
                } catch (const baryonyx::problem_definition_failure& e) {
                    fmt::print(stderr,
                               "definition problem error at {}: {}\n",
                               e.element(),
                               problem_definition_error_format(e.failure()));
                } catch (const baryonyx::solver_failure& e) {
                    fmt::print(stderr,
                               "solver error: {}\n",
                               solver_error_format(e.failure()));
                } catch (const std::exception& e) {
                    fmt::print(stderr, "failure: {}.\n", e.what());
                }
            }
        }
    }

    return EXIT_SUCCESS;
}
