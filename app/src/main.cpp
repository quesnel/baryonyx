/* Copyright (C) 2016-2018 INRA
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

static std::tuple<std::string, std::string>
split_argument(const char* param)
{
    std::string name, value;

    while (*param) {
        if (isalpha(*param) || *param == '_' || *param == '-')
            name += *param;
        else
            break;

        param++;
    }

    if (*param && (*param == ':' || *param == '=')) {
        param++;

        while (*param)
            value += *param++;
    }

    return std::make_tuple(name, value);
}

struct get_param
{
    get_param(int argc_, const char** argv_)
      : argv(argv_)
      , argc(argc_)
      , i(0)
    {}

    const char* operator()(int arg,
                           const char* longp,
                           const char* shortp) noexcept
    {
        i = arg;

        auto longp_length = strlen(longp);
        if (!strncmp(argv[arg], longp, longp_length) &&
            strlen(argv[arg]) + 1 > longp_length &&
            (argv[arg][longp_length] == '=' ||
             argv[arg][longp_length] == ':')) {
            return argv[arg] + longp_length + 1;
        }

        auto shortp_length = strlen(shortp);
        if (!strncmp(argv[arg], shortp, shortp_length) &&
            strlen(argv[arg]) > shortp_length) {
            return argv[arg] + shortp_length;
        }

        if (arg + 1 < argc) {
            if (longp && strcmp(argv[arg], longp) == 0) {
                i = arg + 1;
                return argv[i];
            }

            if (shortp && strcmp(argv[i], shortp) == 0) {
                i = arg + 1;
                return argv[i];
            }
        }

        return nullptr;
    }

    const char** argv;
    int argc;
    int i;
};

static void
help() noexcept
{
    fmt::print(stdout,
               "Baryonyx v{}.{}.{}",
               VERSION_MAJOR,
               VERSION_MINOR,
               VERSION_PATCH);

    if (VERSION_TWEAK)
        fmt::print("-{}", VERSION_TWEAK);

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
      "  - preprocessing: none memory less,greater,equal (or any "
      "combination), p1, p2, p3, p4\n"
      "  - constraint-order: none reversing random-sorting "
      "infeasibility-decr infeasibility-incr lagrangian-decr lagrangian-incr\n"
      "  - theta: real [0, 1]\n"
      "  - delta: real [0, +oo[\n"
      "  - kappa-min: real [0, 1[\n"
      "  - kappa-step: real [0, 1[\n"
      "  - kappa-max: real [0, 1[\n"
      "  - alpha: integer [0, 2]\n"
      "  - w: integer [0, +oo[\n"
      "  - norm: l1 l2 inf none rng\n"
      " * Pushes system parameters\n"
      "  - pushes-limit: integer [0, +oo[\n"
      "  - pushing-objective-amplifier: real [0, +oo[\n"
      "  - pushing-iteration-limit: integer [0, +oo[\n"
      "  - pushing-k-factor: real [0, +oo[\n"
      " * Initialization parameters\n"
      "  - init-policy: bastert random best bastert-cycle random-cycle "
      "best-cycle\n"
      "  - init-random: real [0, 1]\n");
}

static bool
is_equal(std::string name, const char* longf, char shortf = '\0')
{
    if (name.compare(0, std::string::npos, longf) == 0)
        return true;

    if (shortf != '\0' && name.size() == 2 && name[1] == shortf)
        return true;

    return false;
}

static double
assign_0oo(std::string value, double def)
{
    auto ret = ::to_double(value, def);

    return ret < 0 ? 0 : !std::isnormal(ret) ? 1 : ret;
}

static double
assign_01(std::string value, double def)
{
    auto ret = ::to_double(value, def);

    return ret < 0 ? 0 : ret > 1 ? 1 : ret;
}

static int
assign(std::string value, int mindef, int maxdef, int def)
{
    auto ret = ::to_int(value, def);

    return ret < mindef ? mindef : ret > maxdef ? maxdef : ret;
}

static double
assign_d(std::string value, double mindef, double maxdef, double def)
{
    auto ret = ::to_double(value, def);

    return ret < mindef ? mindef : ret > maxdef ? maxdef : ret;
}

static void
assign_parameter(baryonyx::solver_parameters& params,
                 std::string name,
                 std::string value)
{
    if (is_equal(name, "limit", 'l')) {
        params.limit =
          assign(value, -1, std::numeric_limits<int>::max(), params.limit);
    } else if (is_equal(name, "time-limit")) {
        params.time_limit = ::to_double(value, 60);
        if (params.time_limit <= 0)
            params.time_limit = -1;
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
    } else if (is_equal(name, "observer-type")) {
        if (value == "none")
            params.observer = baryonyx::solver_parameters::observer_type::none;
        else if (value == "pnm")
            params.observer = baryonyx::solver_parameters::observer_type::pnm;
        else if (value == "file")
            params.observer = baryonyx::solver_parameters::observer_type::file;
    } else if (is_equal(name, "print-level")) {
        params.print_level = assign(value, 0, 2, params.print_level);
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
    } else if (is_equal(name, "storage-type")) {
        if (value == "five")
            params.storage = baryonyx::solver_parameters::storage_type::five;
        else if (value == "bound")
            params.storage = baryonyx::solver_parameters::storage_type::bound;
        else
            params.storage = baryonyx::solver_parameters::storage_type::one;
    } else if (is_equal(name, "theta")) {
        params.theta = assign_01(value, params.theta);
    } else if (is_equal(name, "delta")) {
        params.delta = assign_0oo(value, params.delta);
    } else if (is_equal(name, "kappa-min")) {
        params.kappa_min = assign_01(value, params.kappa_min);
    } else if (is_equal(name, "kappa-step")) {
        params.kappa_step = assign_01(value, params.kappa_step);
    } else if (is_equal(name, "kappa-max")) {
        params.kappa_max = assign_01(value, params.kappa_max);
    } else if (is_equal(name, "alpha")) {
        params.alpha = assign_d(value, 0, 2, params.alpha);
    } else if (is_equal(name, "w")) {
        params.w = assign(value, 0, std::numeric_limits<int>::max(), params.w);
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
    } else if (is_equal(name, "pushes-limit")) {
        params.pushes_limit = assign(
          value, 0, std::numeric_limits<int>::max(), params.pushes_limit);
    } else if (is_equal(name, "pushing-objective-amplifier")) {
        params.pushing_objective_amplifier =
          assign_0oo(value, params.pushing_objective_amplifier);
    } else if (is_equal(name, "pushing-iteration-limit")) {
        params.pushing_iteration_limit =
          assign(value,
                 0,
                 std::numeric_limits<int>::max(),
                 params.pushing_iteration_limit);
    } else if (is_equal(name, "pushing-k-factor")) {
        params.pushing_k_factor = assign_0oo(value, params.pushing_k_factor);
    } else if (is_equal(name, "init-policy")) {
        if (value == "bastert")
            params.init_policy =
              baryonyx::solver_parameters::init_policy_type::bastert;
        else if (value == "random")
            params.init_policy =
              baryonyx::solver_parameters::init_policy_type::random;
        else if (value == "best")
            params.init_policy =
              baryonyx::solver_parameters::init_policy_type::best;
        else if (value == "bastert-cycle")
            params.init_policy =
              baryonyx::solver_parameters::init_policy_type::bastert_cycle;
        else if (value == "random-cycle")
            params.init_policy =
              baryonyx::solver_parameters::init_policy_type::random_cycle;
        else if (value == "best-cycle")
            params.init_policy =
              baryonyx::solver_parameters::init_policy_type::best_cycle;
    } else if (is_equal(name, "init-random")) {
        params.init_random = assign_01(value, params.init_random);
    } else if (is_equal(name, "thread")) {
        params.thread = assign(value, 0, std::numeric_limits<int>::max(), 0);
    } else if (is_equal(name, "seed")) {
        params.seed = assign(value, 0, std::numeric_limits<int>::max(), 0);
    } else
        fmt::print(
          stderr, "Unknown parameters {} = {}, ignoring.\n", name, value);
}

struct main_parameters
{
    baryonyx::solver_parameters parameters;
    std::vector<std::string> filenames;
    std::string check_filename;
    std::string bench_name;
    int verbose = 6;
    bool check = false;
    bool optimize = false;
    bool quiet = false;
    bool bench = false;
};

static main_parameters
parse(int argc, const char* argv[])
{
    main_parameters ret;
    get_param get(argc, argv);

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        const char* opt;

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

        if ((opt = get(i, "--bench", "-b"))) {
            ret.bench = true;
            ret.bench_name = opt;
            i = get.i;
            continue;
        }

        if (arg == "--disable-preprocessing" || arg == "-np") {
            ret.parameters.preprocessor =
              baryonyx::solver_parameters::preprocessor_options::none;
            continue;
        }

        if ((opt = get(i, "--limit", "-l"))) {
            ret.parameters.limit = ::to_int(opt, ret.parameters.limit);
            i = get.i;
            continue;
        }

        if ((opt = get(i, "--auto", "-a"))) {
            std::string str(opt);

            if (str == "manual")
                ret.parameters.mode =
                  baryonyx::solver_parameters::mode_type::manual;
            else if (str == "nlopt")
                ret.parameters.mode =
                  baryonyx::solver_parameters::mode_type::nlopt;
            else if (str == "branch")
                ret.parameters.mode =
                  baryonyx::solver_parameters::mode_type::branch;
            else if (str == "branch-manual")
                ret.parameters.mode =
                  baryonyx::solver_parameters::mode_type::manual |
                  baryonyx::solver_parameters::mode_type::branch;
            else if (str == "branch-nlopt")
                ret.parameters.mode =
                  baryonyx::solver_parameters::mode_type::nlopt |
                  baryonyx::solver_parameters::mode_type::branch;
            else
                ret.parameters.mode =
                  baryonyx::solver_parameters::mode_type::none;

            i = get.i;
            continue;
        }

        if ((opt = get(i, "--verbose", "-v"))) {
            ret.verbose = ::to_int(opt, 6);
            i = get.i;
            continue;
        }

        if ((opt = get(i, "--check", "-C"))) {
            ret.check_filename = opt;
            i = get.i;
            continue;
        }

        if ((opt = get(i, "--param", "-p"))) {
            std::string name, value;
            std::tie(name, value) = split_argument(opt);
            assign_parameter(ret.parameters, name, value);
            i = get.i;
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

    fmt::print("  - objective : {}\n"
               "  - variables: {}/{}/{} (real/general/binary)\n"
               "  - constraints: {}/{}/{} (equal/greater/less)\n",
               (pb.type == baryonyx::objective_function_type::maximize
                  ? "maximize"
                  : "minimize"),
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
               "\\ variables: {}/{}/{} (real/general/binary)\n"
               "\\ constraints: {}/{}/{} (equal/greater/less)\n",
               (pb.type == baryonyx::objective_function_type::maximize
                  ? "maximize"
                  : "minimize"),
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

    if (params.filenames.empty()) {
        fmt::print(fmt::color::red,
                   "Missing lp/benchmark file(s), see --help\n");
        return EXIT_SUCCESS;
    }

    params.verbose = params.quiet ? 3 : params.verbose;
    auto ctx = baryonyx::make_context(stdout, params.verbose);
    context_set_solver_parameters(ctx, params.parameters);

    if (params.bench) {
        for (auto filename : params.filenames) {
            if (params.bench_name.empty())
                params.bench_name = fmt::format(
                  "bx-{}.{}.{}", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);

            fmt::print("Start benchmark for file {} (name: {})\n",
                       filename,
                       params.bench_name);

            if (!benchmark(ctx, filename, params.bench_name))
                fmt::print(
                  fmt::color::red, "Benchmark for {} failed.\n", filename);
        }
    } else {
        if (params.filenames.size() == 1) {
            try {
                auto pb =
                  baryonyx::make_problem(ctx, params.filenames.front());

                auto filename = fmt::format(
                  "{}-{}.sol", params.filenames.front(), get_pid());
                fmt::print("  - output file: {}\n", filename);

                if (params.check) {
                    auto result =
                      baryonyx::make_result(ctx, params.check_filename);

                    if (result) {
                        auto valid = baryonyx::is_valid_solution(pb, result);
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

                fmt::print(
                  ofs,
                  "\\ solver starts: \n",
                  std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X"));

                auto ret = solve_or_optimize(ctx, pb, params.optimize);

                if (ret.status == baryonyx::result_status::success) {
                    fmt::print("Best solution found: {} in {}s\n",
                               ret.solutions.back().value,
                               ret.duration);
                } else {
                    fmt::print("No solution found. Missing constraints: {}\n",
                               ret.remaining_constraints);
                }

                in_time_t = std::chrono::system_clock::to_time_t(now);

                fmt::print(
                  ofs,
                  "\\ solver finishes: {}\n",
                  std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X"));

                if (ret.status == baryonyx::result_status::success) {
                    fmt::print(ofs,
                               "\\ Solution found: {:f}\n",
                               ret.solutions.back().value);
                    resume(ret, ofs);
                } else {
                    fmt::print(
                      ofs,
                      "\\ Solution not found. Missing constraints: {}\n",
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
        } else {
            auto filename = fmt::format("baryonyx-{}.res", get_pid());
            std::ofstream ofs(filename);

            for (auto& elem : params.filenames) {
                try {
                    auto pb = baryonyx::make_problem(ctx, elem);

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
