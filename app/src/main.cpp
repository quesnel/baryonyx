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

#include <baryonyx/core-out>
#include <baryonyx/core>

#include <fmt/core.h>
#include <fmt/ostream.h>

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>

#include <cmath>
#include <cstring>
#include <utility>

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

static double
to_double(std::string s, double bad_value) noexcept
{
    char* c;
    errno = 0;
    double value = std::strtod(s.c_str(), &c);

    if ((errno == ERANGE and (value == HUGE_VAL or value == -HUGE_VAL)) or
        (value == 0.0 and c == s.c_str()))
        return bad_value;

    return value;
}

static int
to_int(std::string s, int bad_value) noexcept
{
    char* c;
    errno = 0;
    long value = std::strtol(s.c_str(), &c, 10);

    if ((errno == ERANGE and (value == LONG_MIN or value == LONG_MAX)) or
        (value == 0 and c == s.c_str()))
        return bad_value;

    if (value < INT_MIN)
        return INT_MIN;

    if (value > INT_MAX)
        return INT_MAX;

    return static_cast<int>(value);
}

static std::tuple<std::string, std::string>
split_argument(const char* param)
{
    std::string name, value;

    while (*param) {
        if (isalpha(*param) or *param == '_' or *param == '-')
            name += *param;
        else
            break;

        param++;
    }

    if (*param and (*param == ':' or *param == '=')) {
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
        if (not strncmp(argv[arg], longp, longp_length) and
            strlen(argv[arg]) + 1 > longp_length and
            (argv[arg][longp_length] == '=' or
             argv[arg][longp_length] == ':')) {
            return argv[arg] + longp_length + 1;
        }

        auto shortp_length = strlen(shortp);
        if (not strncmp(argv[arg], shortp, shortp_length) and
            strlen(argv[arg]) > shortp_length) {
            return argv[arg] + shortp_length;
        }

        if (arg + 1 < argc) {
            if (longp and strcmp(argv[arg], longp) == 0) {
                i = arg + 1;
                return argv[i];
            }

            if (shortp and strcmp(argv[i], shortp) == 0) {
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

static const char* file_format_error_format(
  baryonyx::file_format_error_tag) noexcept;

static const char* problem_definition_error_format(
  baryonyx::problem_definition_error_tag) noexcept;

static const char* solver_error_format(baryonyx::solver_error_tag) noexcept;

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
      "  --verbose|-v int              Set verbose level\n\n"
      "Parameter list for in the middle heuristic\n"
      " * Global parameters"
      "  - limit: integer ]-oo, +oo[ in loop number\n"
      "  - time-limit: real [0, +oo[ in seconds\n"
      "  - floating-point-type: float double longdouble\n"
      "  - print-level: [0, 2]\n"
      " * In The Middle parameters\n"
      "  - preprocessing: none variables-number variables-weight "
      "constraints-weight implied\n"
      "  - constraint-order: none reversing random-sorting "
      "infeasibility-decr infeasibility-incr\n"
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
      "  - init-policy: bastert random best\n"
      "  - init-random: real [0, 1]\n");
}

static bool
is_equal(std::string name, const char* longf, char shortf = '\0')
{
    if (name.compare(0, std::string::npos, longf) == 0)
        return true;

    if (shortf != '\0' and name.size() == 2 and name[1] == shortf)
        return true;

    return false;
}

static double
assign_0oo(std::string value, double def)
{
    auto ret = ::to_double(value, def);

    return ret < 0 ? 0 : not std::isnormal(ret) ? 1 : ret;
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

void
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
    } else if (is_equal(name, "print-level")) {
        params.print_level = assign(value, 0, 2, params.print_level);
    } else if (is_equal(name, "preprocessing")) {
        if (value == "none")
            params.pre_order =
              baryonyx::solver_parameters::pre_constraint_order::none;
        else if (value == "variables-number")
            params.pre_order = baryonyx::solver_parameters::
              pre_constraint_order::variables_number;
        else if (value == "variables-weight")
            params.pre_order = baryonyx::solver_parameters::
              pre_constraint_order::variables_weight;
        else if (value == "constraints-weight")
            params.pre_order = baryonyx::solver_parameters::
              pre_constraint_order::constraints_weight;
        else if (value == "implied")
            params.pre_order =
              baryonyx::solver_parameters::pre_constraint_order::implied;
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
        else if (value == "random_sorting")
            params.order =
              baryonyx::solver_parameters::constraint_order::random_sorting;
        else if (value == "infeasibility_decr")
            params.order = baryonyx::solver_parameters::constraint_order::
              infeasibility_decr;
        else if (value == "infeasibility_incr")
            params.order = baryonyx::solver_parameters::constraint_order::
              infeasibility_incr;
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
    int verbose = 6;
    bool check = false;
    bool optimize = false;
    bool quiet = false;
    bool preprocessing = true;
};

main_parameters
parse(int argc, const char* argv[])
{
    main_parameters ret;
    get_param get(argc, argv);

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        const char* opt;

        if (arg == "--help" or arg == "-h") {
            ::help();
            continue;
        }

        if (arg == "--quiet" or arg == "-q") {
            ret.quiet = true;
            continue;
        }

        if (arg == "--optimize" or arg == "-O") {
            ret.optimize = true;
            continue;
        }

        if (arg == "--disable-preprocessing" or arg == "-np") {
            ret.preprocessing = false;
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
                ret.parameters.auto_tune =
                  baryonyx::solver_parameters::auto_tune_parameters::manual;
            else if (str == "nlopt")
                ret.parameters.auto_tune =
                  baryonyx::solver_parameters::auto_tune_parameters::nlopt;

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

int
main(int argc, const char* argv[])
{
    auto params = parse(argc, argv);

    if (params.filenames.empty()) {
        fmt::print(stderr, "Missing lp file, see --help\n");
        return EXIT_SUCCESS;
    }

    params.verbose = params.quiet ? 3 : params.verbose;
    auto ctx = baryonyx::make_context(stdout, params.verbose);
    context_set_solver_parameters(ctx, params.parameters);

    if (params.filenames.size() == 1) {
        try {
            auto pb = baryonyx::make_problem(ctx, params.filenames.front());

            auto filename =
              fmt::format("{}-{}.sol", params.filenames.front(), get_pid());
            fmt::print("Output file: {}\n", filename);

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

            if (params.preprocessing)
                pb = baryonyx::preprocess(ctx, pb);
            fmt::print("{}", baryonyx::resume(pb, false));

            std::ofstream ofs(filename);
            ofs << std::boolalpha
                << std::setprecision(static_cast<int>(std::floor(
                     std::numeric_limits<double>::digits * std::log10(2) +
                     2)));

            auto now = std::chrono::system_clock::now();
            auto in_time_t = std::chrono::system_clock::to_time_t(now);

            ofs << baryonyx::resume(pb) << R"(\ solver starts: )"
                << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X")
                << "\n";

            auto ret = params.optimize ? baryonyx::optimize(ctx, pb)
                                       : baryonyx::solve(ctx, pb);

            if (ret.status == baryonyx::result_status::success) {
                fmt::print("Best solution found: {} in {}s\n",
                           ret.solutions.back().value,
                           ret.duration);
            } else {
                fmt::print("No solution found. Missing constraints: {}\n",
                           ret.remaining_constraints);
            }

            in_time_t = std::chrono::system_clock::to_time_t(now);
            ofs << R"(\ solver finishes: )"
                << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X")
                << '\n';

            if (ret.status == baryonyx::result_status::success) {
                ofs << R"(\ Solution found: )" << ret.solutions.back().value
                    << '\n'
                    << ret;
            } else {
                ofs << R"(\ Solution not found. Missing constraints: )"
                    << ret.remaining_constraints << '\n';
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
                       "%s\n",
                       e.line(),
                       e.column(),
                       file_format_error_format(e.failure()));
        } catch (const baryonyx::problem_definition_failure& e) {
            fmt::print(stderr,
                       "definition problem error at {}: {}\n",
                       e.element(),
                       problem_definition_error_format(e.failure()));
        } catch (const baryonyx::solver_failure& e) {
            fmt::print(
              stderr, "solver error: {}\n", solver_error_format(e.failure()));
        } catch (const std::exception& e) {
            fmt::print(stderr, "failure: {}.\n", e.what());
        }
    } else {
        auto filename = fmt::format("baryonyx-{}.res", get_pid());
        std::ofstream ofs(filename);
        ofs << std::boolalpha
            << std::setprecision(static_cast<int>(std::floor(
                 std::numeric_limits<double>::digits * std::log10(2) + 2)));

        for (auto& elem : params.filenames) {
            try {
                auto pb = baryonyx::make_problem(ctx, elem);

                if (params.preprocessing)
                    pb = baryonyx::preprocess(ctx, pb);

                fmt::print("{}", baryonyx::resume(pb, false));
                ofs << elem << " ";

                auto ret = params.optimize ? baryonyx::optimize(ctx, pb)
                                           : baryonyx::solve(ctx, pb);

                if (ret.status == baryonyx::result_status::success) {
                    ofs << ret.solutions.back().value << " " << ret.duration
                        << ' ';

                    for (const auto& sol : ret.solutions)
                        ofs << sol.value << ' ';
                    ofs << "\n";

                } else {
                    ofs << "No solution found.\n";
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
                           "%s\n",
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
    return EXIT_SUCCESS;
}

static const char*
file_format_error_format(baryonyx::file_format_error_tag failure) noexcept
{
    static const char* const tag[] = {
        "end of file",     "unknown",
        "already defined", "incomplete",
        "bad name",        "bad operator",
        "bad integer",     "bad objective function type",
        "bad bound",       "bad function element",
        "bad constraint"
    };

    return tag[static_cast<int>(failure)];
}

static const char*
problem_definition_error_format(
  baryonyx::problem_definition_error_tag failure) noexcept
{
    static const char* const tag[] = {
        "empty variables",
        "empty objective function",
        "variable not used",
        "bad bound",
        "multiple constraints with different value"
    };

    return tag[static_cast<int>(failure)];
}

static const char*
solver_error_format(baryonyx::solver_error_tag failure) noexcept
{
    static const char* const tag[] = { "no solver available",
                                       "unrealisable constraint",
                                       "not enough memory" };

    return tag[static_cast<int>(failure)];
}
