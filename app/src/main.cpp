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

#include <fmt/ostream.h>
#include <fmt/printf.h>

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

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

double
to_double(const char* s, double bad_value) noexcept
{
    char* c;
    errno = 0;
    double value = std::strtod(s, &c);

    if ((errno == ERANGE and (value == HUGE_VAL or value == -HUGE_VAL)) or
        (value == 0.0 and c == s))
        return bad_value;

    return value;
}

int
to_int(const char* s, int bad_value) noexcept
{
    char* c;
    errno = 0;
    long value = std::strtol(s, &c, 10);

    if ((errno == ERANGE and (value == LONG_MIN or value == LONG_MAX)) or
        (value == 0 and c == s))
        return bad_value;

    if (value < INT_MIN)
        return INT_MIN;

    if (value > INT_MAX)
        return INT_MAX;

    return static_cast<int>(value);
}

std::tuple<std::string, baryonyx::parameter>
split_param(const char* param) noexcept
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

    auto valuel = to_int(value.c_str(), INT_MIN);
    auto valued = to_double(value.c_str(), -HUGE_VAL);

    double tmp;
    if (valued != -HUGE_VAL and std::modf(valued, &tmp))
        return std::make_tuple(name, baryonyx::parameter(valued));

    if (valuel != INT_MIN)
        return std::make_tuple(name, baryonyx::parameter(valuel));

    return std::make_tuple(name, baryonyx::parameter(value));
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
            argv[arg][longp_length] == '=') {
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

struct main_parameters
{
    std::unordered_map<std::string, baryonyx::parameter> params;
    std::vector<std::string> filenames;
    std::string check_filename;
    int verbose = 6;
    int limit = 1000;
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
            ret.params["limit"] = ::to_int(opt, 1000);
            i = get.i;
            continue;
        }

        if ((opt = get(i, "--verbose", "-v"))) {
            ret.verbose = ::to_int(opt, 3);
            i = get.i;
            continue;
        }

        if ((opt = get(i, "--check", "-C"))) {
            ret.check_filename = opt;
            i = get.i;
            continue;
        }

        if ((opt = get(i, "--param", "-p"))) {
            std::string name;
            baryonyx::parameter value;
            std::tie(name, value) = ::split_param(opt);
            ret.params[name] = value;
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

    auto ctx =
      baryonyx::make_context(stdout, params.quiet ? 3 : params.verbose);

    context_set_parameters(ctx, std::move(params.params));

    if (params.filenames.size() == 1) {
        try {
            auto pb = baryonyx::make_problem(ctx, params.filenames.front());
            fmt::print("{}", baryonyx::resume(pb, false));

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
                fmt::print(
                  "Best solution found: {} in {}s\n", ret.value, ret.duration);
            } else {
                fmt::print("No solution found. Missing constraints: {}\n",
                           ret.remaining_constraints);
            }

            in_time_t = std::chrono::system_clock::to_time_t(now);
            ofs << R"(\ solver finishes: )"
                << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X")
                << '\n';

            if (ret.status == baryonyx::result_status::success) {
                ofs << R"(\ Solution found: )" << ret.value << '\n' << ret;
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
                fmt::print("{}", baryonyx::resume(pb, false));
                ofs << elem << " ";

                auto ret = params.optimize ? baryonyx::optimize(ctx, pb)
                                           : baryonyx::solve(ctx, pb);

                if (ret.status == baryonyx::result_status::success) {
                    ofs << ret.value << " " << ret.duration << "\n";
                } else {
                    ofs << "No solution found.\n";
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
