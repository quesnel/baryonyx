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

struct write_parameters
{
    std::shared_ptr<baryonyx::context> ctx;

    write_parameters(std::shared_ptr<baryonyx::context> ctx_)
      : ctx(std::move(ctx_))
    {
    }
};

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

static const char* file_format_error_format(
  baryonyx::file_format_error_tag) noexcept;

static const char* problem_definition_error_format(
  baryonyx::problem_definition_error_tag) noexcept;

static const char* solver_error_format(baryonyx::solver_error_tag) noexcept;

static std::ostream&
operator<<(std::ostream& os, const write_parameters& wp);

static baryonyx::result
solve_or_optimize(const std::shared_ptr<baryonyx::context>& ctx,
                  baryonyx::problem& pb);

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

    fmt::print("\nGeneral options:\n"
               "  --help|-h                   This help message\n"
               "  --param|-p [name][:|=][value]   Add a new parameter (name is"
               " [a-z][A-Z]_) value can be a double, an integer otherwise a"
               " string.\n"
               "  --optimize|-O               Optimize model (default "
               "feasibility search only)\n"
               "  --check filename.sol        Check if the solution is correct."
               "\n"
               "  --quiet                     Remove any verbose message\n"
               "  --verbose|-v int            Set verbose level\n\n"
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

static std::tuple<std::unordered_map<std::string, baryonyx::parameter>,
                  std::vector<std::string>>
parse(int argc, char* argv[])
{
    std::vector<std::string> files;
    std::unordered_map<std::string, baryonyx::parameter> params;
    params["verbose"] = baryonyx::parameter(6);
    params["limit"] = baryonyx::parameter(1000);

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);

        if (arg == "--help" or arg == "-h") {
            ::help();
            continue;
        }

        if (arg == "--quiet" or arg == "-q") {
            params["quiet"] = 1;
            continue;
        }

        if (arg == "--optimize" or arg == "-O") {
            params["optimize"] = 1;
            continue;
        }

        if (arg == "--limit" or arg == "-l") {
            if (i + 1 >= argc) {
                fmt::print(stderr, "/!\\ --limit argument required");
            } else {
                params["limit"] = to_int(argv[i + 1], 1000);
                ++i;
            }
            continue;
        }

        if (arg == "--verbose" or arg == "-v") {
            if (i + 1 >= argc) {
                fmt::print(stderr, "/!\\ --verbose argument required\n");
            } else {
                params["verbose"] = ::to_int(argv[i + 1], 3);
                ++i;
            }
            continue;
        }

        if (arg == "--check" or arg == "-C") {
            if (i + 1 >= argc) {
                fmt::print(stderr, "/!\\ --check argument required\n");
            } else {
                params["check-filename"] = argv[i + 1];
                ++i;
            }
            continue;
        }

        if (arg == "--param" or arg == "-p") {
            if (i + 1 >= argc) {
                fmt::print(stderr, "/!\\ --param argument required\n");
            } else {
                std::string name;
                baryonyx::parameter value;
                std::tie(name, value) = ::split_param(argv[i + 1]);
                params[name] = value;
                ++i;
            }
            continue;
        }

        files.emplace_back(argv[i]);
    }

    return std::make_tuple(params, files);
}

int
main(int argc, char* argv[])
{
    std::vector<std::string> files;
    std::unordered_map<std::string, baryonyx::parameter> params;

    std::tie(params, files) = parse(argc, argv);

    if (params.empty()) {
        fmt::print(stderr, "Argument error, see --help\n");
        return EXIT_FAILURE;
    }

    if (files.empty()) {
        fmt::print(stderr, "Missing lp file, see --help\n");
        return EXIT_SUCCESS;
    }

    auto ctx = std::make_shared<baryonyx::context>();
    ctx->set_parameters(std::move(params));

    if (files.size() == 1) {
        try {
            auto pb = baryonyx::make_problem(ctx, files.front());
            fmt::print("{}", baryonyx::resume(pb, false));

            auto filename = fmt::format("{}-{}.sol", files.front(), get_pid());
            fmt::print("Output file: {}\n", filename);

            if (ctx->check()) {
                auto filename =
                  ctx->get_string_parameter("check-filename", std::string());

                auto result = baryonyx::make_result(ctx, filename);
                if (result) {
                    auto valid = baryonyx::is_valid_solution(pb, result);
                    fmt::print("Check {} with {}: {}",
                               files.front(),
                               filename,
                               (valid ? "success" : "failure"));
                }
            }

            std::ofstream ofs(filename);
            ofs << std::boolalpha
                << std::setprecision(static_cast<int>(std::floor(
                     std::numeric_limits<double>::digits * std::log10(2) + 2)));

            auto now = std::chrono::system_clock::now();
            auto in_time_t = std::chrono::system_clock::to_time_t(now);

            ofs << baryonyx::resume(pb) << R"(\ solver starts: )"
                << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X")
                << "\n\\ parameters:\n"
                << ::write_parameters(ctx);

            auto ret = solve_or_optimize(ctx, pb);
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
            fmt::print(stderr, "internal failure\n");
        } catch (const baryonyx::postcondition_failure& e) {
            fmt::print(stderr, "internal failure\n");
        } catch (const baryonyx::numeric_cast_failure& e) {
            fmt::print(stderr, "numeric cast internal failure\n");
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

        for (auto& elem : files) {
            try {
                auto pb = baryonyx::make_problem(ctx, elem);
                fmt::print("{}", baryonyx::resume(pb, false));
                ofs << elem << " ";

                auto ret = solve_or_optimize(ctx, pb);
                if (ret.status == baryonyx::result_status::success) {
                    ofs << ret.value << " " << ret.duration << "\n";
                } else {
                    ofs << "No solution found.\n";
                }
            } catch (const baryonyx::precondition_failure& e) {
                fmt::print(stderr, "internal failure\n");
            } catch (const baryonyx::postcondition_failure& e) {
                fmt::print(stderr, "internal failure\n");
            } catch (const baryonyx::numeric_cast_failure& e) {
                fmt::print(stderr, "numeric cast internal failure\n");
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

static std::ostream&
operator<<(std::ostream& os, const write_parameters& wp)
{
    if (not wp.ctx)
        return os;

    const auto& params = wp.ctx->get_parameters();
    for (const auto& param : params) {
        os << R"(\ )" << param.first << " = ";

        switch (param.second.type) {
        case baryonyx::parameter::tag::string:
            os << param.second.s;
            break;
        case baryonyx::parameter::tag::integer:
            os << param.second.l;
            break;
        case baryonyx::parameter::tag::real:
            os << param.second.d;
            break;
        }

        os << '\n';
    }

    return os;
}

static baryonyx::result
solve_or_optimize(const std::shared_ptr<baryonyx::context>& ctx,
                  baryonyx::problem& pb)
{
    if (ctx->optimize())
        return baryonyx::optimize(ctx, pb);

    return baryonyx::solve(ctx, pb);
}
