/* Copyright (C) 2016 INRA
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

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <lpcore-out>
#include <lpcore>
#include <sstream>

#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <getopt.h>
#include <sys/types.h>
#include <unistd.h>

void help() noexcept;
double to_double(const char* s, double bad_value) noexcept;
long to_long(const char* s, long bad_value) noexcept;
std::tuple<std::string, lp::parameter> split_param(const char* param) noexcept;
const char* file_format_error_format(lp::file_format_error::tag) noexcept;
const char* problem_definition_error_format(
  lp::problem_definition_error::tag) noexcept;
const char* solver_error_format(lp::solver_error::tag) noexcept;

lp::result solve(std::shared_ptr<lp::context> ctx,
                 lp::problem& pb,
                 const std::map<std::string, lp::parameter>& params,
                 bool optimize);

int
main(int argc, char* argv[])
{
    const char* const short_opts = "Ohp:l:qv:";
    const struct option long_opts[] = { { "optimize", 0, nullptr, 'O' },
                                        { "help", 0, nullptr, 'h' },
                                        { "param", 1, nullptr, 'p' },
                                        { "limit", 1, nullptr, 'l' },
                                        { "quiet", 0, nullptr, 0 },
                                        { "verbose", 1, nullptr, 'v' },
                                        { 0, 0, nullptr, 0 } };

    int opt_index;
    int verbose = 1;
    bool fail = false;
    bool optimize = false;
    bool quiet = false;
    std::map<std::string, lp::parameter> parameters;

    while (not fail) {
        const auto opt =
          getopt_long(argc, argv, short_opts, long_opts, &opt_index);
        if (opt == -1)
            break;

        switch (opt) {
        case 0:
            break;
        case 'O':
            optimize = true;
            break;
        case 'l':
            parameters["limit"] = to_long(::optarg, 1000l);
            break;
        case 'h':
            help();
            return EXIT_SUCCESS;
        case 'p': {
            std::string name;
            lp::parameter value;
            std::tie(name, value) = split_param(::optarg);
            parameters[name] = value;
        } break;
        case '?':
        default:
            fail = true;
            std::fprintf(stderr, "Unknown command line option\n");
            break;
        };
    }

    if (fail)
        return EXIT_FAILURE;

    (void)verbose;
    (void)quiet;

    auto ctx = std::make_shared<lp::context>();
    ctx->set_standard_stream_logger();

    for (int i = ::optind; i < argc; ++i) {
        try {
            auto pb = lp::make_problem(ctx, argv[i]);

            std::string filename(argv[i]);
            filename += '-';
            filename += std::to_string(::getpid());
            filename += '_';
            filename += std::to_string(i);
            filename += ".sol";

            ctx->log(lp::context::message_type::info,
                     "solution: %s\n",
                     filename.c_str());

            std::ofstream ofs(filename);
            ofs << std::boolalpha
                << std::setprecision(std::floor(
                     std::numeric_limits<double>::digits * std::log10(2) + 2));

            auto now = std::chrono::system_clock::now();
            auto in_time_t = std::chrono::system_clock::to_time_t(now);

            ofs << "start: "
                << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X")
                << '\n';

            auto ret = solve(ctx, pb, parameters, optimize);

            if (ret.status == lp::result_status::success) {
                ofs << "Solution found: " << ret.value << '\n' << ret;
            } else {
                ofs << "Solution not found. Missing constraints: "
                    << ret.remaining_constraints << '\n';
            }
        } catch (const lp::precondition_error& e) {
            std::fprintf(stderr, "internal failure\n");
        } catch (const lp::postcondition_error& e) {
            std::fprintf(stderr, "internal failure\n");
        } catch (const lp::numeric_cast_error& e) {
            std::fprintf(stderr, "numeric cast interal failure\n");
        } catch (const lp::file_access_error& e) {
            std::fprintf(stderr,
                         "file `%s' fail %d: %s\n",
                         e.file().c_str(),
                         e.error(),
                         std::strerror(e.error()));
        } catch (const lp::file_format_error& e) {
            std::fprintf(stderr,
                         "file format error at line %d column %d "
                         "%s\n",
                         e.line(),
                         e.column(),
                         file_format_error_format(e.failure()));
        } catch (const lp::problem_definition_error& e) {
            std::fprintf(stderr,
                         "definition problem error at %s: %s\n",
                         e.element().c_str(),
                         problem_definition_error_format(e.failure()));
        } catch (const lp::solver_error& e) {
            std::fprintf(
              stderr, "solver error: %s\n", solver_error_format(e.failure()));
        } catch (const std::exception& e) {
            std::fprintf(stderr, "failure: %s.\n", e.what());
        }
    }

    return EXIT_SUCCESS;
}

void
help() noexcept
{
    std::fprintf(stdout,
                 "--help|-h                   This help message\n"
                 "--param|-p [name]:[value]   Add a new parameter (name is"
                 " [a-z][A-Z]_ value can be a double, an integer otherwise a"
                 " string.\n"
                 "--optimize|-O               Optimize model (default "
                 "feasibility search only)\n"
                 "--limit int                 Set limit\n"
                 "--quiet                     Remove any verbose message\n"
                 "--verbose|-v int            Set verbose level\n");
}

const char*
file_format_error_format(lp::file_format_error::tag failure) noexcept
{
    static const char* const tag[] = {
        "end of file",     "unknown",
        "already defined", "incomplete",
        "bad name",        "bad operator",
        "bad integer",     "bad objective function type",
        "bad bound",       "bad function element",
        "bad constraint"
    };

    switch (failure) {
    case lp::file_format_error::tag::end_of_file:
        return tag[0];
    case lp::file_format_error::tag::unknown:
        return tag[1];
    case lp::file_format_error::tag::already_defined:
        return tag[2];
    case lp::file_format_error::tag::incomplete:
        return tag[3];
    case lp::file_format_error::tag::bad_name:
        return tag[4];
    case lp::file_format_error::tag::bad_operator:
        return tag[5];
    case lp::file_format_error::tag::bad_integer:
        return tag[6];
    case lp::file_format_error::tag::bad_objective_function_type:
        return tag[7];
    case lp::file_format_error::tag::bad_bound:
        return tag[8];
    case lp::file_format_error::tag::bad_function_element:
        return tag[9];
    case lp::file_format_error::tag::bad_constraint:
        return tag[10];
    }

    return nullptr;
}

const char*
problem_definition_error_format(
  lp::problem_definition_error::tag failure) noexcept
{
    static const char* const tag[] = {
        "empty variables",
        "empty objective function",
        "variable not used",
        "bad bound",
        "multiple constraints with different value"
    };

    switch (failure) {
    case lp::problem_definition_error::tag::empty_variables:
        return tag[0];
    case lp::problem_definition_error::tag::empty_objective_function:
        return tag[1];
    case lp::problem_definition_error::tag::variable_not_used:
        return tag[2];
    case lp::problem_definition_error::tag::bad_bound:
        return tag[3];
    case lp::problem_definition_error::tag::multiple_constraint:
        return tag[4];
    }

    return nullptr;
}

const char*
solver_error_format(lp::solver_error::tag failure) noexcept
{
    static const char* const tag[] = { "no solver available",
                                       "unrealisable constraint",
                                       "not enough memory" };

    switch (failure) {
    case lp::solver_error::tag::no_solver_available:
        return tag[0];
    case lp::solver_error::tag::unrealisable_constraint:
        return tag[1];
    case lp::solver_error::tag::not_enough_memory:
        return tag[2];
    }

    return nullptr;
}

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

long
to_long(const char* s, long bad_value) noexcept
{
    char* c;
    errno = 0;
    long value = std::strtol(s, &c, 10);

    if ((errno == ERANGE and (value == LONG_MIN or value == LONG_MAX)) or
        (value == 0 and c == s))
        return bad_value;

    return value;
}

std::tuple<std::string, lp::parameter>
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

    auto valuel = to_long(value.c_str(), LONG_MIN);
    auto valued = to_double(value.c_str(), -HUGE_VAL);

    double tmp;
    if (valued != -HUGE_VAL and std::modf(valued, &tmp))
        return std::make_tuple(name, lp::parameter(valued));

    if (valuel != LONG_MIN)
        return std::make_tuple(name, lp::parameter(valuel));

    return std::make_tuple(name, lp::parameter(value));
}

lp::result
solve(std::shared_ptr<lp::context> ctx,
      lp::problem& pb,
      const std::map<std::string, lp::parameter>& params,
      bool optimize)
{
    if (optimize)
        return lp::optimize(ctx, pb, params);

    return lp::solve(ctx, pb, params);
}
