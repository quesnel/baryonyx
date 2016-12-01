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

#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <lpcore>

#include <cstdio>
#include <cerrno>
#include <cstring>
#include <getopt.h>

double to_double(const char* s, double bad_value);
long to_long(const char* s, long bad_value);

std::tuple<std::string, lp::parameter>
split_param(const char *param)
{
    std::string name, value;

    while (*param) {
        if (isalpha(*param) or *param == '_')
            name += *param;
        else
            break;

        param++;
    }

    if (*param and *param == ':') {
        param++;

        while (*param) {
            if (isalnum(*param) or *param == '_' or *param == ':')
                value += *param;
            else
                break;

            param++;
        }
    }

    auto valuel = to_long(::optarg, std::numeric_limits<long>::min());
    if (valuel != std::numeric_limits<long>::min())
        return std::make_tuple(name, lp::parameter(valuel));

    auto valued = to_double(::optarg, std::numeric_limits<double>::min());
    if (valued != std::numeric_limits<double>::min())
        return std::make_tuple(name, lp::parameter(valued));

    return std::make_tuple(name, lp::parameter(value));
}

const char* file_format_error_format(lp::file_format_error::tag);
const char* problem_definition_error_format(lp::problem_definition_error::tag);
const char* solver_error_format(lp::solver_error::tag);

void help()
{
    std::fprintf(stdout,
            "--help|-h                   This help message\n"
            "--param|-p [name]:[value]   Add a new parameter (name is"
            " [a-z][A-Z]_ value can be a double, an integer otherwise a"
            " string.\n"
            "--limit int                 Set limit\n"
            "--quiet                     Remove any verbose message\n"
            "--verbose|-v int            Set verbose level\n");
}


int main(int argc, char *argv[])
{
    const char* const short_opts = "hp:l:qv:";
    const struct option long_opts[] = {
        {"help", 0, nullptr, 'h'},
        {"param", 1, nullptr, 'p'},
        {"limit", 1, nullptr, 'l'},
        {"quiet", 0, nullptr, 0},
        {"verbose", 1, nullptr, 'v'},
        {0, 0, nullptr, 0}};

    int opt_index;
    int verbose = 1;
    bool fail = false;
    bool quiet = false;
    std::map<std::string, lp::parameter> parameters;

    while (not fail) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts,
                                     &opt_index);
        if (opt == -1)
            break;

        switch (opt) {
        case 0:
            break;
        case 'h':
            help();
            return EXIT_SUCCESS;
        case 'p':
            {
                std::string name;
                lp::parameter value;
                std::tie(name, value) = split_param(::optarg);

                parameters[name] = value;
            }
            break;
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

    for (int i = ::optind; i < argc; ++i) {
        try {
            auto pb = lp::make_problem(argv[i]);
            auto ret = lp::solve(pb, parameters);

            for (std::size_t i {0}, e {ret.variable_name.size()};
                 i != e; ++i)
                std::fprintf(stdout, "%s = %d\n",
                             ret.variable_name[i].c_str(),
                             ret.variable_value[i]);
        } catch(const lp::precondition_error& e) {
            std::fprintf(stderr, "internal failure\n");
        } catch(const lp::postcondition_error& e) {
            std::fprintf(stderr, "internal failure\n");
        } catch(const lp::numeric_cast_error& e) {
            std::fprintf(stderr, "numeric cast interal failure\n");
        } catch(const lp::file_access_error& e) {
            std::fprintf(stderr, "file `%s' fail %d: %s\n",
                e.file().c_str(), e.error(), std::strerror(e.error()));
        } catch(const lp::file_format_error& e) {
            std::fprintf(stderr, "file format error at line %d column %d "
                "%s\n", e.line(), e.column(),
                file_format_error_format(e.failure()));
        } catch(const lp::problem_definition_error& e) {
            std::fprintf(stderr, "definition problem error at %s: %s\n",
                e.element().c_str(),
                problem_definition_error_format(e.failure()));
        } catch(const lp::solver_error& e) {
            std::fprintf(stderr, "solver error: %s\n",
                solver_error_format(e.failure()));
        } catch(const std::exception& e) {
            std::fprintf(stderr, "failure: %s.\n", e.what());
        }
    }

    return EXIT_SUCCESS;
}

const char*
file_format_error_format(lp::file_format_error::tag failure)
{
    static const char *const tag[] = {
        "end of file",
        "unknown",
        "already defined",
        "incomplete",
        "bad name", "bad operator", "bad integer",
        "bad objective function type",
        "bad bound",
        "bad function element",
        "bad constraint"
    };

    switch(failure) {
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
    case lp::file_format_error::tag::bad_constraint :
        return tag[10];
    }

    return nullptr;
}

const char*
problem_definition_error_format(lp::problem_definition_error::tag failure)
{
    static const char *const tag[] = {
        "empty variables",
        "empty objective function",
        "variable not used", "bad bound"
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
    }

    return nullptr;
}

const char*
solver_error_format(lp::solver_error::tag failure)
{
    static const char *const tag[] = {
        "no solver available",
        "not enough memory"
    };

    switch (failure) {
    case lp::solver_error::tag::no_solver_available:
        return tag[0];
    case lp::solver_error::tag::not_enough_memory:
        return tag[1];
    }

    return nullptr;
}

double to_double(const char* s, double bad_value)
{
    char *c;
    errno = 0;
    double value = std::strtof(s, &c);

    if ((errno == ERANGE and (value == std::numeric_limits<double>::lowest()
                              or value == -std::numeric_limits<double>::max()))
        or (errno != 0 and value == 0)
        or (c == ::optarg))
        return value;

    return bad_value;
}

long to_long(const char* s, long bad_value)
{
    char *c;
    errno = 0;
    long value = std::strtol(s, &c, 10);

    if ((errno == ERANGE and (value == std::numeric_limits<long>::lowest()
                              or value == std::numeric_limits<long>::max()))
        or (errno != 0 and value == 0)
        or (c == ::optarg))
        return value;

    return bad_value;
}
