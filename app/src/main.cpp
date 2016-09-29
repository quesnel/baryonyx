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
#include <experimental/optional>

#include <lpcore>

#include <cstdio>
#include <cerrno>
#include <getopt.h>

namespace std {

using experimental::optional;
using experimental::make_optional;

};

void help()
{
    std::fprintf(stdout,
            "--help|-h        This help message\n"
            "--kappa|-k real  Set Kappa parameter\n"
            "--delta|-d real  Set Delta parameter\n"
            "--theta|-t real  Set Theta parameter\n"
            "--limit int      Set limit\n"
            "--quiet          Remove any verbose message\n"
            "--verbose|-v int Set verbose level\n");
}

std::optional<double> to_double(const char* s, double min, double max)
{
    char *c;
    errno = 0;
    double value = std::strtof(s, &c);

    if ((errno == ERANGE and (value == std::numeric_limits<double>::lowest()
                              or value == -std::numeric_limits<double>::max()))
        or (errno != 0 and value == 0)
        or (c == ::optarg)
        or (value < min)
        or (value > max))
        return std::make_optional(value);

    return std::optional<double>{};
}

std::optional<long> to_long(const char* s, long min, long max)
{
    char *c;
    errno = 0;
    long value = std::strtol(s, &c, 10);

    if ((errno == ERANGE and (value == std::numeric_limits<long>::lowest()
                              or value == std::numeric_limits<long>::max()))
        or (errno != 0 and value == 0)
        or (c == ::optarg)
        or (value < min)
        or (value > max))
        return std::make_optional(value);

    return std::optional<long>{};
}

int main(int argc, char *argv[])
{
    const char* const short_opts = "hk:d:t:l:qv:";
    const struct option long_opts[] = {
        {"help", 0, nullptr, 'h'},
        {"kappa", 1, nullptr, 'k'},
        {"delta", 1, nullptr, 'd'},
        {"theta", 1, nullptr, 't'},
        {"limit", 1, nullptr, 'l'},
        {"quiet", 0, nullptr, 0},
        {"verbose", 1, nullptr, 'v'},
        {0, 0, nullptr, 0}};

    int opt_index;
    double kappa = 0.001;
    double delta = 0.0001;
    double theta = 0.001;
    int verbose = 1;
    long limit = 1000;
    bool fail = false;
    bool quiet = false;

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
        case 'k':
            {
                auto val = to_double(::optarg, 0, 1);
                if (not val)
                    std::fprintf(stderr, "fail to convert parameter `%s' for "
                                 " parameter kappa (or k)\n", ::optarg);
                else
                    kappa = *val;
            }
            break;
        case 'd':
            {
                auto val = to_double(::optarg, 0, std::numeric_limits<double>::max());
                if (not val)
                    std::fprintf(stderr, "fail to convert parameter `%s' for parameter"
                                 " delta (or d)\n", ::optarg);
                else
                    delta = *val;
            }
            break;
        case 't':
            {
                auto val = to_double(::optarg, 0, 1);
                if (not val)
                    std::fprintf(stderr, "fail to convert parameter `%s' for parameter"
                                 " theta (or t)\n", ::optarg);
                else
                    theta = *val;
            }
            break;
        case 'l':
            {
                auto val = to_long(::optarg, 0, std::numeric_limits<long>::max());
                if (not val)
                    std::fprintf(stderr, "fail to convert parameter `%s' for parameter"
                                 " limit (or l)\n", ::optarg);
                else
                    limit = *val;
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

            std::map<std::string, lp::parameter> params;
            params["kappa"] = kappa;
            params["theta"] = theta;
            params["delta"] = delta;
            params["limit"] = limit;

            auto ret = lp::solve(pb, params);

            for (std::size_t i {0}, e {ret.variable_name.size()};
                 i != e; ++i)
                std::fprintf(stdout, "%s = %d\n",
                             ret.variable_name[i].c_str(),
                             ret.variable_value[i]);
        } catch(const std::exception& e) {
            std::fprintf(stderr, "failure: %s.\n", e.what());
        }
    }

    return EXIT_SUCCESS;
}
