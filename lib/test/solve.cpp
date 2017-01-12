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

#include <lpcore>
#include <lpcore-compare>
#include <lpcore-out>
#include <iostream>
#include <fstream>
#include <numeric>
#include <map>
#include <sstream>
#include "unit-test.hpp"

void test_general_lp()
{
    auto pb = lp::make_problem(EXAMPLES_DIR "/general.lp");
    std::cout << __func__ << '\n' << lp::resume(pb) << '\n';

    Ensures(pb.type == lp::objective_function_type::minimize);
    Ensures(pb.vars.names.size() == 3);
    Ensures(pb.vars.values.size() == 3);

    lp::index nb {0};
    for (auto& elem : pb.vars.values)
        if (elem.type == lp::variable_type::general)
            ++nb;

    Ensures(nb == 3);

    std::map<std::string, lp::parameter> params;
    params["kappa"] = 0.5;
    params["theta"] = 0.5;
    params["delta"] = 0.5;
    params["limit"] = 1000l;

    auto result = lp::solve(pb, params);

    std::cout << result << '\n';
}

int main(int /* argc */, char */* argv */[])
{
    test_general_lp();

    return unit_test::report_errors();
}
