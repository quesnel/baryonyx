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

#include "problem.hpp"
#include "resume.hpp"

#include <boost/ut.hpp>

#include <fstream>
#include <map>
#include <numeric>
#include <streambuf>
#include <string_view>

#include <baryonyx/core-compare>
#include <baryonyx/core>
#include <utility>

class membuf : public std::streambuf
{
private:
    using ios_base = std::ios_base;

protected:
    pos_type seekoff(off_type off,
                     ios_base::seekdir dir,
                     [[maybe_unused]] ios_base::openmode which =
                       ios_base::in | ios_base::out) override
    {
        if (dir == ios_base::cur)
            gbump(static_cast<int>(off));
        else if (dir == ios_base::end)
            setg(eback(), egptr() + off, egptr());
        else if (dir == ios_base::beg)
            setg(eback(), eback() + off, egptr());
        return gptr() - eback();
    }

    pos_type seekpos(pos_type sp, ios_base::openmode which) override
    {
        return seekoff(sp - pos_type(off_type(0)), ios_base::beg, which);
    }

public:
    membuf(const char* data, size_t size) noexcept
    {
        auto p = const_cast<char*>(data);
        this->setg(p, p, p + size);
    }

    membuf(std::string_view buffer) noexcept
      : membuf(buffer.data(), buffer.size())
    {
    }
};

struct imemstream
  : private virtual membuf
  , public std::istream
{
    imemstream(std::string_view buffer) noexcept
      : membuf(buffer)
      , std::istream(static_cast<std::streambuf*>(this))
    {
    }
};

int
main()
{
    using namespace boost::ut;

    "named objective"_test = [] {
        const char* example = "maximize\n"
                              "x0: +x1 + 2x2 + 3x3 - 100\n"
                              "end\n";

        auto ctx = baryonyx::make_context(4);
        std::istringstream iss(example);

        auto pb = baryonyx::make_problem(*ctx, iss);
        expect(pb);
        expect(pb.type == baryonyx::objective_function_type::maximize);
        expect(pb.objective.elements.size() == 3_ul);
        expect(pb.objective.elements[0].factor == 1_i);
        expect(pb.objective.elements[0].variable_index == 0_i);
        expect(pb.objective.elements[1].factor == 2_i);
        expect(pb.objective.elements[1].variable_index == 1_i);
        expect(pb.objective.elements[2].factor == 3_i);
        expect(pb.objective.elements[2].variable_index == 2_i);
        expect(pb.objective.value == -100.0_d);
    };

    "no-named objective"_test = [] {
        const char* example = "maximize\n"
                              "st: x1 + x2 + x3 = 1\n"
                              "end\n";

        auto ctx = baryonyx::make_context(4);
        std::istringstream iss(example);

        auto pb = baryonyx::make_problem(*ctx, iss);
        expect(pb);
        expect(pb.type == baryonyx::objective_function_type::maximize);
        expect(pb.objective.elements.size() == 0_ul);
        expect(pb.vars.names.size() == 3_ul);
        expect(pb.vars.values.size() == 3_ul);
        expect(pb.less_constraints.size() == 0_ul);
        expect(pb.greater_constraints.size() == 0_ul);
        expect(pb.equal_constraints.size() == 1_ul);
    };

    "small lp"_test = [] {
        const char* example_1 = "maximize\n"
                                "obj: x1 + 2x2 + 3x3 - 100\n"
                                "st\n"
                                "time:  -x1 + x2 + x3 <= 20\n"
                                "labor:  x1 - 3x2 + x3 <= 30\n"
                                "test: x1 - 3x2 + x3 <= -5\n"
                                "bounds\n"
                                "x1 <= 40\n"
                                "end\n";

        auto ctx = baryonyx::make_context(4);
        std::istringstream iss(example_1);

        auto pb = baryonyx::make_problem(*ctx, iss);
        expect(pb);
        expect(pb.type == baryonyx::objective_function_type::maximize);
        expect(pb.objective.elements.size() == 3_ul);
        expect(pb.objective.elements[0].factor == 1_i);
        expect(pb.objective.elements[0].variable_index == 0_i);
        expect(pb.objective.elements[1].factor == 2_i);
        expect(pb.objective.elements[1].variable_index == 1_i);
        expect(pb.objective.elements[2].factor == 3_i);
        expect(pb.objective.elements[2].variable_index == 2_i);
        expect(pb.objective.value == -100.0_d);

        expect(pb.vars.names.size() == 3_ul);
        expect(pb.vars.values.size() == 3_ul);
        expect(pb.less_constraints.size() == 3_ul);

        expect(pb.less_constraints[0].elements.size() == 3_ul);
        expect(pb.less_constraints[0].elements[0].factor == -1);
        expect(pb.less_constraints[0].elements[0].variable_index == 0_i);
        expect(pb.less_constraints[0].elements[1].factor == 1_i);
        expect(pb.less_constraints[0].elements[1].variable_index == 1_i);
        expect(pb.less_constraints[0].elements[2].factor == 1_i);
        expect(pb.less_constraints[0].elements[2].variable_index == 2_i);
        expect(pb.less_constraints[0].value == 20);

        expect(pb.less_constraints[1].elements.size() == 3_ul);
        expect(pb.less_constraints[1].elements[0].factor == 1_i);
        expect(pb.less_constraints[1].elements[0].variable_index == 0_i);
        expect(pb.less_constraints[1].elements[1].factor == -3_i);
        expect(pb.less_constraints[1].elements[1].variable_index == 1_i);
        expect(pb.less_constraints[1].elements[2].factor == 1_i);
        expect(pb.less_constraints[1].elements[2].variable_index == 2_i);
        expect(pb.less_constraints[1].value == 30);

        expect(pb.less_constraints[2].elements.size() == 3_ul);
        expect(pb.less_constraints[1].elements[0].factor == 1_i);
        expect(pb.less_constraints[2].elements[0].variable_index == 0_i);
        expect(pb.less_constraints[2].elements[1].factor == -3_i);
        expect(pb.less_constraints[2].elements[1].variable_index == 1_i);
        expect(pb.less_constraints[2].elements[2].factor == 1_i);
        expect(pb.less_constraints[2].elements[2].variable_index == 2_i);
        expect(pb.less_constraints[2].value == -5);

        expect(pb.vars.names[0] == "x1");
        expect(pb.vars.names[1] == "x2");
        expect(pb.vars.names[2] == "x3");

        expect(pb.vars.values[0].min == 0_i);
        expect(pb.vars.values[1].min == 0_i);
        expect(pb.vars.values[2].min == 0_i);

        expect(pb.vars.values[0].max == 40_i);
        expect(pb.vars.values[1].max == std::numeric_limits<int>::max());
        expect(pb.vars.values[2].max == std::numeric_limits<int>::max());
    };

    "quadratic lp"_test = [] {
        const char* example_1 = "maximize\n"
                                "obj: x1 + 2x2 + 3x3 "
                                "- [] /2 - 100\n"
                                "st\n"
                                "time:  -x1 + x2 + x3 <= 20\n"
                                "labor:  x1 - 3x2 + x3 <= 30\n"
                                "test: x1 - 3x2 + x3 <= -5\n"
                                "bounds\n"
                                "x1 <= 40\n"
                                "end\n";

        auto ctx = baryonyx::make_context(4);
        std::istringstream iss(example_1);

        auto pb = baryonyx::make_problem(*ctx, iss);

        expect(pb);
        expect(pb.type == baryonyx::objective_function_type::maximize);
        expect(pb.objective.elements.size() == 3);
        expect(pb.objective.elements[0].factor == 1);
        expect(pb.objective.elements[0].variable_index == 0);
        expect(pb.objective.elements[1].factor == 2);
        expect(pb.objective.elements[1].variable_index == 1);
        expect(pb.objective.elements[2].factor == 3);
        expect(pb.objective.elements[2].variable_index == 2);

        expect(pb.objective.qelements.size() == 0);

        expect(pb.objective.value == -100);

        expect(pb.vars.names.size() == 3);
        expect(pb.vars.values.size() == 3);
        expect(pb.less_constraints.size() == 3);

        expect(pb.less_constraints[0].elements.size() == 3);
        expect(pb.less_constraints[0].elements[0].factor == -1);
        expect(pb.less_constraints[0].elements[0].variable_index == 0);
        expect(pb.less_constraints[0].elements[1].factor == 1);
        expect(pb.less_constraints[0].elements[1].variable_index == 1);
        expect(pb.less_constraints[0].elements[2].factor == 1);
        expect(pb.less_constraints[0].elements[2].variable_index == 2);
        expect(pb.less_constraints[0].value == 20);

        expect(pb.less_constraints[1].elements.size() == 3);
        expect(pb.less_constraints[1].elements[0].factor == 1);
        expect(pb.less_constraints[1].elements[0].variable_index == 0);
        expect(pb.less_constraints[1].elements[1].factor == -3);
        expect(pb.less_constraints[1].elements[1].variable_index == 1);
        expect(pb.less_constraints[1].elements[2].factor == 1);
        expect(pb.less_constraints[1].elements[2].variable_index == 2);
        expect(pb.less_constraints[1].value == 30);

        expect(pb.less_constraints[2].elements.size() == 3);
        expect(pb.less_constraints[1].elements[0].factor == 1);
        expect(pb.less_constraints[2].elements[0].variable_index == 0);
        expect(pb.less_constraints[2].elements[1].factor == -3);
        expect(pb.less_constraints[2].elements[1].variable_index == 1);
        expect(pb.less_constraints[2].elements[2].factor == 1);
        expect(pb.less_constraints[2].elements[2].variable_index == 2);
        expect(pb.less_constraints[2].value == -5);

        expect(pb.vars.names[0] == "x1");
        expect(pb.vars.names[1] == "x2");
        expect(pb.vars.names[2] == "x3");

        expect(pb.vars.values[0].min == 0);
        expect(pb.vars.values[1].min == 0);
        expect(pb.vars.values[2].min == 0);

        expect(pb.vars.values[0].max == 40);
        expect(pb.vars.values[1].max == std::numeric_limits<int>::max());
        expect(pb.vars.values[2].max == std::numeric_limits<int>::max());
    };

    "quadratic 2"_test = [] {
        const char* example_1 = "maximize\n"
                                "obj: x1 + 2x2 + 3x3 "
                                "- [ -x1^2 - 17x1*x2 + 5x2^2] /2 - 100\n"
                                "st\n"
                                "time:  -x1 + x2 + x3 <= 20\n"
                                "labor:  x1 - 3x2 + x3 <= 30\n"
                                "test: x1 - 3x2 + x3 <= -5\n"
                                "bounds\n"
                                "x1 <= 40\n"
                                "end\n";

        auto ctx = baryonyx::make_context(4);
        std::istringstream iss(example_1);

        auto pb = baryonyx::make_problem(*ctx, iss);

        expect(pb);
        expect(pb.type == baryonyx::objective_function_type::maximize);
        expect(pb.objective.elements.size() == 3);
        expect(pb.objective.elements[0].factor == 1);
        expect(pb.objective.elements[0].variable_index == 0);
        expect(pb.objective.elements[1].factor == 2);
        expect(pb.objective.elements[1].variable_index == 1);
        expect(pb.objective.elements[2].factor == 3);
        expect(pb.objective.elements[2].variable_index == 2);

        expect(pb.objective.qelements.size() == 3);

        expect(pb.objective.qelements[0].factor == 1.0 / 2.0);
        expect(pb.objective.qelements[0].variable_index_a == 0);
        expect(pb.objective.qelements[0].variable_index_b == 0);
        expect(pb.objective.qelements[1].factor == 17.0 / 2.0);
        expect(pb.objective.qelements[1].variable_index_a == 0);
        expect(pb.objective.qelements[1].variable_index_b == 1);
        expect(pb.objective.qelements[2].factor == -5.0 / 2.0);
        expect(pb.objective.qelements[2].variable_index_a == 1);
        expect(pb.objective.qelements[2].variable_index_b == 1);

        expect(pb.objective.value == -100);

        expect(pb.vars.names.size() == 3);
        expect(pb.vars.values.size() == 3);
        expect(pb.less_constraints.size() == 3);

        expect(pb.less_constraints[0].elements.size() == 3);
        expect(pb.less_constraints[0].elements[0].factor == -1);
        expect(pb.less_constraints[0].elements[0].variable_index == 0);
        expect(pb.less_constraints[0].elements[1].factor == 1);
        expect(pb.less_constraints[0].elements[1].variable_index == 1);
        expect(pb.less_constraints[0].elements[2].factor == 1);
        expect(pb.less_constraints[0].elements[2].variable_index == 2);
        expect(pb.less_constraints[0].value == 20);

        expect(pb.less_constraints[1].elements.size() == 3);
        expect(pb.less_constraints[1].elements[0].factor == 1);
        expect(pb.less_constraints[1].elements[0].variable_index == 0);
        expect(pb.less_constraints[1].elements[1].factor == -3);
        expect(pb.less_constraints[1].elements[1].variable_index == 1);
        expect(pb.less_constraints[1].elements[2].factor == 1);
        expect(pb.less_constraints[1].elements[2].variable_index == 2);
        expect(pb.less_constraints[1].value == 30);

        expect(pb.less_constraints[2].elements.size() == 3);
        expect(pb.less_constraints[1].elements[0].factor == 1);
        expect(pb.less_constraints[2].elements[0].variable_index == 0);
        expect(pb.less_constraints[2].elements[1].factor == -3);
        expect(pb.less_constraints[2].elements[1].variable_index == 1);
        expect(pb.less_constraints[2].elements[2].factor == 1);
        expect(pb.less_constraints[2].elements[2].variable_index == 2);
        expect(pb.less_constraints[2].value == -5);

        expect(pb.vars.names[0] == "x1");
        expect(pb.vars.names[1] == "x2");
        expect(pb.vars.names[2] == "x3");

        expect(pb.vars.values[0].min == 0);
        expect(pb.vars.values[1].min == 0);
        expect(pb.vars.values[2].min == 0);

        expect(pb.vars.values[0].max == 40);
        expect(pb.vars.values[1].max == std::numeric_limits<int>::max());
        expect(pb.vars.values[2].max == std::numeric_limits<int>::max());
    };

    "quadratic 3"_test = [] {
        const char* example_1 = "maximize\n"
                                "obj: x1 + 2x2 + 3x3 "
                                "- [ -x1*x2 - 17x1*x2 + 5x1*x2] /2 - 100\n"
                                "st\n"
                                "time:  -x1 + x2 + x3 <= 20\n"
                                "labor:  x1 - 3x2 + x3 <= 30\n"
                                "test: x1 - 3x2 + x3 <= -5\n"
                                "bounds\n"
                                "x1 <= 40\n"
                                "end\n";

        auto ctx = baryonyx::make_context(4);
        std::istringstream iss(example_1);

        auto pb = baryonyx::make_problem(*ctx, iss);

        expect(pb);
        expect(pb.type == baryonyx::objective_function_type::maximize);
        expect(pb.objective.elements.size() == 3);
        expect(pb.objective.elements[0].factor == 1);
        expect(pb.objective.elements[0].variable_index == 0);
        expect(pb.objective.elements[1].factor == 2);
        expect(pb.objective.elements[1].variable_index == 1);
        expect(pb.objective.elements[2].factor == 3);
        expect(pb.objective.elements[2].variable_index == 2);

        expect(pb.objective.qelements.size() == 1);

        expect(pb.objective.qelements[0].factor == 13.0 / 2.0);
        expect(pb.objective.qelements[0].variable_index_a == 0);
        expect(pb.objective.qelements[0].variable_index_b == 1);

        expect(pb.objective.value == -100);

        expect(pb.vars.names.size() == 3);
        expect(pb.vars.values.size() == 3);
        expect(pb.less_constraints.size() == 3);

        expect(pb.less_constraints[0].elements.size() == 3);
        expect(pb.less_constraints[0].elements[0].factor == -1);
        expect(pb.less_constraints[0].elements[0].variable_index == 0);
        expect(pb.less_constraints[0].elements[1].factor == 1);
        expect(pb.less_constraints[0].elements[1].variable_index == 1);
        expect(pb.less_constraints[0].elements[2].factor == 1);
        expect(pb.less_constraints[0].elements[2].variable_index == 2);
        expect(pb.less_constraints[0].value == 20);

        expect(pb.less_constraints[1].elements.size() == 3);
        expect(pb.less_constraints[1].elements[0].factor == 1);
        expect(pb.less_constraints[1].elements[0].variable_index == 0);
        expect(pb.less_constraints[1].elements[1].factor == -3);
        expect(pb.less_constraints[1].elements[1].variable_index == 1);
        expect(pb.less_constraints[1].elements[2].factor == 1);
        expect(pb.less_constraints[1].elements[2].variable_index == 2);
        expect(pb.less_constraints[1].value == 30);

        expect(pb.less_constraints[2].elements.size() == 3);
        expect(pb.less_constraints[1].elements[0].factor == 1);
        expect(pb.less_constraints[2].elements[0].variable_index == 0);
        expect(pb.less_constraints[2].elements[1].factor == -3);
        expect(pb.less_constraints[2].elements[1].variable_index == 1);
        expect(pb.less_constraints[2].elements[2].factor == 1);
        expect(pb.less_constraints[2].elements[2].variable_index == 2);
        expect(pb.less_constraints[2].value == -5);

        expect(pb.vars.names[0] == "x1");
        expect(pb.vars.names[1] == "x2");
        expect(pb.vars.names[2] == "x3");

        expect(pb.vars.values[0].min == 0);
        expect(pb.vars.values[1].min == 0);
        expect(pb.vars.values[2].min == 0);

        expect(pb.vars.values[0].max == 40);
        expect(pb.vars.values[1].max == std::numeric_limits<int>::max());
        expect(pb.vars.values[2].max == std::numeric_limits<int>::max());
    };

    "quadratic 4"_test = [] {
        const char* example_1 = "maximize\n"
                                "obj: x1 + 2x2 + 3x3 "
                                "+ [] /2 - 100\n"
                                "st\n"
                                "time:  -x1 + x2 + x3 <= 20\n"
                                "labor:  x1 - 3x2 + x3 <= 30\n"
                                "test: x1 - 3x2 + x3 <= -5\n"
                                "bounds\n"
                                "x1 <= 40\n"
                                "end\n";

        auto ctx = baryonyx::make_context(4);
        std::istringstream iss(example_1);

        auto pb = baryonyx::make_problem(*ctx, iss);

        expect(pb);
        expect(pb.type == baryonyx::objective_function_type::maximize);
        expect(pb.objective.elements.size() == 3);
        expect(pb.objective.elements[0].factor == 1);
        expect(pb.objective.elements[0].variable_index == 0);
        expect(pb.objective.elements[1].factor == 2);
        expect(pb.objective.elements[1].variable_index == 1);
        expect(pb.objective.elements[2].factor == 3);
        expect(pb.objective.elements[2].variable_index == 2);

        expect(pb.objective.qelements.size() == 0);

        expect(pb.objective.value == -100);

        expect(pb.vars.names.size() == 3);
        expect(pb.vars.values.size() == 3);
        expect(pb.less_constraints.size() == 3);

        expect(pb.less_constraints[0].elements.size() == 3);
        expect(pb.less_constraints[0].elements[0].factor == -1);
        expect(pb.less_constraints[0].elements[0].variable_index == 0);
        expect(pb.less_constraints[0].elements[1].factor == 1);
        expect(pb.less_constraints[0].elements[1].variable_index == 1);
        expect(pb.less_constraints[0].elements[2].factor == 1);
        expect(pb.less_constraints[0].elements[2].variable_index == 2);
        expect(pb.less_constraints[0].value == 20);

        expect(pb.less_constraints[1].elements.size() == 3);
        expect(pb.less_constraints[1].elements[0].factor == 1);
        expect(pb.less_constraints[1].elements[0].variable_index == 0);
        expect(pb.less_constraints[1].elements[1].factor == -3);
        expect(pb.less_constraints[1].elements[1].variable_index == 1);
        expect(pb.less_constraints[1].elements[2].factor == 1);
        expect(pb.less_constraints[1].elements[2].variable_index == 2);
        expect(pb.less_constraints[1].value == 30);

        expect(pb.less_constraints[2].elements.size() == 3);
        expect(pb.less_constraints[1].elements[0].factor == 1);
        expect(pb.less_constraints[2].elements[0].variable_index == 0);
        expect(pb.less_constraints[2].elements[1].factor == -3);
        expect(pb.less_constraints[2].elements[1].variable_index == 1);
        expect(pb.less_constraints[2].elements[2].factor == 1);
        expect(pb.less_constraints[2].elements[2].variable_index == 2);
        expect(pb.less_constraints[2].value == -5);

        expect(pb.vars.names[0] == "x1");
        expect(pb.vars.names[1] == "x2");
        expect(pb.vars.names[2] == "x3");

        expect(pb.vars.values[0].min == 0);
        expect(pb.vars.values[1].min == 0);
        expect(pb.vars.values[2].min == 0);

        expect(pb.vars.values[0].max == 40);
        expect(pb.vars.values[1].max == std::numeric_limits<int>::max());
        expect(pb.vars.values[2].max == std::numeric_limits<int>::max());
    };

    "quadratic 5"_test = [] {
        const char* example_1 = "maximize\n"
                                "obj: x1 + 2x2 + 3x3 "
                                "+ [ -x1^2 - 17x1*x2 + 5x2^2] /2 - 100\n"
                                "st\n"
                                "time:  -x1 + x2 + x3 <= 20\n"
                                "labor:  x1 - 3x2 + x3 <= 30\n"
                                "test: x1 - 3x2 + x3 <= -5\n"
                                "bounds\n"
                                "x1 <= 40\n"
                                "end\n";

        auto ctx = baryonyx::make_context(4);
        std::istringstream iss(example_1);

        auto pb = baryonyx::make_problem(*ctx, iss);

        expect(pb);
        expect(pb.type == baryonyx::objective_function_type::maximize);
        expect(pb.objective.elements.size() == 3);
        expect(pb.objective.elements[0].factor == 1);
        expect(pb.objective.elements[0].variable_index == 0);
        expect(pb.objective.elements[1].factor == 2);
        expect(pb.objective.elements[1].variable_index == 1);
        expect(pb.objective.elements[2].factor == 3);
        expect(pb.objective.elements[2].variable_index == 2);

        expect(pb.objective.qelements.size() == 3);

        expect(pb.objective.qelements[0].factor == -1.0 / 2.0);
        expect(pb.objective.qelements[0].variable_index_a == 0);
        expect(pb.objective.qelements[0].variable_index_b == 0);
        expect(pb.objective.qelements[1].factor == -17.0 / 2.0);
        expect(pb.objective.qelements[1].variable_index_a == 0);
        expect(pb.objective.qelements[1].variable_index_b == 1);
        expect(pb.objective.qelements[2].factor == 5.0 / 2.0);
        expect(pb.objective.qelements[2].variable_index_a == 1);
        expect(pb.objective.qelements[2].variable_index_b == 1);

        expect(pb.objective.value == -100);

        expect(pb.vars.names.size() == 3);
        expect(pb.vars.values.size() == 3);
        expect(pb.less_constraints.size() == 3);

        expect(pb.less_constraints[0].elements.size() == 3);
        expect(pb.less_constraints[0].elements[0].factor == -1);
        expect(pb.less_constraints[0].elements[0].variable_index == 0);
        expect(pb.less_constraints[0].elements[1].factor == 1);
        expect(pb.less_constraints[0].elements[1].variable_index == 1);
        expect(pb.less_constraints[0].elements[2].factor == 1);
        expect(pb.less_constraints[0].elements[2].variable_index == 2);
        expect(pb.less_constraints[0].value == 20);

        expect(pb.less_constraints[1].elements.size() == 3);
        expect(pb.less_constraints[1].elements[0].factor == 1);
        expect(pb.less_constraints[1].elements[0].variable_index == 0);
        expect(pb.less_constraints[1].elements[1].factor == -3);
        expect(pb.less_constraints[1].elements[1].variable_index == 1);
        expect(pb.less_constraints[1].elements[2].factor == 1);
        expect(pb.less_constraints[1].elements[2].variable_index == 2);
        expect(pb.less_constraints[1].value == 30);

        expect(pb.less_constraints[2].elements.size() == 3);
        expect(pb.less_constraints[1].elements[0].factor == 1);
        expect(pb.less_constraints[2].elements[0].variable_index == 0);
        expect(pb.less_constraints[2].elements[1].factor == -3);
        expect(pb.less_constraints[2].elements[1].variable_index == 1);
        expect(pb.less_constraints[2].elements[2].factor == 1);
        expect(pb.less_constraints[2].elements[2].variable_index == 2);
        expect(pb.less_constraints[2].value == -5);

        expect(pb.vars.names[0] == "x1");
        expect(pb.vars.names[1] == "x2");
        expect(pb.vars.names[2] == "x3");

        expect(pb.vars.values[0].min == 0);
        expect(pb.vars.values[1].min == 0);
        expect(pb.vars.values[2].min == 0);

        expect(pb.vars.values[0].max == 40);
        expect(pb.vars.values[1].max == std::numeric_limits<int>::max());
        expect(pb.vars.values[2].max == std::numeric_limits<int>::max());
    };

    "quadratic 6"_test = [] {
        const char* example_1 = "maximize\n"
                                "obj: x1 + 2x2 + 3x3 "
                                "+ [ -x1*x2 - 17x1*x2 + 5x1*x2] /2 - 100\n"
                                "st\n"
                                "time:  -x1 + x2 + x3 <= 20\n"
                                "labor:  x1 - 3x2 + x3 <= 30\n"
                                "test: x1 - 3x2 + x3 <= -5\n"
                                "bounds\n"
                                "x1 <= 40\n"
                                "end\n";

        auto ctx = baryonyx::make_context(4);
        std::istringstream iss(example_1);

        auto pb = baryonyx::make_problem(*ctx, iss);

        expect(pb);
        expect(pb.type == baryonyx::objective_function_type::maximize);
        expect(pb.objective.elements.size() == 3);
        expect(pb.objective.elements[0].factor == 1);
        expect(pb.objective.elements[0].variable_index == 0);
        expect(pb.objective.elements[1].factor == 2);
        expect(pb.objective.elements[1].variable_index == 1);
        expect(pb.objective.elements[2].factor == 3);
        expect(pb.objective.elements[2].variable_index == 2);

        expect(pb.objective.qelements.size() == 1);

        expect(pb.objective.qelements[0].factor == -13.0 / 2.0);
        expect(pb.objective.qelements[0].variable_index_a == 0);
        expect(pb.objective.qelements[0].variable_index_b == 1);

        expect(pb.objective.value == -100);

        expect(pb.vars.names.size() == 3);
        expect(pb.vars.values.size() == 3);
        expect(pb.less_constraints.size() == 3);

        expect(pb.less_constraints[0].elements.size() == 3);
        expect(pb.less_constraints[0].elements[0].factor == -1);
        expect(pb.less_constraints[0].elements[0].variable_index == 0);
        expect(pb.less_constraints[0].elements[1].factor == 1);
        expect(pb.less_constraints[0].elements[1].variable_index == 1);
        expect(pb.less_constraints[0].elements[2].factor == 1);
        expect(pb.less_constraints[0].elements[2].variable_index == 2);
        expect(pb.less_constraints[0].value == 20);

        expect(pb.less_constraints[1].elements.size() == 3);
        expect(pb.less_constraints[1].elements[0].factor == 1);
        expect(pb.less_constraints[1].elements[0].variable_index == 0);
        expect(pb.less_constraints[1].elements[1].factor == -3);
        expect(pb.less_constraints[1].elements[1].variable_index == 1);
        expect(pb.less_constraints[1].elements[2].factor == 1);
        expect(pb.less_constraints[1].elements[2].variable_index == 2);
        expect(pb.less_constraints[1].value == 30);

        expect(pb.less_constraints[2].elements.size() == 3);
        expect(pb.less_constraints[1].elements[0].factor == 1);
        expect(pb.less_constraints[2].elements[0].variable_index == 0);
        expect(pb.less_constraints[2].elements[1].factor == -3);
        expect(pb.less_constraints[2].elements[1].variable_index == 1);
        expect(pb.less_constraints[2].elements[2].factor == 1);
        expect(pb.less_constraints[2].elements[2].variable_index == 2);
        expect(pb.less_constraints[2].value == -5);

        expect(pb.vars.names[0] == "x1");
        expect(pb.vars.names[1] == "x2");
        expect(pb.vars.names[2] == "x3");

        expect(pb.vars.values[0].min == 0);
        expect(pb.vars.values[1].min == 0);
        expect(pb.vars.values[2].min == 0);

        expect(pb.vars.values[0].max == 40);
        expect(pb.vars.values[1].max == std::numeric_limits<int>::max());
        expect(pb.vars.values[2].max == std::numeric_limits<int>::max());
    };

    "assignment_problem 0"_test = [] {
        auto ctx = baryonyx::make_context(4);
        std::ifstream ifs;

        std::vector<std::vector<int>> values(3);

        values[0] = { 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };
        values[1] = { 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0 };
        values[2] = { 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0 };

        for (int i = 1; i != 4; ++i) {
            std::string filepath{ EXAMPLES_DIR "/assignment_problem_" };
            filepath += std::to_string(i);
            filepath += ".lp";

            auto pb = baryonyx::make_problem(ctx, filepath);

            expect(pb);
            expect(pb.status == baryonyx::file_format_error_tag::success);
            expect(pb.vars.names.size() == 16);
            expect(pb.vars.values.size() == 16);

            fmt::memory_buffer buffer;
            fmt::format_to(std::back_inserter(buffer), "{}", pb);

            imemstream ss(std::string_view(buffer.data(), buffer.size()));
            auto pb2 = baryonyx::make_problem(*ctx, ss);

            expect(pb == pb2);
        }
    };

    "geom-30a-3-ext_1000_support"_test = [] {
        auto ctx = baryonyx::make_context(4);
        auto pb = baryonyx::make_problem(
          ctx, EXAMPLES_DIR "/geom-30a-3-ext_1000_support.lp");

        expect(pb.type == baryonyx::objective_function_type::minimize);
        expect(pb.vars.names.size() == 819);
        expect(pb.vars.values.size() == 819);

        int nb = 0;
        for (auto& elem : pb.vars.values)
            if (elem.type == baryonyx::variable_type::binary)
                ++nb;

        expect(nb == 90);
    };

    "general"_test = [] {
        auto ctx = baryonyx::make_context(4);
        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/general.lp");

        expect(pb.type == baryonyx::objective_function_type::minimize);
        expect(pb.vars.names.size() == 3);
        expect(pb.vars.values.size() == 3);

        int nb = 0;
        for (auto& elem : pb.vars.values)
            if (elem.type == baryonyx::variable_type::general)
                ++nb;

        expect(nb == 3);
    };

    "sudoku"_test = [] {
        auto ctx = baryonyx::make_context(4);
        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/sudoku.lp");

        expect(pb.vars.names.size() == 81);
        expect(pb.vars.values.size() == 81);

        for (auto& vv : pb.vars.values) {
            expect(vv.min == 1);
            expect(vv.max == 9);
            expect(vv.type == baryonyx::variable_type::general);
        }
    };

    "8_queens_puzzle"_test = [] {
        auto ctx = baryonyx::make_context(4);
        auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/8_queens_puzzle.lp");

        expect(pb.vars.names.size() == 64);
        expect(pb.vars.values.size() == 64);

        for (auto& vv : pb.vars.values) {
            expect(vv.min == 0);
            expect(vv.max == 1);
            expect(vv.type == baryonyx::variable_type::binary);
        }
    };

    "vm"_test = [] {
        auto ctx = baryonyx::make_context(4);
        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/vm.lp");
    };
}
