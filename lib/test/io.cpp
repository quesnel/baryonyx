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

#include <algorithm>
#include <fstream>
#include <map>
#include <numeric>
#include <streambuf>
#include <string_view>
#include <tuple>
#include <vector>

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

struct canonical_problem
{
    using constraint_key = std::tuple<std::string,
                                      std::vector<std::pair<std::string, int>>,
                                      char,
                                      int>;

    baryonyx::objective_function_type type;
    double constant;
    std::map<std::string, double> linear;
    std::map<std::pair<std::string, std::string>, double> quadratic;
    std::map<std::string, std::tuple<int, int, baryonyx::variable_type>>
      variables;
    std::vector<constraint_key> constraints;

    bool operator==(const canonical_problem&) const = default;
};

/// Builds an index-free view of a problem: the writer omits null
/// coefficients, so a re-read problem may number its variables differently.
static canonical_problem
make_canonical(const baryonyx::raw_problem& pb)
{
    canonical_problem ret{ pb.type, pb.objective.value, {}, {}, {}, {} };

    const auto name = [&pb](int index) {
        return std::string(pb.vars.names[static_cast<std::size_t>(index)]);
    };

    for (const auto& elem : pb.objective.elements)
        if (elem.factor != 0.0)
            ret.linear[name(elem.variable_index)] += elem.factor;

    for (const auto& elem : pb.objective.qelements) {
        if (elem.factor == 0.0)
            continue;

        auto a = name(elem.variable_index_a);
        auto b = name(elem.variable_index_b);
        if (b < a)
            std::swap(a, b);

        ret.quadratic[{ std::move(a), std::move(b) }] += elem.factor;
    }

    for (std::size_t i = 0, e = pb.vars.names.size(); i != e; ++i)
        ret.variables[std::string(pb.vars.names[i])] = { pb.vars.values[i].min,
                                                         pb.vars.values[i].max,
                                                         pb.vars.values[i].type };

    const auto add_constraints = [&](const auto& constraints, char op) {
        for (const auto& cst : constraints) {
            std::vector<std::pair<std::string, int>> terms;
            for (const auto& elem : cst.elements)
                terms.emplace_back(name(elem.variable_index), elem.factor);

            std::sort(terms.begin(), terms.end());
            ret.constraints.emplace_back(
              std::string(cst.label), std::move(terms), op, cst.value);
        }
    };

    add_constraints(pb.equal_constraints, '=');
    add_constraints(pb.greater_constraints, '>');
    add_constraints(pb.less_constraints, '<');
    std::sort(ret.constraints.begin(), ret.constraints.end());

    return ret;
}

static std::string
write_problem(const baryonyx::raw_problem& pb)
{
    return fmt::format("{}", pb);
}

int
main()
{
    using namespace boost::ut;

    "missing solution file"_test = [] {
        auto ctx = baryonyx::make_context(0);
        auto result =
          baryonyx::make_result(ctx, EXAMPLES_DIR "/does-not-exist.sol");

        expect(!result);
        expect(result.status == baryonyx::result_status::uninitialized);
    };

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

    const auto parse = [](std::string_view str) {
        auto ctx = baryonyx::make_context(0);
        std::istringstream iss{ std::string(str) };
        return baryonyx::make_problem(*ctx, iss);
    };

    const auto status = [&parse](std::string_view str) {
        return parse(str).status;
    };

    using baryonyx::file_format_error_tag;

    "parser long variable name"_test = [&parse] {
        const std::string name(2048, 'x');
        const auto pb = parse("minimize\n" + name + "\nst\n" + name +
                              " >= 1\nend\n");
        expect(pb);
        expect((pb.vars.names.size() == 1_ul) >> fatal);
        expect(pb.vars.names[0] == name);
    };

    "parser comments"_test = [&parse] {
        const auto pb = parse("\\ header comment\n"
                              "minimize\n"
                              "x + y \\ trailing comment\n"
                              "st\n"
                              "c1: x + y >= 1\\glued comment\n"
                              "end\n"
                              "\\ last comment");
        expect(pb);
        expect(pb.objective.elements.size() == 2_ul);
        expect(pb.greater_constraints.size() == 1_ul);
    };

    "parser scientific notation"_test = [&parse] {
        const auto pb = parse("minimize\n"
                              "1e-1 x + 2.5E+1 y - 3e2\n"
                              "end\n");
        expect(pb);
        expect((pb.objective.elements.size() == 2_ul) >> fatal);
        expect(pb.objective.elements[0].factor == 0.1_d);
        expect(pb.objective.elements[1].factor == 25.0_d);
        expect(pb.objective.value == -300.0_d);
    };

    "parser names looking like numbers"_test = [&parse] {
        const auto pb = parse("minimize\n"
                              "information + nancy + 2 infinite\n"
                              "end\n");
        expect(pb);
        expect((pb.objective.elements.size() == 3_ul) >> fatal);
        expect(pb.vars.names[0] == "information");
        expect(pb.vars.names[1] == "nancy");
        expect(pb.vars.names[2] == "infinite");
        expect(pb.objective.elements[2].factor == 2.0_d);
        expect(pb.objective.value == 0.0_d);
    };

    "parser spaced quadratic"_test = [&parse] {
        const auto pb = parse("minimize\n"
                              "+ [ 2 alpha * beta + 4 gamma ^ 2 "
                              "- 6 delta * epsilon ] / 2\n"
                              "end\n");
        expect(pb);
        expect((pb.objective.qelements.size() == 3_ul) >> fatal);
        expect((pb.vars.names.size() == 5_ul) >> fatal);
        expect(pb.vars.names[0] == "alpha");
        expect(pb.vars.names[1] == "beta");
        expect(pb.vars.names[2] == "gamma");
        expect(pb.vars.names[3] == "delta");
        expect(pb.vars.names[4] == "epsilon");
        expect(pb.objective.qelements[0].factor == 1.0_d);
        expect(pb.objective.qelements[1].factor == 2.0_d);
        expect(pb.objective.qelements[1].variable_index_a == 2_i);
        expect(pb.objective.qelements[1].variable_index_b == 2_i);
        expect(pb.objective.qelements[2].factor == -3.0_d);
    };

    "parser bad quadratic"_test = [&status] {
        expect(status("minimize\n[ x ] /2\nend\n") ==
               file_format_error_tag::bad_objective_quadratic);
        expect(status("minimize\n[ x ^ 3 ] /2\nend\n") ==
               file_format_error_tag::bad_objective_quadratic);
        expect(status("minimize\n[ x * y \nend\n") ==
               file_format_error_tag::bad_objective_quadratic);
    };

    "parser objective type with colon"_test = [&parse] {
        const auto pb = parse("max: 3 x\nend\n");
        expect(pb);
        expect(pb.type == baryonyx::objective_function_type::maximize);
        expect(pb.objective.elements.size() == 1_ul);
    };

    "parser subject to variants"_test = [&parse] {
        for (const auto* st : { "st", "s.t", "s.t.", "st.", "subject to",
                                "such that", "Such That:" }) {
            const auto pb =
              parse(fmt::format("minimize\nx\n{}\nx >= 1\nend\n", st));
            expect(pb) << st;
            expect(pb.vars.names.size() == 1_ul) << st;
            expect(pb.greater_constraints.size() == 1_ul) << st;
        }
    };

    "parser unterminated objective"_test = [&status] {
        expect(status("minimize\nx + y") == file_format_error_tag::bad_end);
        expect(status("minimize\nx + y <= 3\nend\n") ==
               file_format_error_tag::bad_objective);
        expect(status("minimize\nx - \nend\n") ==
               file_format_error_tag::bad_objective);
    };

    "parser bounds"_test = [&parse] {
        const auto pb = parse("minimize\n"
                              "a + b + c + d + e + f + g + h\n"
                              "bounds\n"
                              "a >= 2\n"
                              "b = 3\n"
                              "c free\n"
                              "2 <= d\n"
                              "e <= 4\n"
                              "e >= -1\n"
                              "-inf <= f <= 1e30\n"
                              "5 >= g\n"
                              "0.5 <= h <= 7.5\n"
                              "end\n");
        expect(pb);
        expect((pb.vars.values.size() == 8_ul) >> fatal);

        const auto& v = pb.vars.values;
        constexpr int int_min = std::numeric_limits<int>::min();
        constexpr int int_max = std::numeric_limits<int>::max();

        expect(v[0].min == 2_i && v[0].max == int_max);
        expect(v[1].min == 3_i && v[1].max == 3_i);
        expect(v[2].min == int_min && v[2].max == int_max);
        expect(v[3].min == 2_i && v[3].max == int_max);
        expect(v[4].min == -1_i && v[4].max == 4_i);
        expect(v[5].min == int_min && v[5].max == int_max);
        expect(v[6].min == 0_i && v[6].max == 5_i);
        expect(v[7].min == 1_i && v[7].max == 7_i);
    };

    "parser bad bounds"_test = [&status] {
        expect(status("minimize\nx\nbounds\n3 <= x <= 1\nend\n") ==
               file_format_error_tag::bad_bound);
        expect(status("minimize\nx\nbounds\nx <= \nend\n") ==
               file_format_error_tag::bad_bound);
        expect(status("minimize\nx\nbounds\nx >= inf\nend\n") ==
               file_format_error_tag::bad_bound);
    };

    "parser constant in constraint"_test = [&parse] {
        const auto pb = parse("minimize\nx + y\nst\nx + 2 + y - 1 <= 5\nend\n");
        expect(pb);
        expect(pb.vars.names.size() == 2_ul);
        expect((pb.less_constraints.size() == 1_ul) >> fatal);
        expect(pb.less_constraints[0].elements.size() == 2_ul);
        expect(pb.less_constraints[0].value == 4_i);
    };

    "parser non integer constraint"_test = [&parse, &status] {
        expect(status("minimize\nx\nst\n2.5 x >= 1\nend\n") ==
               file_format_error_tag::bad_constraint);
        expect(status("minimize\nx\nst\nx = 2.5\nend\n") ==
               file_format_error_tag::bad_constraint);
        expect(status("minimize\nx\nst\nx >= 1e20\nend\n") ==
               file_format_error_tag::bad_constraint);

        const auto pb =
          parse("minimize\nx\nst\nx >= 2.5\nx <= 2.5\n0.5 x + 0.5 x = 1\nend\n");
        expect(pb);
        expect((pb.greater_constraints.size() == 1_ul) >> fatal);
        expect(pb.greater_constraints[0].value == 3_i);
        expect((pb.less_constraints.size() == 1_ul) >> fatal);
        expect(pb.less_constraints[0].value == 2_i);
        expect((pb.equal_constraints.size() == 1_ul) >> fatal);
        expect((pb.equal_constraints[0].elements.size() == 1_ul) >> fatal);
        expect(pb.equal_constraints[0].elements[0].factor == 1_i);
    };

    "parser bad constraints"_test = [&status] {
        expect(status("minimize\nx\nst\n-5 <= x + y <= 5\nend\n") ==
               file_format_error_tag::bad_constraint);
        expect(status("minimize\nx\nst\nx + y <= 5 <= 6\nend\n") ==
               file_format_error_tag::bad_constraint);
        expect(status("minimize\nx\nst\nx <= y\nend\n") ==
               file_format_error_tag::bad_constraint);
        expect(status("minimize\nx\nst\nc1: 2 <= 3\nend\n") ==
               file_format_error_tag::bad_constraint);
        expect(status("minimize\nx\nst\nx + - <= 3\nend\n") ==
               file_format_error_tag::bad_constraint);
    };

    "parser sections in any order"_test = [&parse] {
        const auto pb = parse("minimize\n"
                              "x + y + z\n"
                              "st\n"
                              "x + y + z >= 1\n"
                              "general\n"
                              "y\n"
                              "binary\n"
                              "x w\n"
                              "bounds\n"
                              "y <= 4\n"
                              "end\n");
        expect(pb);
        expect((pb.vars.values.size() == 4_ul) >> fatal);
        expect(pb.vars.values[0].type == baryonyx::variable_type::binary);
        expect(pb.vars.values[1].type == baryonyx::variable_type::general);
        expect(pb.vars.values[1].max == 4_i);
        expect(pb.vars.values[2].type == baryonyx::variable_type::real);
        expect(pb.vars.names[3] == "w");
        expect(pb.vars.values[3].type == baryonyx::variable_type::binary);
    };

    "parser bad end"_test = [&status] {
        expect(status("minimize\nx\nend\nx\n") ==
               file_format_error_tag::bad_end_of_file);
        expect(status("minimize\nx\nbinary\n<=\nend\n") ==
               file_format_error_tag::bad_binary);
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

    "QPLIB_3852 read write read"_test = [] {
        auto ctx = baryonyx::make_context(0);
        const auto pb =
          baryonyx::make_problem(ctx, EXAMPLES_DIR "/QPLIB_3852.lp");

        expect((pb.status == file_format_error_tag::success) >> fatal);
        expect(pb.type == baryonyx::objective_function_type::maximize);
        expect(pb.vars.names.size() == 231_ul);
        expect(pb.objective.elements.size() == 231_ul);
        expect(pb.objective.qelements.size() == 440_ul);
        expect(pb.equal_constraints.empty());
        expect(pb.greater_constraints.empty());
        expect(pb.less_constraints.empty());

        for (const auto& value : pb.vars.values) {
            expect(value.type == baryonyx::variable_type::binary);
            expect(value.min == 0_i);
            expect(value.max == 1_i);
        }

        const auto written = write_problem(pb);
        std::istringstream iss(written);
        const auto pb2 = baryonyx::make_problem(*ctx, iss);

        expect((pb2.status == file_format_error_tag::success) >> fatal);
        expect(make_canonical(pb) == make_canonical(pb2));

        const auto rewritten = write_problem(pb2);
        std::istringstream iss2(rewritten);
        const auto pb3 = baryonyx::make_problem(*ctx, iss2);

        expect((pb3.status == file_format_error_tag::success) >> fatal);
        expect(pb2 == pb3);
        expect(make_canonical(pb2) == make_canonical(pb3));
        expect(write_problem(pb3) == rewritten);
    };

    "vm"_test = [] {
        auto ctx = baryonyx::make_context(4);
        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/vm.lp");
    };
}
