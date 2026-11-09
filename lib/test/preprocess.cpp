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
#include "utils.hpp"

#include <algorithm>
#include <functional>
#include <optional>
#include <random>
#include <sstream>
#include <string>

#include <boost/ut.hpp>

#include <fmt/format.h>

#include <baryonyx/core>

static int
get_variable(const std::vector<std::string_view>& variable_names,
             std::string_view name)
{
    for (int i = 0, e = baryonyx::length(variable_names); i != e; ++i)
        if (variable_names[i] == name)
            return i;

    return -1;
}

// Returns the best objective value over every 0-1 assignment, or nullopt when
// no assignment satisfies the constraints.
template<typename Problem>
static std::optional<double>
brute_force(const Problem& pb)
{
    const int n = baryonyx::length(pb.vars.names);
    std::optional<double> best;

    auto valid = [](const auto& csts, const std::vector<int>& x, auto op) {
        for (const auto& cst : csts) {
            int sum = 0;
            for (const auto& elem : cst.elements)
                sum += elem.factor * x[elem.variable_index];
            if (!op(sum, cst.value))
                return false;
        }
        return true;
    };

    std::vector<int> x(n);
    for (unsigned mask = 0; mask != (1u << n); ++mask) {
        bool in_bounds = true;
        for (int i = 0; i != n; ++i) {
            x[i] = (mask >> i) & 1u;
            in_bounds = in_bounds && pb.vars.values[i].min <= x[i] &&
                        x[i] <= pb.vars.values[i].max;
        }

        if (!in_bounds)
            continue;

        if (!valid(pb.equal_constraints, x, std::equal_to<>{}) ||
            !valid(pb.greater_constraints, x, std::greater_equal<>{}) ||
            !valid(pb.less_constraints, x, std::less_equal<>{}))
            continue;

        double value = pb.objective.value;
        for (const auto& elem : pb.objective.elements)
            value += elem.factor * x[elem.variable_index];
        for (const auto& elem : pb.objective.qelements)
            value += elem.factor * x[elem.variable_index_a] *
                     x[elem.variable_index_b];

        if (!best || value < *best)
            best = value;
    }

    return best;
}

static std::string
random_problem(std::mt19937& rng)
{
    std::uniform_int_distribution<int> factor(-2, 2);
    std::uniform_int_distribution<int> rhs(-2, 3);
    std::uniform_int_distribution<int> kind(0, 2);
    std::uniform_int_distribution<int> var(0, 5);
    std::uniform_int_distribution<int> length(1, 4);
    std::bernoulli_distribution fixed(0.1);

    auto term = [](int f, int v) {
        return fmt::format("{}{} x{} ", f < 0 ? "- " : "+ ", std::abs(f), v);
    };

    std::string ret = "minimize\n";
    for (int v = 0; v != 6; ++v)
        if (int f = factor(rng); f != 0)
            ret += term(f, v);
    ret += "+ [ 2 x0 * x1 ] / 2\nst\n";

    for (int c = 0; c != 4; ++c) {
        std::vector<int> used;
        for (int i = 0, e = length(rng); i != e; ++i) {
            int v = var(rng);
            if (std::ranges::find(used, v) != used.end())
                continue;
            used.push_back(v);

            int f = factor(rng);
            ret += term(f == 0 ? 1 : f, v);
        }
        constexpr const char* ops[] = { "=", ">=", "<=" };
        ret += fmt::format("{} {}\n", ops[kind(rng)], rhs(rng));
    }

    ret += "bounds\n";
    for (int v = 0; v != 6; ++v)
        if (fixed(rng))
            ret += fmt::format("x{} = {}\n", v, v % 2);

    ret += "general\nx0 x1 x2 x3 x4 x5\nend\n";
    return ret;
}

int
main()
{
    using namespace boost::ut;

    "test_bqp_clean"_test = [] {
        auto ctx = baryonyx::make_context(4);

        {
            const char* str_pb = "minimize\n"
                                 "-5a + [ 2a * b + 6c * d ] /2 + 10b\n"
                                 "Subject to:\n"
                                 "- a + b >= 0\n"
                                 "b + c - d >= 0\n"
                                 "Binaries\n"
                                 "a b c d\n"
                                 "End\n";

            std::istringstream iss(str_pb);
            auto pb = baryonyx::make_problem(*ctx, iss);
            expect(pb);

            auto pb_pp = baryonyx::preprocess(*ctx, pb);

            expect(pb_pp.vars.names.size() == 4u);
            expect(pb_pp.affected_vars.names.size() == 0u);
            expect(pb_pp.equal_constraints.size() == 0u);
            expect(pb_pp.less_constraints.size() == 0u);
            expect(pb_pp.greater_constraints.size() == 2u);

            auto ida = get_variable(pb.vars.names, "a");
            auto idb = get_variable(pb.vars.names, "b");
            auto idc = get_variable(pb.vars.names, "c");
            auto idd = get_variable(pb.vars.names, "d");

            expect(ida >= 0 && ida < baryonyx::length(pb.vars.names));
            expect(idb >= 0 && idb < baryonyx::length(pb.vars.names));
            expect(idc >= 0 && idc < baryonyx::length(pb.vars.names));
            expect(idd >= 0 && idd < baryonyx::length(pb.vars.names));

            expect(pb_pp.objective.elements.size() == 2u);
            expect(pb_pp.objective.elements[0].variable_index == ida);
            expect(pb_pp.objective.elements[0].factor == -5);
            expect(pb_pp.objective.elements[1].variable_index == idb);
            expect(pb_pp.objective.elements[1].factor == 10);

            expect(pb_pp.objective.qelements.size() == 2u);
            expect(pb_pp.objective.qelements[0].variable_index_a == ida);
            expect(pb_pp.objective.qelements[0].variable_index_b == idb);
            expect(pb_pp.objective.qelements[0].factor == 1);
            expect(pb_pp.objective.qelements[1].variable_index_a == idc);
            expect(pb_pp.objective.qelements[1].variable_index_b == idd);
            expect(pb_pp.objective.qelements[1].factor == 3);
        }

        {
            const char* str_pb = "minimize\n"
                                 "-5a + [ 2a * b + 6c * d ] /2 + 10b\n"
                                 "Subject to:\n"
                                 "- a + b >= 0\n"
                                 "b + c - d >= 0\n"
                                 "Binaries\n"
                                 "a b c d\n"
                                 "End\n";

            std::istringstream iss(str_pb);
            auto pb = baryonyx::make_problem(*ctx, iss);
            expect(pb);

            auto pb_pp = baryonyx::preprocess(*ctx, pb);

            expect(pb_pp.vars.names.size() == 4u);
            expect(pb_pp.affected_vars.names.size() == 0u);
            expect(pb_pp.equal_constraints.size() == 0u);
            expect(pb_pp.less_constraints.size() == 0u);
            expect(pb_pp.greater_constraints.size() == 2u);

            auto ida = get_variable(pb.vars.names, "a");
            auto idb = get_variable(pb.vars.names, "b");
            auto idc = get_variable(pb.vars.names, "c");
            auto idd = get_variable(pb.vars.names, "d");

            expect(ida >= 0 && ida < baryonyx::length(pb.vars.names));
            expect(idb >= 0 && idb < baryonyx::length(pb.vars.names));
            expect(idc >= 0 && idc < baryonyx::length(pb.vars.names));
            expect(idd >= 0 && idd < baryonyx::length(pb.vars.names));

            expect(pb_pp.objective.elements.size() == 2u);
            expect(pb_pp.objective.elements[0].variable_index == ida);
            expect(pb_pp.objective.elements[0].factor == -5);
            expect(pb_pp.objective.elements[1].variable_index == idb);
            expect(pb_pp.objective.elements[1].factor == 10);

            expect(pb_pp.objective.qelements.size() == 2u);
            expect(pb_pp.objective.qelements[0].variable_index_a == ida);
            expect(pb_pp.objective.qelements[0].variable_index_b == idb);
            expect(pb_pp.objective.qelements[0].factor == 1);
            expect(pb_pp.objective.qelements[1].variable_index_a == idc);
            expect(pb_pp.objective.qelements[1].variable_index_b == idd);
            expect(pb_pp.objective.qelements[1].factor == 3);
        }

        {
            const char* str_pb =
              "minimize\n"
              "-5a + 10a + [ 2a * b + 6c * d ] /2 + 10b - b -b "
              "-b + [2 a * b + 6 c * d] / 2\n"
              "Subject to:\n"
              "- a + b >= 0\n"
              "b + c - d >= 0\n"
              "Binaries\n"
              "a b c d\n"
              "End\n";

            std::istringstream iss(str_pb);
            auto pb = baryonyx::make_problem(*ctx, iss);
            expect(pb);

            auto pb_pp = baryonyx::preprocess(*ctx, pb);

            expect(pb_pp.vars.names.size() == 4u);
            expect(pb_pp.affected_vars.names.size() == 0u);
            expect(pb_pp.equal_constraints.size() == 0u);
            expect(pb_pp.less_constraints.size() == 0u);
            expect(pb_pp.greater_constraints.size() == 2u);

            auto ida = get_variable(pb.vars.names, "a");
            auto idb = get_variable(pb.vars.names, "b");
            auto idc = get_variable(pb.vars.names, "c");
            auto idd = get_variable(pb.vars.names, "d");

            expect(ida >= 0 && ida < baryonyx::length(pb.vars.names));
            expect(idb >= 0 && idb < baryonyx::length(pb.vars.names));
            expect(idc >= 0 && idc < baryonyx::length(pb.vars.names));
            expect(idd >= 0 && idd < baryonyx::length(pb.vars.names));

            expect(pb_pp.objective.elements.size() == 2u);
            expect(pb_pp.objective.elements[0].variable_index == ida);
            expect(pb_pp.objective.elements[0].factor == 5);
            expect(pb_pp.objective.elements[1].variable_index == idb);
            expect(pb_pp.objective.elements[1].factor == 7);

            expect(pb_pp.objective.qelements.size() == 2u);
            expect(pb_pp.objective.qelements[0].variable_index_a == ida);
            expect(pb_pp.objective.qelements[0].variable_index_b == idb);
            expect(pb_pp.objective.qelements[0].factor == 2);
            expect(pb_pp.objective.qelements[1].variable_index_a == idc);
            expect(pb_pp.objective.qelements[1].variable_index_b == idd);
            expect(pb_pp.objective.qelements[1].factor == 6);
        }
    };

    "test_bound_affectation"_test = [] {
        auto ctx = baryonyx::make_context(4);

        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/bound.lp");

        expect(pb);
        expect(pb.vars.names.size() == static_cast<size_t>(6));

        // The preprocessing can not remove any constraint or variable from the
        // raw_problem.

        auto pb_pp = baryonyx::preprocess(*ctx, pb);
        expect(pb_pp.vars.names.size() == static_cast<size_t>(6));
        expect(pb_pp.affected_vars.names.size() == static_cast<size_t>(0));
        expect(pb_pp.equal_constraints.size() == static_cast<size_t>(0));
        expect(pb_pp.less_constraints.size() == static_cast<size_t>(2));
        expect(pb_pp.greater_constraints.size() == static_cast<size_t>(0));

        // We add an "affectation bound" for the variable 'f' i.e. in bound
        // section of the lp file format:
        //
        // bounds:
        // f = 0

        auto id = get_variable(pb.vars.names, "f");
        expect(id >= 0);
        expect(id < baryonyx::length(pb.vars.names));

        expect(pb.vars.names.size() == pb.vars.values.size());

        pb.vars.values[id].min = 0;
        pb.vars.values[id].max = 0;

        auto pb_pp2 = baryonyx::preprocess(*ctx, pb);
        expect(pb_pp2.vars.names.size() == static_cast<size_t>(5));
        expect(pb_pp2.affected_vars.names.size() == static_cast<size_t>(1));
        expect(pb_pp2.equal_constraints.size() == static_cast<size_t>(0));
        expect(pb_pp2.less_constraints.size() == static_cast<size_t>(2));
        expect(pb_pp2.greater_constraints.size() == static_cast<size_t>(0));
    };

    "test_cleaning_affected_variables"_test = [] {
        auto ctx = baryonyx::make_context(4);
        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/prepro.lp");

        expect(pb);
        expect(pb.vars.names.size() == static_cast<size_t>(23));

        auto pb_pp = baryonyx::preprocess(*ctx, pb);

        expect(pb_pp.vars.names.size() == static_cast<size_t>(2));
        expect(pb_pp.affected_vars.names.size() == static_cast<size_t>(21));
    };

    "test_affect_variable"_test = [] {
        auto ctx = baryonyx::make_context(4);
        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/prepro.lp");

        expect(pb);
        auto pb_pp = baryonyx::preprocess(*ctx, pb);

        expect(pb_pp.vars.names.size() == static_cast<size_t>(2));
        expect(pb_pp.affected_vars.names.size() == static_cast<size_t>(21));

        {
            // + d + e <= 1
            // d = 1 implies e = 0

            auto pb_af = baryonyx::affect(*ctx, pb_pp, 0, true);
            expect(pb_af.vars.names.size() == static_cast<size_t>(0));
            expect(pb_af.affected_vars.names.size() ==
                   static_cast<size_t>(23));
        }

        {
            // + d + e <= 1
            // d = 0 implies e = {0, 1} but minimize so => 0.

            auto pb_af = baryonyx::affect(*ctx, pb_pp, 0, false);
            expect(pb_af.vars.names.size() == static_cast<size_t>(0));
            expect(pb_af.affected_vars.names.size() ==
                   static_cast<size_t>(23));
        }

        {
            // + d + e <= 1
            // e = 1 implies d = 0

            auto pb_af = baryonyx::affect(*ctx, pb_pp, 1, true);
            expect(pb_af.vars.names.size() == static_cast<size_t>(0));
            expect(pb_af.affected_vars.names.size() ==
                   static_cast<size_t>(23));
        }

        {
            // + d + e <= 1
            // e = 0 implies d = {0, 1} but minimize so => 0.

            auto pb_af = baryonyx::affect(*ctx, pb_pp, 1, false);
            expect(pb_af.vars.names.size() == static_cast<size_t>(0));
            expect(pb_af.affected_vars.names.size() ==
                   static_cast<size_t>(23));
        }
    };

    "test_split"_test = [] {
        auto ctx = baryonyx::make_context(4);
        auto pb = baryonyx::make_problem(ctx, EXAMPLES_DIR "/prepro.lp");
        expect(pb);
        auto pb_pp = baryonyx::preprocess(*ctx, pb);

        expect(pb_pp.vars.names.size() == static_cast<size_t>(2));
        expect(pb_pp.affected_vars.names.size() == static_cast<size_t>(21));

        baryonyx::problem p0, p1;
        std::tie(p0, p1) = baryonyx::split(*ctx, pb_pp, 0);

        {
            // + d + e <= 1
            // d = 1 implies e = 0

            expect(p0.vars.names.size() == static_cast<size_t>(0));
            expect(p0.affected_vars.names.size() == static_cast<size_t>(23));
        }

        {
            // + d + e <= 1
            // d = 0 implies e = {0, 1}

            expect(p1.vars.names.size() == static_cast<size_t>(0));
            expect(p1.affected_vars.names.size() == static_cast<size_t>(23));
        }
    };

    "test_unrealisable_problems"_test = [] {
        auto unrealisable = [](const char* str_pb) {
            auto ctx = baryonyx::make_context(4);
            std::istringstream iss(str_pb);
            auto pb = baryonyx::make_problem(*ctx, iss);

            try {
                baryonyx::preprocess(*ctx, pb);
            } catch (const baryonyx::solver_failure& e) {
                return e.failure() ==
                       baryonyx::solver_error_tag::unrealisable_constraint;
            }

            return false;
        };

        // The first constraint fixes every variable to 1, which violates the
        // second one.
        expect(unrealisable("minimize\n-x1 + x3\nst\n"
                            "x1 + x2 + x3 = 3\nx1 + x2 = 1\n"
                            "binary\nx1 x2 x3\nend\n"));

        expect(unrealisable("minimize\nx + y\nst\nx + y <= -1\n"
                            "binary\nx y\nend\n"));

        expect(unrealisable("minimize\nx + y\nst\nx + y >= 3\n"
                            "binary\nx y\nend\n"));

        expect(unrealisable(
          "minimize\nx + y\nst\nx + y <= 0\n"
          "bounds\nx >= 1\nx <= 1\nbinary\ny\ngeneral\nx\nend\n"));

        expect(unrealisable("minimize\nx + y\nst\nx + y >= 1\n"
                            "bounds\nx >= 2\nbinary\ny\ngeneral\nx\nend\n"));

        expect(unrealisable("minimize\nx + y + z\nst\n2 x = 1\ny + z >= 1\n"
                            "binary\nx y z\nend\n"));
    };

    "test_propagation"_test = [] {
        auto ctx = baryonyx::make_context(4);

        const char* str_pb = "minimize\n"
                             "a + b + c + d + e + f\n"
                             "st\n"
                             "c1: 2 a >= 1\n"
                             "c2: b - c = 1\n"
                             "c3: d + e + f >= 1\n"
                             "c4: d + e - f <= 1\n"
                             "bounds\n"
                             "d >= 1\n"
                             "binary\n"
                             "a b c e f\n"
                             "general\n"
                             "d\n"
                             "end\n";

        std::istringstream iss(str_pb);
        auto pb = baryonyx::make_problem(*ctx, iss);
        auto pb_pp = baryonyx::preprocess(*ctx, pb);

        auto value_of = [&](std::string_view name) {
            for (std::size_t i = 0; i != pb_pp.affected_vars.names.size(); ++i)
                if (pb_pp.affected_vars.names[i] == name)
                    return static_cast<int>(pb_pp.affected_vars.values[i]);
            return -1;
        };

        expect(value_of("a") == 1);
        expect(value_of("b") == 1);
        expect(value_of("c") == 0);
        expect(value_of("d") == 1);

        // d = 1 satisfies c3 and leaves e - f <= 0 in c4.
        expect(pb_pp.vars.names.size() == 2u);
        expect(pb_pp.less_constraints.size() == 1u);
        expect(pb_pp.greater_constraints.size() == 0u);
        expect(pb_pp.equal_constraints.size() == 0u);
        expect(pb_pp.objective.value == 3.0);
    };

    "test_quadratic_objective_cleaning"_test = [] {
        auto ctx = baryonyx::make_context(4);

        const char* str_pb = "minimize\n"
                             "x + y + [ 2 a * b + 4 c * d ] / 2\n"
                             "st\n"
                             "x + y >= 1\n"
                             "c <= 0\n"
                             "binary\n"
                             "x y a b c d\n"
                             "end\n";

        std::istringstream iss(str_pb);
        auto pb = baryonyx::make_problem(*ctx, iss);
        auto pb_pp = baryonyx::preprocess(*ctx, pb);

        // a and b only appear in a product: their best value depends on each
        // other, so they stay in the problem. c = 0 cancels the c * d product.
        expect(pb_pp.vars.names.size() == 5u);
        expect(pb_pp.affected_vars.names.size() == 1u);
        expect(pb_pp.affected_vars.names[0] == "c");

        expect(pb_pp.objective.qelements.size() == 1u);
        expect(pb_pp.objective.qelements[0].factor == 1.0);
        expect(pb_pp.objective.elements.size() == 2u);
        expect(pb_pp.objective.value == 0.0);
    };

    "test_split_on_unrealisable_side"_test = [] {
        auto ctx = baryonyx::make_context(4);

        const char* str_pb = "minimize\n"
                             "x + y + z\n"
                             "st\n"
                             "x + y + z >= 2\n"
                             "x + y + z <= 2\n"
                             "y + z <= 1\n"
                             "binary\n"
                             "x y z\n"
                             "end\n";

        std::istringstream iss(str_pb);
        auto raw = baryonyx::make_problem(*ctx, iss);
        auto pb = baryonyx::unpreprocess(*ctx, raw);

        auto x = get_variable(pb.vars.names, "x");
        expect(x >= 0);

        auto pb_1 = baryonyx::affect(*ctx, pb, x, true);
        expect(pb_1.vars.names.size() == 2u);

        expect(throws<baryonyx::solver_failure>(
          [&] { baryonyx::affect(*ctx, pb, x, false); }));
    };

    "test_random_problems_against_brute_force"_test = [] {
        auto ctx = baryonyx::make_context(2);
        std::mt19937 rng(20261109);

        for (int i = 0; i != 2000; ++i) {
            const auto str = random_problem(rng);
            std::istringstream iss(str);
            auto pb = baryonyx::make_problem(*ctx, iss);
            for (auto& v : pb.vars.values)
                v.max = std::min(v.max, 1);

            const auto expected = brute_force(pb);

            try {
                const auto pp = baryonyx::preprocess(*ctx, pb);
                const auto got = brute_force(pp);
                expect(got == expected)
                  << str << "expected:" << expected.value_or(-999)
                  << "got:" << got.value_or(-999);

                if (pp.vars.names.empty())
                    continue;

                const auto origin =
                  get_variable(pb.vars.names, pp.vars.names[0]);

                for (int value : { 0, 1 }) {
                    auto fixed = pb;
                    fixed.vars.values[origin].min = value;
                    fixed.vars.values[origin].max = value;
                    const auto expected_fixed = brute_force(fixed);

                    try {
                        const auto af =
                          baryonyx::affect(*ctx, pp, 0, value == 1);
                        expect(brute_force(af) == expected_fixed) << str;
                    } catch (const baryonyx::solver_failure&) {
                        expect(!expected_fixed.has_value()) << str;
                    }
                }
            } catch (const baryonyx::solver_failure&) {
                expect(!expected.has_value()) << str;
            }
        }
    };
}
