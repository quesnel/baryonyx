/* Copyright (C) 2016-2026 INRAE
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

#include <algorithm>
#include <utility>
#include <vector>

namespace baryonyx {

void
clear(raw_problem& pb)
{
    std::vector<objective_function_element>().swap(pb.objective.elements);
    std::vector<objective_quadratic_element>().swap(pb.objective.qelements);

    std::vector<constraint>().swap(pb.equal_constraints);
    std::vector<constraint>().swap(pb.greater_constraints);
    std::vector<constraint>().swap(pb.less_constraints);

    std::vector<std::string_view>().swap(pb.vars.names);
    std::vector<variable_value>().swap(pb.vars.values);

    pb.type = objective_function_type::maximize;
}

static std::vector<std::int8_t>
unconstrained_assignment(const problem& pb)
{
    const auto n = pb.vars.values.size();
    const double sign =
      pb.type == objective_function_type::maximize ? -1.0 : 1.0;

    // Minimized cost of setting each variable to 1, the others fixed: the
    // linear part plus the quadratic terms with the variables set to 1.
    std::vector<double> linear(n, 0.0);
    std::vector<std::vector<std::pair<int, double>>> neighbours(n);

    for (const auto& elem : pb.objective.elements)
        linear[static_cast<std::size_t>(elem.variable_index)] +=
          sign * elem.factor;

    for (const auto& elem : pb.objective.qelements) {
        const auto a = static_cast<std::size_t>(elem.variable_index_a);
        const auto b = static_cast<std::size_t>(elem.variable_index_b);

        if (a == b) {
            linear[a] += sign * elem.factor;
        } else {
            neighbours[a].emplace_back(elem.variable_index_b,
                                       sign * elem.factor);
            neighbours[b].emplace_back(elem.variable_index_a,
                                       sign * elem.factor);
        }
    }

    std::vector<std::int8_t> x(n);
    for (std::size_t i = 0; i != n; ++i)
        x[i] = linear[i] < 0.0;

    // Without quadratic terms the linear choice is optimal.
    if (pb.objective.qelements.empty())
        return x;

    // Otherwise, a tabu search flipping one variable at a time, the usual
    // heuristic for unconstrained binary quadratic problems: the best
    // non-tabu flip is applied even when it degrades, and a flipped variable
    // stays fixed for `tenure` iterations unless it improves the best value.
    std::vector<double> field = linear;
    for (std::size_t i = 0; i != n; ++i)
        if (x[i])
            for (const auto& [j, factor] : neighbours[i])
                field[static_cast<std::size_t>(j)] += factor;

    const auto delta = [&](std::size_t i) {
        return x[i] ? -field[i] : field[i];
    };

    const long long size = static_cast<long long>(n);
    const long long tenure = std::min<long long>(20, size / 4) + 1;
    const long long iterations =
      std::clamp<long long>(20'000'000 / size, 1'000, 100'000);

    std::vector<long long> tabu_until(n, 0);
    std::vector<std::int8_t> best = x;
    double current = 0.0;
    double best_value = 0.0;

    for (long long iter = 1; iter <= iterations; ++iter) {
        std::size_t chosen = n;
        double chosen_delta = 0.0;

        for (std::size_t i = 0; i != n; ++i) {
            const double d = delta(i);
            const bool allowed =
              tabu_until[i] < iter || current + d < best_value - 1e-9;

            if (allowed && (chosen == n || d < chosen_delta)) {
                chosen = i;
                chosen_delta = d;
            }
        }

        if (chosen == n)
            break;

        const double sign_flip = x[chosen] ? -1.0 : 1.0;
        x[chosen] = 1 - x[chosen];
        for (const auto& [j, factor] : neighbours[chosen])
            field[static_cast<std::size_t>(j)] += sign_flip * factor;

        current += chosen_delta;
        tabu_until[chosen] = iter + tenure;

        if (current < best_value - 1e-9) {
            best_value = current;
            best = x;
        }
    }

    return best;
}

solution
unconstrained_solution(const problem& pb)
{
    if (pb.vars.values.empty())
        return solution({}, pb.objective.value);

    auto x = unconstrained_assignment(pb);
    const double value = compute_solution(pb, x);

    return solution(std::move(x), value);
}

} // baryonyx namespace
