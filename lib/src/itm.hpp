/* Copyright (C) 2017 INRA
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublipnse, and/or sell copies of the Software, and to
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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_ITM_HPP

#include <baryonyx/core-compare>
#include <baryonyx/core-out>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <future>
#include <iterator>
#include <random>
#include <set>
#include <thread>
#include <unordered_map>

#include "fixed_array.hpp"
#include "matrix.hpp"
#include "private.hpp"
#include "utils.hpp"

#include <cassert>

namespace baryonyx {
namespace itm {

result
inequalities_1coeff_wedelin_solve(std::shared_ptr<context> ctx, problem& pb);

result
inequalities_1coeff_wedelin_optimize(std::shared_ptr<context> ctx,
                                     problem& pb,
                                     int thread);

enum class constraint_order
{
    none,
    reversing,
    random_sorting,
    infeasibility_decr,
    infeasibility_incr,
};

inline const char*
constraint_order_to_string(constraint_order type) noexcept
{
    static const char* ret[] = {
        "none",
        "reversing",
        "random-sorting",
        "infeasibility-decr",
        "infeasibility-incr",
    };

    switch (type) {
    case constraint_order::none:
        return ret[0];
    case constraint_order::reversing:
        return ret[1];
    case constraint_order::random_sorting:
        return ret[2];
    case constraint_order::infeasibility_decr:
        return ret[3];
    case constraint_order::infeasibility_incr:
        return ret[4];
    default:
        return ret[0];
    }
}

inline constraint_order
get_constraint_order(std::shared_ptr<context> ctx) noexcept
{
    auto str = ctx->get_string_parameter("constraint-order", "none");

    if (str == "none")
        return constraint_order::none;
    if (str == "reversing")
        return constraint_order::reversing;
    if (str == "random-sorting")
        return constraint_order::random_sorting;
    if (str == "infeasibility-decr")
        return constraint_order::infeasibility_decr;
    if (str == "infeasibility-incr")
        return constraint_order::infeasibility_incr;

    return constraint_order::none;
}

struct parameters
{
    parameters(std::shared_ptr<context> ctx)
      : time_limit(ctx->get_real_parameter("time-limit", -1.0))
      , theta(ctx->get_real_parameter("theta", 0.5))
      , delta(ctx->get_real_parameter("delta", 0.01))
      , kappa_min(ctx->get_real_parameter("kappa-min", 0.0))
      , kappa_step(ctx->get_real_parameter("kappa-step", 1.e-3))
      , kappa_max(ctx->get_real_parameter("kappa-max", 0.6))
      , alpha(ctx->get_real_parameter("alpha", 1.0))
      , pushing_k_factor(ctx->get_real_parameter("pushing-k-factor", 0.9))
      , pushes_limit(ctx->get_integer_parameter("pushes-limit", 10))
      , pushing_objective_amplifier(
          ctx->get_real_parameter("pushing-objective-amplifier", 5))
      , pushing_iteration_limit(
          ctx->get_integer_parameter("pushing-iteration-limit", 20))
      , limit(ctx->get_integer_parameter("limit", 1000))
      , w(ctx->get_integer_parameter("w", 20))
      , order(get_constraint_order(ctx))
      , preprocessing(ctx->get_string_parameter("preprocessing", "none"))
      , norm(ctx->get_string_parameter("norm", "none"))
      , serialize(ctx->get_integer_parameter("serialize", 0))
    {
        if (limit < 0)
            limit = std::numeric_limits<int>::max();

        ctx->info("solver: inequalities_1coeff_wedelin\n"
                  "solver parameters:\n"
                  "  - preprocessing: %s\n"
                  "  - constraint-order: %s\n"
                  "  - time_limit: %.10g\n"
                  "  - theta: %.10g\n"
                  "  - delta: %.10g\n"
                  "  - limit: %d\n"
                  "  - kappa: %.10g %.10g %.10g\n"
                  "  - alpha: %.10g\n"
                  "  - w: %d\n"
                  "  - norm: %s\n"
                  "  - serialise: %d\n"
                  "optimizer parameters:\n"
                  "  - pushed limit: %d\n"
                  "  - pushing objective amplifier: %.10g\n"
                  "  - pushing iteration limit: %d\n"
                  "  - pushing k factor: %.10g\n",
                  preprocessing.c_str(),
                  constraint_order_to_string(order),
                  time_limit,
                  theta,
                  delta,
                  limit,
                  kappa_min,
                  kappa_step,
                  kappa_max,
                  alpha,
                  w,
                  norm.c_str(),
                  serialize,
                  pushes_limit,
                  pushing_objective_amplifier,
                  pushing_iteration_limit,
                  pushing_k_factor);
    }

    double time_limit;
    double theta;
    double delta;
    double kappa_min;
    double kappa_step;
    double kappa_max;
    double alpha;
    double pushing_k_factor;
    int pushes_limit;
    double pushing_objective_amplifier;
    int pushing_iteration_limit;
    int limit;
    int w;
    constraint_order order;
    std::string preprocessing;
    std::string norm;
    int serialize;
};

struct merged_constraint
{
    merged_constraint(const std::vector<function_element>& elements_,
                      int min_,
                      int max_)
      : elements(elements_)
      , min(min_)
      , max(max_)
    {
    }

    std::vector<function_element> elements;
    int min;
    int max;
};

struct merged_constraint_hash
{
    inline size_t operator()(const std::vector<function_element>& fct) const
      noexcept
    {
        std::size_t seed{ fct.size() };

        for (auto& elem : fct)
            seed ^=
              elem.variable_index + 0x9e3779b9 + (seed << 6) + (seed >> 2);

        return seed;
    }
};

inline void
sort_merged_constraints(std::shared_ptr<context> ctx,
                        const problem& pb,
                        const parameters& params,
                        std::vector<merged_constraint>& ret)
{
    ctx->info("  - static sort merged constraints: %s\n",
              params.preprocessing.c_str());

    if (params.preprocessing.empty() or params.preprocessing == "none")
        return;

    std::vector<int> vars(pb.vars.values.size(), 0);
    std::vector<std::set<int>> linkvars(pb.vars.values.size());
    std::vector<std::set<int>> linkcst(pb.vars.values.size());

    for (auto it{ ret.begin() }, et{ ret.end() }; it != et; ++it) {
        for (std::size_t i{ 0 }, e{ it->elements.size() }; i != e; ++i) {

            linkcst[it->elements[i].variable_index].emplace(
              std::distance(it, ret.end()));

            for (std::size_t j{ 0 }; j != e; ++j) {
                if (i != j)
                    linkvars[it->elements[i].variable_index].emplace(
                      it->elements[j].variable_index);
            }
        }

        for (auto& f : it->elements) {
            vars[f.variable_index]++;
        }
    }

    std::vector<std::pair<merged_constraint, int>> tosort;

    if (params.preprocessing == "variables-number") {
        for (auto& cst : ret) {
            tosort.emplace_back(cst, 0);
            for (auto& elem : cst.elements) {
                for (auto& s : linkcst[elem.variable_index])
                    tosort.back().second += linkvars[s].size();
            }
        }

        std::sort(
          tosort.begin(), tosort.end(), [](const auto& lhs, const auto& rhs) {
              return rhs.second < lhs.second;
          });

        ret.clear();
        std::transform(tosort.begin(),
                       tosort.end(),
                       std::back_inserter(ret),
                       [](const auto& elem) { return elem.first; });

        return;
    }

    if (params.preprocessing == "variables-weight") {
        std::sort(
          ret.begin(), ret.end(), [vars](const auto& lhs, const auto& rhs) {
              int sumlhs{ 0 };
              int sumrhs{ 0 };

              for (auto& f : lhs.elements)
                  sumlhs += vars[f.variable_index];

              for (auto& f : rhs.elements)
                  sumrhs += vars[f.variable_index];

              return sumlhs < sumrhs;
          });

        return;
    }

    if (params.preprocessing == "constraints-weight") {
        std::sort(ret.begin(),
                  ret.end(),
                  [linkvars, linkcst](const auto& lhs, const auto& rhs) {
                      int sumlhs{ 1 };
                      int sumrhs{ 1 };

                      for (auto& f : lhs.elements)
                          sumlhs *= linkcst[f.variable_index].size();

                      for (auto& f : rhs.elements)
                          sumrhs *= linkcst[f.variable_index].size();

                      return sumrhs > sumlhs;
                  });

        return;
    }

    if (params.preprocessing == "implied") {
        // Algorithm to build a tosort vector according to constraints of type:
        // -x1 -x2 -x3 +x4 <= 0

        for (auto& cst : ret) {
            tosort.emplace_back(cst, 0);

            int nbneg{ 0 }, nbpos{ 0 };

            for (auto& elem : cst.elements)
                if (elem.factor < 0)
                    nbneg++;
                else
                    nbpos++;

            if (((nbneg > 1 and nbpos == 1) or (nbpos > 1 and nbneg == 1)) and
                cst.min == cst.max) {
                tosort.back().second = 1;
            }
        }

        std::sort(
          tosort.begin(), tosort.end(), [](const auto& lhs, const auto& rhs) {
              return rhs.second < lhs.second;
          });

        ret.clear();
        std::transform(tosort.begin(),
                       tosort.end(),
                       std::back_inserter(ret),
                       [](const auto& elem) { return elem.first; });

        return;
    }
}

inline std::vector<merged_constraint>
make_merged_constraints(std::shared_ptr<context> ctx,
                        const problem& pb,
                        const parameters& p)
{
    std::unordered_map<std::vector<function_element>,
                       std::size_t,
                       merged_constraint_hash>
      cache;

    //
    // Fills the cache and uses it to search same function and merge less and
    // greater equal constraints.
    //

    std::vector<merged_constraint> ret;

    for (const auto& elem : pb.equal_constraints) {
        cache.emplace(elem.elements, ret.size());
        ret.emplace_back(elem.elements,
                         numeric_cast<int>(std::lround(elem.value)),
                         numeric_cast<int>(std::lround(elem.value)));
    }

    for (const auto& elem : pb.less_constraints) {
        auto it = cache.find(elem.elements);
        if (it == cache.end()) {
            cache.emplace(elem.elements, ret.size());
            ret.emplace_back(elem.elements,
                             std::numeric_limits<int>::min(),
                             numeric_cast<int>(std::lround(elem.value)));
        } else {
            ret[it->second].max = std::min(
              ret[it->second].max, numeric_cast<int>(std::lround(elem.value)));
        }
    }

    for (const auto& elem : pb.greater_constraints) {
        auto it = cache.find(elem.elements);
        if (it == cache.end()) {
            cache.emplace(elem.elements, ret.size());
            ret.emplace_back(
              elem.elements, elem.value, std::numeric_limits<int>::max());
        } else {
            ret[it->second].min = std::max(
              ret[it->second].min, numeric_cast<int>(std::lround(elem.value)));
        }
    }

    ctx->info("  - merged constraints: %d\n",
              numeric_cast<int>(pb.equal_constraints.size() +
                                pb.less_constraints.size() +
                                pb.greater_constraints.size() - ret.size()));

    sort_merged_constraints(ctx, pb, p, ret);

    return ret;
}
}
}

#endif
