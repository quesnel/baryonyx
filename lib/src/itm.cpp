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

#include "itm.hpp"
#include "itm-common.hpp"

#include "itm-solver-equalities-01.hpp"
#include "itm-solver-inequalities-Z.hpp"

#include <set>

#include <cassert>

namespace baryonyx {
namespace itm {

struct merged_constraint_hash
{
    inline size_t operator()(
      const std::vector<baryonyx::function_element>& fct) const noexcept
    {
        std::size_t seed{ fct.size() };

        for (auto& elem : fct)
            seed ^=
              elem.variable_index + 0x9e3779b9 + (seed << 6) + (seed >> 2);

        return seed;
    }
};

//
// Fills the cache and uses it to search same function and merge equal, less
// and greater equal constraints. To be use in any order, the fill operation
// depends of the operator type passed in parameter.
//
template<typename cacheT, typename retT>
static void
fill_merged_constraints(const baryonyx::context_ptr& ctx,
                        cacheT& cache,
                        baryonyx::operator_type op,
                        const baryonyx::problem& pb,
                        retT& ret)
{
    switch (op) {
    case baryonyx::operator_type::equal:
        for (const auto& elem : pb.equal_constraints) {
            auto it = cache.find(elem.elements);
            if (it == cache.end()) {
                cache.emplace(elem.elements, ret.size());
                ret.emplace_back(elem.elements,
                                 baryonyx::numeric_cast<int>(elem.value),
                                 baryonyx::numeric_cast<int>(elem.value),
                                 elem.id);
            } else {
                if (ret[it->second].min <= elem.value and
                    elem.value <= ret[it->second].max) {
                    ret[it->second].min = elem.value;
                    ret[it->second].max = elem.value;
                } else {
                    error(ctx,
                          "Constraint {}: is inconsistent with previous {}.\n",
                          elem.id,
                          ret[it->second].id);
                }
            }
        }
        break;
    case baryonyx::operator_type::less:
        for (const auto& elem : pb.less_constraints) {
            auto it = cache.find(elem.elements);
            if (it == cache.end()) {
                cache.emplace(elem.elements, ret.size());
                ret.emplace_back(elem.elements,
                                 std::numeric_limits<int>::min(),
                                 baryonyx::numeric_cast<int>(elem.value),
                                 elem.id);
            } else {
                ret[it->second].max =
                  std::min(ret[it->second].max,
                           baryonyx::numeric_cast<int>(elem.value));
            }
        }
        break;
    case baryonyx::operator_type::greater:
        for (const auto& elem : pb.greater_constraints) {
            auto it = cache.find(elem.elements);
            if (it == cache.end()) {
                cache.emplace(elem.elements, ret.size());
                ret.emplace_back(elem.elements,
                                 elem.value,
                                 std::numeric_limits<int>::max(),
                                 elem.id);
            } else {
                ret[it->second].min =
                  std::max(ret[it->second].min,
                           baryonyx::numeric_cast<int>(elem.value));
            }
        }
        break;
    default:
        break;
    }
}

static inline baryonyx::operator_type
operator_from_string(const std::string& op) noexcept
{
    if (op == "equal")
        return baryonyx::operator_type::equal;
    if (op == "less")
        return baryonyx::operator_type::less;
    if (op == "greater")
        return baryonyx::operator_type::greater;

    return baryonyx::operator_type::undefined;
}

//
// Restore original order: reorder merged constraints according to the
// merged_constraint::id (initially initialized in the lp format reading
// process).
//
static std::vector<baryonyx::itm::merged_constraint>
make_unsorted_merged_constraints(const baryonyx::context_ptr& ctx,
                                 const baryonyx::problem& pb)
{
    info(ctx, "  - merge constraints without any sort\n");

    std::unordered_map<std::vector<baryonyx::function_element>,
                       std::size_t,
                       merged_constraint_hash>
      cache;

    std::vector<baryonyx::itm::merged_constraint> ret;
    ret.reserve(pb.equal_constraints.size() + pb.less_constraints.size() +
                pb.greater_constraints.size());

    fill_merged_constraints(
      ctx, cache, baryonyx::operator_type::equal, pb, ret);
    fill_merged_constraints(
      ctx, cache, baryonyx::operator_type::less, pb, ret);
    fill_merged_constraints(
      ctx, cache, baryonyx::operator_type::greater, pb, ret);

    std::sort(ret.begin(), ret.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.id < rhs.id;
    });

    return ret;
}

//
// Reorder constraints according to the type of operator proposed into the
// preprocessing parameter. If parameter is badly defined, we fall back to the
// unsorted_merged_constraints.
//
static std::vector<baryonyx::itm::merged_constraint>
make_ordered_merged_constraints(const baryonyx::context_ptr& ctx,
                                const baryonyx::problem& pb,
                                const baryonyx::itm::parameters& p)
{
    info(ctx, "  - merge constraints with ordered sort\n");

    auto first = p.preprocessing.find(',', 0);
    if (first != std::string::npos) {
        auto second = p.preprocessing.find(',', first + 1);
        if (second != std::string::npos) {
            auto op0 = p.preprocessing.substr(0, first);
            auto op1 = p.preprocessing.substr(first + 1, second - first - 1);
            auto op2 = p.preprocessing.substr(second + 1);

            baryonyx::operator_type op[3] = { baryonyx::operator_type::less,
                                              baryonyx::operator_type::greater,
                                              baryonyx::operator_type::equal };

            op[0] = operator_from_string(op0);
            op[1] = operator_from_string(op1);
            op[2] = operator_from_string(op2);

            if (op[0] != baryonyx::operator_type::undefined and
                op[1] != baryonyx::operator_type::undefined and
                op[2] != baryonyx::operator_type::undefined and
                op[0] != op[1] and op[1] != op[2] and op[0] != op[2]) {

                std::unordered_map<std::vector<baryonyx::function_element>,
                                   std::size_t,
                                   merged_constraint_hash>
                  cache;

                std::vector<baryonyx::itm::merged_constraint> ret;
                ret.reserve(pb.equal_constraints.size() +
                            pb.less_constraints.size() +
                            pb.greater_constraints.size());

                fill_merged_constraints(ctx, cache, op[0], pb, ret);
                fill_merged_constraints(ctx, cache, op[1], pb, ret);
                fill_merged_constraints(ctx, cache, op[2], pb, ret);

                return ret;
            }
        }

        warning(ctx,
                "    Bad preprocessing order mode: `{}`. Falling back "
                "to the unsorted merged constraint\n",
                p.preprocessing);
    }

    return make_unsorted_merged_constraints(ctx, pb);
}

//
// Reorder constraints according to the special operator proposed into the
// preprocessing parameter. If parameter is badly defined, we fall back to the
// unsorted_merged_constraints.
//
static std::vector<baryonyx::itm::merged_constraint>
make_special_merged_constraints(const baryonyx::context_ptr& ctx,
                                const baryonyx::problem& pb,
                                const baryonyx::itm::parameters& p)
{
    info(ctx,
         "  - merge constraint according to the `{}` algorithm\n",
         p.preprocessing);

    auto ret = make_unsorted_merged_constraints(ctx, pb);

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

    std::vector<std::pair<baryonyx::itm::merged_constraint, int>> tosort;

    if (p.preprocessing == "variables-number") {
        for (auto& cst : ret) {
            tosort.emplace_back(cst, 0);
            for (auto& elem : cst.elements) {
                for (auto& s : linkcst[elem.variable_index])
                    tosort.back().second +=
                      static_cast<int>(linkvars[s].size());
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

        return ret;
    }

    if (p.preprocessing == "variables-weight") {
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

        return ret;
    }

    if (p.preprocessing == "constraints-weight") {
        std::sort(ret.begin(),
                  ret.end(),
                  [linkvars, linkcst](const auto& lhs, const auto& rhs) {
                      std::size_t sumlhs{ 1 };
                      std::size_t sumrhs{ 1 };

                      for (auto& f : lhs.elements)
                          sumlhs *= linkcst[f.variable_index].size();

                      for (auto& f : rhs.elements)
                          sumrhs *= linkcst[f.variable_index].size();

                      return sumrhs > sumlhs;
                  });

        return ret;
    }

    if (p.preprocessing == "implied") {
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

        return ret;
    }

    return ret;
}

std::vector<merged_constraint>
make_merged_constraints(const context_ptr& ctx,
                        const problem& pb,
                        const parameters& p)
{
    auto original_nb = static_cast<int>(pb.equal_constraints.size() +
                                        pb.less_constraints.size() +
                                        pb.greater_constraints.size());

    std::vector<merged_constraint> ret;
    ret.reserve(original_nb);

    if (p.preprocessing.empty() or p.preprocessing == "none") {
        ret = make_unsorted_merged_constraints(ctx, pb);
    } else if (p.preprocessing.find(',') != std::string::npos) {
        ret = make_ordered_merged_constraints(ctx, pb, p);
    } else {
        ret = make_special_merged_constraints(ctx, pb, p);
    }

    info(ctx,
         "  - merged constraints removed: {}\n",
         original_nb - static_cast<int>(ret.size()));

    return ret;
}

template<typename floatingpointT, typename modeT, typename randomT>
struct solver_inequalities_Zcoeff;

template<typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT>
inline result
dispatch_solver(const context_ptr& ctx, problem& pb, const itm::parameters& p)
{
    switch (pb.problem_type) {
    case problem_solver_type::equalities_01:
        return solve_problem<
          solver_equalities_01coeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb, p);
    case problem_solver_type::equalities_101:
        return solve_problem<
          solver_inequalities_Zcoeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb, p);
    case problem_solver_type::equalities_Z:
        return solve_problem<
          solver_inequalities_Zcoeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb, p);
    case problem_solver_type::inequalities_01:
        return solve_problem<
          solver_inequalities_Zcoeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb, p);
    case problem_solver_type::inequalities_101:
        return solve_problem<
          solver_inequalities_Zcoeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb, p);
    case problem_solver_type::inequalities_Z:
        return solve_problem<
          solver_inequalities_Zcoeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb, p);
    }

    return result(result_status::internal_error);
}

template<typename floatingpointT,
         typename modeT,
         typename constraintOrderT,
         typename randomT>
inline result
dispatch_optimizer(const context_ptr& ctx,
                   problem& pb,
                   const itm::parameters& p,
                   int thread)
{
    switch (pb.problem_type) {
    case problem_solver_type::equalities_01:
        return optimize_problem<
          solver_equalities_01coeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb, p, thread);
    case problem_solver_type::equalities_101:
        return optimize_problem<
          solver_inequalities_Zcoeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb, p, thread);
    case problem_solver_type::equalities_Z:
        return optimize_problem<
          solver_inequalities_Zcoeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb, p, thread);
    case problem_solver_type::inequalities_01:
        return optimize_problem<
          solver_inequalities_Zcoeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb, p, thread);
    case problem_solver_type::inequalities_101:
        return optimize_problem<
          solver_inequalities_Zcoeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb, p, thread);
    case problem_solver_type::inequalities_Z:
        return optimize_problem<
          solver_inequalities_Zcoeff<floatingpointT, modeT, randomT>,
          floatingpointT,
          modeT,
          constraintOrderT,
          randomT>(ctx, pb, p, thread);
    }

    return result(result_status::internal_error);
}

template<typename realT, typename modeT, typename randomT>
inline result
dispatch_solver_parameters(const context_ptr& ctx,
                           problem& pb,
                           const itm::parameters& p)
{
    switch (p.order) {
    case itm::constraint_order::none:
        return dispatch_solver<realT,
                               modeT,
                               compute_none<realT, randomT>,
                               randomT>(ctx, pb, p);
    case itm::constraint_order::reversing:
        return dispatch_solver<realT,
                               modeT,
                               compute_reversing<realT, randomT>,
                               randomT>(ctx, pb, p);
    case itm::constraint_order::random_sorting:
        return dispatch_solver<realT,
                               modeT,
                               compute_random<realT, randomT>,
                               randomT>(ctx, pb, p);
    case itm::constraint_order::infeasibility_decr:
        return dispatch_solver<
          realT,
          modeT,
          compute_infeasibility<realT, randomT, compute_infeasibility_decr>,
          randomT>(ctx, pb, p);
    case itm::constraint_order::infeasibility_incr:
        return dispatch_solver<
          realT,
          modeT,
          compute_infeasibility<realT, randomT, compute_infeasibility_incr>,
          randomT>(ctx, pb, p);
    }

    return result(result_status::internal_error);
}

template<typename realT, typename modeT, typename randomT>
inline result
dispatch_optimizer_parameters(const context_ptr& ctx,
                              problem& pb,
                              const itm::parameters& p,
                              int thread)
{
    switch (p.order) {
    case itm::constraint_order::none:
        return dispatch_optimizer<realT,
                                  modeT,
                                  compute_none<realT, randomT>,
                                  randomT>(ctx, pb, p, thread);
    case itm::constraint_order::reversing:
        return dispatch_optimizer<realT,
                                  modeT,
                                  compute_reversing<realT, randomT>,
                                  randomT>(ctx, pb, p, thread);
    case itm::constraint_order::random_sorting:
        return dispatch_optimizer<realT,
                                  modeT,
                                  compute_random<realT, randomT>,
                                  randomT>(ctx, pb, p, thread);
    case itm::constraint_order::infeasibility_decr:
        return dispatch_optimizer<
          realT,
          modeT,
          compute_infeasibility<realT, randomT, compute_infeasibility_decr>,
          randomT>(ctx, pb, p, thread);
    case itm::constraint_order::infeasibility_incr:
        return dispatch_optimizer<

          realT,
          modeT,
          compute_infeasibility<realT, randomT, compute_infeasibility_incr>,
          randomT>(ctx, pb, p, thread);
    }

    return result(result_status::internal_error);
}

result
solve(const baryonyx::context_ptr& ctx, problem& pb)
{
    info(ctx, "Solve mode\n");
    parameters p(ctx);

    using random_type = std::default_random_engine;

    if (pb.type == baryonyx::objective_function_type::maximize) {
        switch (p.float_type) {
        case floating_point_type::float_type:
            return dispatch_solver_parameters<float,
                                              maximize_tag,
                                              random_type>(ctx, pb, p);
        case floating_point_type::double_type:
            return dispatch_solver_parameters<double,
                                              maximize_tag,
                                              random_type>(ctx, pb, p);
        case floating_point_type::longdouble_type:
            return dispatch_solver_parameters<long double,
                                              maximize_tag,
                                              random_type>(ctx, pb, p);
        }
    } else {
        switch (p.float_type) {
        case floating_point_type::float_type:
            return dispatch_solver_parameters<float,
                                              minimize_tag,
                                              random_type>(ctx, pb, p);
        case floating_point_type::double_type:
            return dispatch_solver_parameters<double,
                                              minimize_tag,
                                              random_type>(ctx, pb, p);
        case floating_point_type::longdouble_type:
            return dispatch_solver_parameters<long double,
                                              minimize_tag,
                                              random_type>(ctx, pb, p);
        }
    }

    return result(result_status::internal_error);
}

result
optimize(const baryonyx::context_ptr& ctx, problem& pb, int thread)
{
    info(ctx, "Optimize mode\n");
    parameters p(ctx);

    using random_type = std::default_random_engine;

    if (pb.type == baryonyx::objective_function_type::maximize) {
        switch (p.float_type) {
        case floating_point_type::float_type:
            return dispatch_optimizer_parameters<float,
                                                 maximize_tag,
                                                 random_type>(
              ctx, pb, p, thread);
        case floating_point_type::double_type:
            return dispatch_optimizer_parameters<double,
                                                 maximize_tag,
                                                 random_type>(
              ctx, pb, p, thread);
        case floating_point_type::longdouble_type:
            return dispatch_optimizer_parameters<long double,
                                                 maximize_tag,
                                                 random_type>(
              ctx, pb, p, thread);
        }
    } else {
        switch (p.float_type) {
        case floating_point_type::float_type:
            return dispatch_optimizer_parameters<float,
                                                 minimize_tag,
                                                 random_type>(
              ctx, pb, p, thread);
        case floating_point_type::double_type:
            return dispatch_optimizer_parameters<double,
                                                 minimize_tag,
                                                 random_type>(
              ctx, pb, p, thread);
        case floating_point_type::longdouble_type:
            return dispatch_optimizer_parameters<long double,
                                                 minimize_tag,
                                                 random_type>(
              ctx, pb, p, thread);
        }
    }

    return result(result_status::internal_error);
}

} // namespace itm
} // namespace baryonyx
