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

#ifndef ORG_VLEPROJECT_LP_INEQUALITIES_1COEFF_HPP
#define ORG_VLEPROJECT_LP_INEQUALITIES_1COEFF_HPP

#include "lpcore-compare"
#include "lpcore-out"
#include "matrix.hpp"
#include "mitm.hpp"
#include "utils.hpp"

#include <Eigen/Core>

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

namespace lp {
namespace inequalities_1coeff {

using AP_type = lp::SparseArray<std::int8_t, double>;
using b_type = Eigen::Matrix<int, 2, Eigen::Dynamic>;
using c_type = Eigen::VectorXf;
using x_type = Eigen::VectorXi;
using pi_type = Eigen::VectorXf;

enum class constraint_order
{
    none,
    reversing,
    random_sorting,
    infeasibility_decr,
    infeasibility_incr,
};

const char*
constraint_order_to_string(constraint_order type)
{
    static const char* ret[] = {
        "none",
        "reversing",
        "random-sorting",
        "infeasibility-decr",
        "infeasibility-incr",
    };

    switch (type) {
    case inequalities_1coeff::constraint_order::none:
        return ret[0];
    case inequalities_1coeff::constraint_order::reversing:
        return ret[1];
    case inequalities_1coeff::constraint_order::random_sorting:
        return ret[2];
    case inequalities_1coeff::constraint_order::infeasibility_decr:
        return ret[3];
    case inequalities_1coeff::constraint_order::infeasibility_incr:
        return ret[4];
    }

    return nullptr;
}

double
get_real(std::shared_ptr<context> ctx,
         const std::map<std::string, parameter>& params,
         std::string param,
         double def)
{
    auto it = params.find(param);
    if (it == params.cend())
        return def;

    if (it->second.type == parameter::tag::real)
        return it->second.d;

    if (it->second.type == parameter::tag::integer)
        return static_cast<double>(it->second.l);

    ctx->warning("fail to convert parameter %s\n", param.c_str());

    return def;
}

long int
get_integer(std::shared_ptr<context> ctx,
            const std::map<std::string, parameter>& params,
            std::string param,
            long int def)
{
    auto it = params.find(param);
    if (it == params.cend())
        return def;

    if (it->second.type == parameter::tag::integer)
        return it->second.l;

    if (it->second.type == parameter::tag::real)
        return static_cast<long>(it->second.d);

    ctx->warning("fail to convert parameter %s\n", param.c_str());

    return def;
}

inequalities_1coeff::constraint_order
get_constraint_order(std::shared_ptr<context> ctx,
                     const std::map<std::string, parameter>& params,
                     std::string param,
                     inequalities_1coeff::constraint_order def)
{
    auto it = params.find(param);
    if (it == params.cend())
        return def;

    if (it->second.type != parameter::tag::string)
        return def;

    if (it->second.s == "none")
        return inequalities_1coeff::constraint_order::none;
    if (it->second.s == "reversing")
        return inequalities_1coeff::constraint_order::reversing;
    if (it->second.s == "random-sorting")
        return inequalities_1coeff::constraint_order::random_sorting;
    if (it->second.s == "infeasibility-decr")
        return inequalities_1coeff::constraint_order::infeasibility_decr;
    if (it->second.s == "infeasibility-incr")
        return inequalities_1coeff::constraint_order::infeasibility_incr;

    ctx->warning("fail to convert parameter %s\n", param.c_str());

    return def;
}

std::string
get_pre_constraint_order(std::shared_ptr<context> /*ctx*/,
                         const std::map<std::string, parameter>& params,
                         std::string param)
{
    auto it = params.find(param);
    if (it == params.cend())
        return {};

    return it->second.s;
}

struct parameters
{
    parameters(std::shared_ptr<context> ctx,
               const std::map<std::string, parameter>& params)
      : time_limit(get_real(ctx, params, "time-limit", -1.0))
      , theta(get_real(ctx, params, "theta", 0.5))
      , delta(get_real(ctx, params, "delta", 0.01))
      , kappa_min(get_real(ctx, params, "kappa-min", 0.0))
      , kappa_step(get_real(ctx, params, "kappa-step", 1.e-3))
      , kappa_max(get_real(ctx, params, "kappa-max", 0.6))
      , alpha(get_real(ctx, params, "alpha", 1.0))
      , pushing_k_factor(get_real(ctx, params, "pushing-k-factor", 0.9))
      , pushes_limit(get_integer(ctx, params, "pushes-limit", 10l))
      , pushing_objective_amplifier(
          get_integer(ctx, params, "pushing-objective-amplifier", 5l))
      , pushing_iteration_limit(
          get_integer(ctx, params, "pushing-iteration-limit", 20l))
      , limit(get_integer(ctx, params, "limit", 1000l))
      , w(get_integer(ctx, params, "w", 20l))
      , order(get_constraint_order(ctx,
                                   params,
                                   "constraint-order",
                                   constraint_order::none))
      , preprocessing(get_pre_constraint_order(ctx, params, "preprocessing"))
      , serialize(get_integer(ctx, params, "serialize", 0l))
    {
        ctx->info("solver: inequalities_1coeff_wedelin\n"
                  "solver parameters:\n"
                  "  - preprocessing: %s\n"
                  "  - constraint-order: %s\n"
                  "  - time_limit: %.10g\n"
                  "  - theta: %.10g\n"
                  "  - delta: %.10g\n"
                  "  - limit: %ld\n"
                  "  - kappa: %.10g %.10g %.10g\n"
                  "  - alpha: %.10g\n"
                  "  - w: %ld\n"
                  "  - serialise: %d\n"
                  "optimizer parameters:\n"
                  "  - pushed limit: %ld\n"
                  "  - pushing objective amplifier: %ld\n"
                  "  - pushing iteration limit: %ld\n"
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
    long int pushes_limit;
    long int pushing_objective_amplifier;
    long int pushing_iteration_limit;
    long int limit;
    long int w;
    constraint_order order;
    std::string preprocessing;
    bool serialize;
};

struct maximize_tag
{
};

struct minimize_tag
{
};

template <typename iteratorT, typename randomT>
void
random_shuffle_unique(iteratorT begin, iteratorT end, randomT& rng)
{
    auto ret = begin++;
    for (; begin != end; ++begin) {
        if (ret->value != begin->value) {
            std::shuffle(ret, begin, rng);
            ret = begin;
        }
    }

    std::shuffle(ret, begin, rng);
}

template <typename iteratorT, typename randomT>
void
calculator_sort(iteratorT begin, iteratorT end, randomT& rng, minimize_tag)
{
    if (std::distance(begin, end) > 1) {

        //
        // TODO We can use stable_sort to be sure reproducible run:
        // std::stable_sort(begin, end, [](const auto& lhs, const auto& rhs)
        //

        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return lhs.value < rhs.value;
        });

        random_shuffle_unique(begin, end, rng);
    }
}

template <typename iteratorT, typename randomT>
void
calculator_sort(iteratorT begin, iteratorT end, randomT& rng, maximize_tag)
{
    if (std::distance(begin, end) > 1) {

        //
        // TODO We can use stable_sort to be sure reproducible run:
        // std::stable_sort(begin, end, [](const auto& lhs, const auto& rhs)
        //

        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return rhs.value < lhs.value;
        });

        random_shuffle_unique(begin, end, rng);
    }
}

bool
stop_iterating(double value, minimize_tag) noexcept
{
    return value > 0;
}

bool
stop_iterating(double value, maximize_tag) noexcept
{
    return value < 0;
}

bool
is_better_solution(double lhs, double rhs, minimize_tag) noexcept
{
    return lhs < rhs;
}

bool
init_x(int v, minimize_tag) noexcept
{
    return v <= 0;
}

bool
init_x(int v, maximize_tag) noexcept
{
    return v >= 0;
}

bool
is_better_solution(double lhs, double rhs, maximize_tag) noexcept
{
    return lhs > rhs;
}

bool
is_time_limit(double limit,
              std::chrono::steady_clock::time_point begin,
              std::chrono::steady_clock::time_point end) noexcept
{
    if (limit <= 0)
        return false;

    return std::chrono::duration_cast<std::chrono::duration<double>>(end -
                                                                     begin)
             .count() > limit;
}

template <typename modeT, typename randomT>
struct constraint_calculator
{
    using mode_type = modeT;
    using random_generator_type = randomT;

    struct r_data
    {
        r_data(double value_, index index_)
          : value(value_)
          , id(index_)
        {
        }

        double value;
        index id;
    };

    random_generator_type& rng;
    AP_type& ap;
    b_type& b;
    const c_type& cost;
    x_type& x;
    pi_type& pi;
    std::vector<r_data> r;
    std::vector<index> C; // Stores variables with negative coefficient.
    index m;
    index n;

    constraint_calculator(random_generator_type& rng_,
                          index k,
                          index m_,
                          index n_,
                          AP_type& ap_,
                          b_type& b_,
                          const c_type& c_,
                          x_type& x_,
                          pi_type& pi_)
      : rng(rng_)
      , ap(ap_)
      , b(b_)
      , cost(c_)
      , x(x_)
      , pi(pi_)
      , m(m_)
      , n(n_)
    {
        const auto& ak{ ap.row(k) };
        const auto& va{ ap.A() };

        r.reserve(ak.size());

        for (std::size_t i{ 0 }, e{ ak.size() }; i != e; ++i) {

            //
            // We don't need to represent the I vector since the SparseArray
            // stores non null coefficient in each row (and column).
            // I.emplace_back(ak[i].position);
            //

            if (va[ak[i].value] != 0)
                r.emplace_back(0, ak[i].position);

            if (va[ak[i].value] < 0)
                C.emplace_back(ak[i].position);
        }
    }

    void update_row(index k, double kappa, double delta, double theta)
    {
        std::uniform_real_distribution<> dst(0.0, 1e-4);
        const auto& ak{ ap.row(k) };
        const auto& va{ ap.A() };

        //
        // Decrease influence of local preferences. 0 will completely reset the
        // preference values for the current row. > 0 will keep former decision
        // in mind.
        //

        ap.mult_row_p(k, theta);

        //
        // Calculate reduced costs
        //

        for (std::size_t i{ 0 }, e{ ak.size() }; i != e; ++i) {
            double sum_a_pi{ 0 };
            double sum_a_p{ 0 };

            const auto& H{ ap.column(ak[i].position) };

            for (std::size_t h{ 0 }, eh{ H.size() }; h != eh; ++h) {
                sum_a_pi +=
                  ap.A(H[h].position, ak[i].position) * pi(H[h].position);
                sum_a_p += ap.A(H[h].position, ak[i].position) *
                           ap.P(H[h].position, ak[i].position);
            }

            r[i].id = ak[i].position;
            r[i].value = cost(ak[i].position) - sum_a_pi - sum_a_p;

            if (is_essentially_equal(r[i].value, 0.0, 1e-10))
                r[i].value += (std::signbit(r[i].value) ? -1. : 1.) * dst(rng);
        }

        //
        // Negate reduced costs and coefficients of these variables. We need to
        // parse the row Ak[i] because we need to use r[i] not available in C.
        //

        if (not C.empty()) {
            for (std::size_t i{ 0 }, e{ ak.size() }; i != e; ++i) {
                if (va[ak[i].value] < 0) {
                    r[i].value = -r[i].value;
                    ap.invert(k, ak[i].position);
                }
            }

            b(0, k) += C.size();
            b(1, k) += C.size();
        }

        calculator_sort(r.begin(), r.end(), rng, mode_type());

        //
        // The bkmin and bkmax constraint bounds are not equal and can be
        // assigned to -infinity or +infinity. We have to scan the r vector and
        // search a value j such as b(0, k) <= Sum A(k, r[j]) < b(1, k).
        //

        index i{ 0 }, selected{ -1 }, first, second;
        const index endi{ numeric_cast<index>(r.size()) };
        int sum{ 0 };

        for (; i != endi; ++i) {
            sum += ap.A(k, r[i].id);

            if (b(0, k) <= sum)
                break;
        }

        //
        // If the b(0, k) can not be reached, this is an error of the
        // preprocessing step.
        //

        if (b(0, k) > sum)
            throw solver_error(solver_error::tag::unrealisable_constraint);

        //
        // If all variable must be assigned to 0, we let selected assigned to 0
        // and we go to the next part otherwise, we continue to scan.
        //

        if (b(0, k) <= sum and sum <= b(1, k)) {
            selected = i;
            for (; i != endi; ++i) {
                sum += ap.A(k, r[i].id);

                if (sum <= b(1, k)) {
                    if (stop_iterating(r[i].value, mode_type()))
                        break;
                    ++selected;
                } else
                    break;
            }

            if (i == endi)
                throw solver_error(solver_error::tag::unrealisable_constraint);

            first = selected;
            second = selected + 1;
        }

        if (selected < 0) {
            for (index j{ 0 }; j < endi; ++j) {
                x(r[j].id) = 0;
                ap.add_p(k, r[j].id, -delta);
            }
        } else if (second >= endi) {
            for (index j{ 0 }; j < endi; ++j) {
                x(r[j].id) = 1;
                ap.add_p(k, r[j].id, +delta);
            }
        } else {
            pi(k) += ((r[first].value + r[second].value) / 2.0);

            double d = delta + ((kappa / (1.0 - kappa)) *
                                (r[second].value - r[first].value));

            index j{ 0 };
            for (; j <= selected; ++j) {
                x(r[j].id) = 1;
                ap.add_p(k, r[j].id, +d);
            }

            for (; j != endi; ++j) {
                x(r[j].id) = 0;
                ap.add_p(k, r[j].id, -d);
            }
        }

        //
        // Clean up: correct negated costs and adjust value of negated
        // variables.
        //

        if (not C.empty()) {
            b(0, k) -= C.size();
            b(1, k) -= C.size();

            for (std::size_t i{ 0 }, e{ C.size() }; i != e; ++i) {
                ap.invert(k, C[i]);
                x[C[i]] = 1 - x[C[i]];
            }
        }
    }
};

c_type
make_objective_function(const lp::objective_function& obj, index n)
{
    c_type ret = c_type::Zero(n);

    for (const auto& elem : obj.elements)
        ret(elem.variable_index) += elem.factor;

    return ret;
}

struct merged_constraint
{
    merged_constraint(const std::vector<lp::function_element>& elements_,
                      int min_,
                      int max_)
      : elements(elements_)
      , min(min_)
      , max(max_)
    {
    }

    std::vector<lp::function_element> elements;
    int min;
    int max;
};

struct merged_constraint_hash
{
    inline size_t operator()(
      const std::vector<lp::function_element>& fct) const noexcept
    {
        std::size_t seed{ fct.size() };

        for (auto& elem : fct)
            seed ^=
              elem.variable_index + 0x9e3779b9 + (seed << 6) + (seed >> 2);

        return seed;
    }
};

/**
 * If a variable is assigned, we remove then from the list of variable, we
 * remove the constraint and we remove all reference to this variable in all
 * constraints.
 */
template <typename constraintsT>
long
remove_small_constraints(constraintsT& csts) noexcept
{
    for (auto& elem : csts) {
        if (elem.elements.size() == 1) {
        }
    }

    return 0;
}

/**
 * For each constraints, this function removes element from constraint
 * functions where the factor equal 0. Constraints are remove for all
 * empty.
 *
 * \param csts The constraints to process.
 * \return a \c std::tuple of two integers. First is the number of 0 factor
 * removes, second is the number of constraints removed due to the 0 factor
 * removed.
 */
template <typename constraintsT>
std::tuple<long, long>
remove_element_with_factor_0(constraintsT& csts) noexcept
{
    std::size_t element_removed{ 0 }, constraint_removed{ 0 };
    std::size_t size;

    for (auto& elem : csts) {
        size = elem.elements.size();
        elem.elements.erase(
          std::remove_if(elem.elements.begin(),
                         elem.elements.end(),
                         [](const auto& e) { return e.factor == 0; }),
          elem.elements.end());

        element_removed += (size - elem.elements.size());
    }

    size = csts.size();
    csts.erase(
      std::remove_if(csts.begin(),
                     csts.end(),
                     [](const auto& e) { return e.elements.empty(); }),
      csts.end());

    constraint_removed = size - csts.size();

    return std::make_tuple(numeric_cast<long>(element_removed),
                           numeric_cast<long>(constraint_removed));
}

std::vector<merged_constraint>
make_merged_constraints(std::shared_ptr<context> ctx,
                        const lp::problem& pb,
                        const parameters& params)
{
    std::vector<merged_constraint> ret;
    std::unordered_map<std::vector<lp::function_element>,
                       std::size_t,
                       merged_constraint_hash>
      cache;

    const std::size_t origin_constraints_number{
        pb.equal_constraints.size() + pb.less_equal_constraints.size() +
        pb.greater_equal_constraints.size()
    };

    //
    // Merge less and greater equal constraints if function elements are the
    // same.
    //

    for (const auto& elem : pb.equal_constraints) {
        cache.emplace(elem.elements, ret.size());
        ret.emplace_back(elem.elements,
                         numeric_cast<int>(std::lround(elem.value)),
                         numeric_cast<int>(std::lround(elem.value)));
    }

    for (const auto& elem : pb.less_equal_constraints) {
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

    for (const auto& elem : pb.greater_equal_constraints) {
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

    ctx->info(
      "  - removed constraints (merged less and greater operator): %ld\n",
      numeric_cast<long>(origin_constraints_number - ret.size()));

    //
    // Remove element from constraint functions where the factor equal 0.
    //

    {
        auto removed = remove_element_with_factor_0(ret);

        ctx->info("  - removed elements in constraints: %ld\n"
                  "  - removed empty functions in constraints: %ld\n",
                  std::get<0>(removed),
                  std::get<1>(removed));
    }

    {
        auto removed = remove_small_constraints(ret);

        ctx->info("  - removed small constraints: %ld\n", removed);
    }

    /* Compute metrics from constraints and variables uses.
       - vars counts the number of use of the variables
       - linkvars the link between vars
       - linkcst ?? */

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

    std::vector<std::pair<merged_constraint, long>> tosort;

    if (params.preprocessing == "variables-number") {
        // Algorithm to a tosort vector according to the number variables used
        // by each constraints.

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

        return ret;
    } else if (params.preprocessing == "variables-weight") {
        std::sort(
          ret.begin(), ret.end(), [vars](const auto& lhs, const auto& rhs) {
              long sumlhs{ 0 };
              long sumrhs{ 0 };

              for (auto& f : lhs.elements)
                  sumlhs += vars[f.variable_index];

              for (auto& f : rhs.elements)
                  sumrhs += vars[f.variable_index];

              return sumlhs < sumrhs;
          });

        return ret;
    } else if (params.preprocessing == "constraints-weight") {
        std::sort(ret.begin(),
                  ret.end(),
                  [linkvars, linkcst](const auto& lhs, const auto& rhs) {
                      long sumlhs{ 1 };
                      long sumrhs{ 1 };
                      // std::size_t i, e;

                      for (auto& f : lhs.elements)
                          sumlhs *= linkcst[f.variable_index].size();

                      for (auto& f : rhs.elements)
                          sumrhs *= linkcst[f.variable_index].size();

                      return sumrhs > sumlhs;
                  });

        return ret;
    } else if (params.preprocessing == "implied") {
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
    } else {
        if (not params.preprocessing.empty()) {
            ctx->warning("Unknown preprocessing `%s'\n",
                         params.preprocessing.c_str());
        }

        return ret;
    }

    // return ret;
}

template <typename modeT, typename randomT>
struct solver
{
    using mode_type = modeT;
    using random_generator_type = randomT;

    random_generator_type& rng;
    std::vector<constraint_calculator<modeT, randomT>> row_updaters;
    index m;
    index n;
    AP_type ap;
    b_type b;
    c_type c;
    x_type x;
    pi_type pi;

    solver(random_generator_type& rng_,
           index n_,
           const c_type& c_,
           const std::vector<merged_constraint>& csts)
      : rng(rng_)
      , m(csts.size())
      , n(n_)
      , ap(m, n)
      , b(b_type::Zero(2, m))
      , c(c_)
      , x(x_type::Zero(n))
      , pi(pi_type::Zero(m))
    {
        for (std::size_t i{ 0 }, e{ csts.size() }; i != e; ++i) {
            int lower{ 0 }, upper{ 0 };

            for (const auto& cst : csts[i].elements) {
                ap.set(i, cst.variable_index, cst.factor, 0.0);

                if (cst.factor > 0)
                    ++upper;
                if (cst.factor < 0)
                    --lower;
            }

            if (csts[i].min == std::numeric_limits<int>::min())
                b(0, i) = lower;
            else {
                if (csts[i].min < 0)
                    b(0, i) = std::max(csts[i].min, lower);
                else
                    b(0, i) = csts[i].min;
            }

            if (csts[i].max == std::numeric_limits<int>::max())
                b(1, i) = upper;
            else {
                if (csts[i].max > 0)
                    b(1, i) = std::min(csts[i].max, upper);
                else
                    b(1, i) = csts[i].max;
            }
        }

        ap.sort();

        init(rng_);
    }

    void reinit(random_generator_type& rng_)
    {
        std::fill(ap.P().begin(), ap.P().end(), 0.0);

        pi = pi_type::Zero(m);

        row_updaters.clear();
        init(rng_);

        std::bernoulli_distribution d(0.5);

        for (index i = 0; i != n; ++i)
            x(i) = d(rng_);
    }

    void reinit(random_generator_type& rng_, const x_type& best_previous)
    {
        std::fill(ap.P().begin(), ap.P().end(), 0.0);

        pi = pi_type::Zero(m);

        row_updaters.clear();
        init(rng_);

        x = best_previous;
        std::bernoulli_distribution d(0.5);

        for (index i = 0; i != n; ++i)
            x(i) = d(rng_);
    }

    void init(random_generator_type& rng_)
    {
        for (index i = 0; i != n; ++i)
            x(i) = init_x(c(i), mode_type());

        for (index k = 0; k != m; ++k)
            row_updaters.emplace_back(rng_, k, m, n, ap, b, c, x, pi);
    }

    void serialize(std::ostream& os) const
    {
        std::vector<index> vec;

        for (index k{ 0 }; k != m; ++k) {
            int v = 0;

            const auto& ak{ ap.row(k) };
            const auto& values{ ap.A() };

            for (std::size_t i{ 0 }, e{ ak.size() }; i != e; ++i)
                v += values[ak[i].value] * x(ak[i].position);

            if (not(b(0, k) <= v and v <= b(1, k)))
                vec.push_back(k);
        }

        for (auto it{ vec.cbegin() }, et{ vec.cend() }; it != et; ++it)
            row_updaters[*it].serialize(*it, os);
    }

    bool is_valid_solution() const noexcept
    {
        for (index k{ 0 }; k != m; ++k) {
            int v{ 0 };
            const auto& ak{ ap.row(k) };
            const auto& values{ ap.A() };

            for (std::size_t i{ 0 }, e{ ak.size() }; i != e; ++i)
                v += values[ak[i].value] * x(ak[i].position);

            if (not(b(0, k) <= v and v <= b(1, k)))
                return false;
        }

        return true;
    }

    result results(const c_type& original_costs, const int cost_constant) const
    {
        result ret;

        if (is_valid_solution()) {
            ret.status = result_status::success;
            int value = cost_constant;
            for (index i = 0; i != n; ++i)
                value += original_costs[i] * x[i];

            ret.value = static_cast<double>(value);

            ret.variable_value.resize(n, 0);

            for (index i = 0; i != n; ++i)
                ret.variable_value[i] = x(i);
        }

        ret.variables = n;
        ret.constraints = m;

        return ret;
    }
};

template <typename T>
std::size_t
cycle_avoidance_hash(const std::vector<T>& vec)
{
    std::size_t seed{ vec.size() };

    for (auto& elem : vec)
        seed ^= elem + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    return seed;
}

std::size_t
cycle_avoidance_hash(const x_type& vec)
{
    std::size_t seed{ numeric_cast<std::size_t>(vec.rows()) };

    for (std::size_t i{ 0 }, e{ seed }; i != e; ++i)
        seed ^= vec(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    return seed;
}

std::size_t
cycle_avoidance_hash(const std::vector<std::pair<index, index>>& vec)
{
    std::size_t seed{ vec.size() };

    for (auto& elem : vec)
        seed ^= elem.first + 0x9e3779b9 + (seed << 6) + (seed >> 2);

    return seed;
}

template <typename T>
struct no_cycle_avoidance
{
    using value_type = T;

    no_cycle_avoidance(std::size_t limit_ = 0) noexcept
    {
        (void)limit_;
    }

    bool have_cycle(const x_type&) const noexcept
    {
        return false;
    }

    bool have_cycle(const std::vector<T>&) const noexcept
    {
        return false;
    }

    bool have_cycle(const std::vector<std::pair<index, index>>&) const noexcept
    {
        return false;
    }
};

template <typename T>
struct cycle_avoidance
{
    using value_type = T;

    std::vector<std::size_t> history;
    std::size_t limit;
    index nb;

    cycle_avoidance(std::size_t limit_ = 48l)
      : history(limit_)
      , limit(limit_)
      , nb(0)
    {
        history.clear();
    }

    ~cycle_avoidance()
    {
    }

    bool have_cycle()
    {
        long distance{ 2 };
        const long end{ numeric_cast<long>(history.size()) };

        if (end < distance)
            return false;

        for (; distance != end; ++distance) {
            index i{ end - 1 };
            index j{ i - distance };
            long cycle_size{ 0 };

            while (j >= 0 and history[i] == history[j]) {
                cycle_size++;
                --j;
                --i;
            }

            if (cycle_size > 1) {
                history.clear();
                nb++;
                return true;
            }
        }

        if (history.size() > limit)
            history.erase(history.begin(),
                          history.begin() + (history.size() - limit));

        return false;
    }

    bool have_cycle(const x_type& x)
    {
        history.push_back(cycle_avoidance_hash(x));
        return have_cycle();
    }

    bool have_cycle(const std::vector<T>& R)
    {
        history.push_back(cycle_avoidance_hash(R));
        return have_cycle();
    }

    bool have_cycle(const std::vector<std::pair<index, index>>& R)
    {
        history.push_back(cycle_avoidance_hash(R));
        return have_cycle();
    }
};

template <typename solverT>
std::size_t
compute_missing_constraint(solverT& solver, std::vector<index>& R)
{
    R.clear();

    for (index k{ 0 }; k != solver.m; ++k) {
        int v = 0;
        const auto& ak{ solver.ap.row(k) };
        const auto& values{ solver.ap.A() };

        for (std::size_t i{ 0 }, e{ ak.size() }; i != e; ++i)
            v += values[ak[i].value] * solver.x(ak[i].position);

        if (not(solver.b(0, k) <= v and v <= solver.b(1, k)))
            R.push_back(k);
    }

    return R.size();
}

template <typename randomT>
struct compute_reversing
{
    std::vector<index> R;
    no_cycle_avoidance<index> detect_infeasability_cycle;

    compute_reversing(randomT&)
    {
    }

    template <typename solverT>
    std::size_t run_all(solverT& solver,
                        double kappa,
                        double delta,
                        double theta)
    {
        for (index i{ solver.m - 1 }; i >= 0; --i)
            solver.row_updaters[i].update_row(i, kappa, delta, theta);

        return compute_missing_constraint(solver, R);
    }

    template <typename solverT>
    std::size_t run(solverT& solver, double kappa, double delta, double theta)
    {
        auto ret = compute_missing_constraint(solver, R);
        if (ret == 0)
            return 0;

        if (detect_infeasability_cycle.have_cycle(R))
            return run_all(solver, 0.8 * kappa, delta, theta);

        for (auto it{ R.crbegin() }, et{ R.crend() }; it != et; ++it)
            solver.row_updaters[*it].update_row(*it, kappa, delta, theta);

        return compute_missing_constraint(solver, R);
    }
};

template <typename randomT>
struct compute_none
{
    std::vector<index> R;
    no_cycle_avoidance<index> detect_infeasability_cycle;

    compute_none(randomT&)
    {
    }

    template <typename solverT>
    std::size_t run_all(solverT& solver,
                        double kappa,
                        double delta,
                        double theta)
    {
        for (index i{ 0 }; i != solver.m; ++i)
            solver.row_updaters[i].update_row(i, kappa, delta, theta);

        return compute_missing_constraint(solver, R);
    }

    template <typename solverT>
    std::size_t run(solverT& solver, double kappa, double delta, double theta)
    {
        auto ret = compute_missing_constraint(solver, R);
        if (ret == 0)
            return 0;

        if (detect_infeasability_cycle.have_cycle(R))
            return run_all(solver, 0.8 * kappa, delta, theta);

        for (auto it{ R.cbegin() }, et{ R.cend() }; it != et; ++it)
            solver.row_updaters[*it].update_row(*it, kappa, delta, theta);

        return compute_missing_constraint(solver, R);
    }
};

template <typename randomT>
struct compute_random
{
    using random_generator_type = randomT;

    no_cycle_avoidance<index> detect_infeasability_cycle;
    random_generator_type& rng;
    std::vector<index> R;

    compute_random(random_generator_type& rng_)
      : rng(rng_)
    {
    }

    template <typename solverT>
    std::size_t run_all(solverT& solver,
                        double kappa,
                        double delta,
                        double theta)
    {
        R.resize(solver.m);
        std::iota(R.begin(), R.end(), 0);
        std::shuffle(R.begin(), R.end(), rng);

        for (auto it{ R.cbegin() }, et{ R.cend() }; it != et; ++it)
            solver.row_updaters[*it].update_row(*it, kappa, delta, theta);

        return compute_missing_constraint(solver, R);
    }

    template <typename solverT>
    std::size_t run(solverT& solver, double kappa, double delta, double theta)
    {
        auto ret = compute_missing_constraint(solver, R);
        if (ret == 0)
            return 0;

        if (detect_infeasability_cycle.have_cycle(R))
            return run_all(solver, 0.8 * kappa, delta, theta);

        std::shuffle(R.begin(), R.end(), rng);
        for (auto it{ R.cbegin() }, et{ R.cend() }; it != et; ++it)
            solver.row_updaters[*it].update_row(*it, kappa, delta, theta);

        return compute_missing_constraint(solver, R);
    }
};

struct compute_infeasibility_incr
{
};

struct compute_infeasibility_decr
{
};

template <typename randomT, typename directionT>
struct compute_infeasibility
{
    using random_generator_type = randomT;
    using direction_type = directionT;

    random_generator_type& rng;
    std::vector<std::pair<index, index>> R;
    no_cycle_avoidance<index> detect_infeasability_cycle;

    template <typename iteratorT>
    void sort(iteratorT begin, iteratorT end, compute_infeasibility_incr) const
    {
        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return lhs.second < rhs.second;
        });
    }

    template <typename iteratorT>
    void sort(iteratorT begin, iteratorT end, compute_infeasibility_decr) const
    {
        std::sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return rhs.second < lhs.second;
        });
    }

    compute_infeasibility(random_generator_type& rng_)
      : rng(rng_)
    {
    }

    template <typename solverT>
    std::size_t compute_missing_constraint(solverT& solver)
    {
        R.clear();

        for (index k{ 0 }; k != solver.m; ++k) {
            int v = 0;

            const auto& ak{ solver.ap.row(k) };
            const auto& values{ solver.ap.A() };

            for (std::size_t i{ 0 }, e{ ak.size() }; i != e; ++i)
                v += values[ak[i].value] * solver.x(ak[i].position);

            if (solver.b(0, k) > v)
                R.emplace_back(k, solver.b(0, k) - v);

            if (solver.b(1, k) < v)
                R.emplace_back(k, v - solver.b(1, k));
        }

        return R.size();
    }

    template <typename solverT>
    std::size_t run_all(solverT& solver,
                        double kappa,
                        double delta,
                        double theta)
    {
        R.clear();

        for (index k{ 0 }; k != solver.m; ++k) {
            int v = 0;

            const auto& ak{ solver.ap.row(k) };
            const auto& values{ solver.ap.A() };

            for (std::size_t i{ 0 }, e{ ak.size() }; i != e; ++i)
                v += values[ak[i].value] * solver.x(ak[i].position);

            if (solver.b(0, k) > v)
                R.emplace_back(k, solver.b(0, k) - v);
            else if (solver.b(1, k) < v)
                R.emplace_back(k, v - solver.b(1, k));
            else
                R.emplace_back(k, 0);
        }

        sort(R.begin(), R.end(), direction_type());

        auto ret = R.begin();
        auto it = ret + 1;
        for (; it != R.end(); ++it) {
            if (ret->second != it->second) {
                std::shuffle(ret, it, rng);
                ret = it;
            }
        }

        std::shuffle(ret, R.end(), rng);
        for (std::size_t i{ 0 }, e = { R.size() }; i != e; ++i)
            solver.row_updaters[R[i].first].update_row(
              R[i].first, kappa, delta, theta);

        return compute_missing_constraint(solver);
    }

    template <typename solverT>
    std::size_t run(solverT& solver, double kappa, double delta, double theta)
    {
        auto ret = compute_missing_constraint(solver);
        if (ret == 0)
            return 0;

        if (detect_infeasability_cycle.have_cycle(R))
            return run_all(solver, 0.8 * kappa, delta, theta);

        sort(R.begin(), R.end(), direction_type());

        auto first = R.begin();
        auto it = first + 1;
        for (; it != R.end(); ++it) {
            if (first->second != it->second) {
                std::shuffle(first, it, rng);
                first = it;
            }
        }

        std::shuffle(first, R.end(), rng);
        for (std::size_t i{ 0 }, e = { R.size() }; i != e; ++i)
            solver.row_updaters[R[i].first].update_row(
              R[i].first, kappa, delta, theta);

        return compute_missing_constraint(solver);
    }
};

template <typename modeT, typename constraintOrderT, typename randomT>
struct solver_functor
{
    using mode_type = modeT;
    using constraint_order_type = constraintOrderT;
    using random_generator_type = randomT;

    std::chrono::time_point<std::chrono::steady_clock> m_begin;
    std::chrono::time_point<std::chrono::steady_clock> m_end;

    std::shared_ptr<context> m_ctx;
    x_type m_best_x;
    result m_best;

    solver_functor(std::shared_ptr<context> ctx)
      : m_ctx(ctx)
    {
    }

    result operator()(const std::vector<merged_constraint>& constraints,
                      index variables,
                      const c_type& original_costs,
                      const c_type& norm_costs,
                      int cost_constant,
                      const parameters& p,
                      randomT& rng)
    {
        m_begin = std::chrono::steady_clock::now();
        m_end = m_begin;

        long int i{ 0 };
        long int i2{ 0 };
        double kappa_old{ 0 };
        double kappa = p.kappa_min;
        index best_remaining{ -1 };

        solver<mode_type, random_generator_type> slv(
          rng, variables, norm_costs, constraints);

        constraint_order_type compute(rng);

        m_ctx->info("* solver starts:\n");

        for (;;) {
            index remaining = compute.run(slv, kappa, p.delta, p.theta);
            auto current = slv.results(original_costs, cost_constant);
            current.loop = i;
            current.remaining_constraints = remaining;
            current.duration =
              std::chrono::duration_cast<std::chrono::duration<double>>(
                m_end - m_begin)
                .count();

            if (best_remaining == -1 or remaining < best_remaining) {
                best_remaining = remaining;
                m_best = current;

                m_ctx->info("  - constraints remaining: %ld/%ld at %fs\n",
                            remaining,
                            current.constraints,
                            current.duration);
            }

            if (current.status == result_status::success) {
                m_ctx->info("  - Solution found: %f\n", current.value);
                m_best = current;
                return m_best;
            }

            if (i2 <= p.w) {
                kappa = p.kappa_min;
                i2++;
            } else {
                i2 = 0;
                kappa = kappa_old +
                        p.kappa_step *
                          std::pow(static_cast<double>(remaining) /
                                     static_cast<double>(current.constraints),
                                   p.alpha);
                kappa_old = kappa;
            }

            if (++i > p.limit) {
                m_ctx->info("  - Loop reached: %f\n", kappa);
                m_best.status = result_status::limit_reached;
                return m_best;
            }

            if (kappa > p.kappa_max) {
                m_ctx->info("  - Kappa reached: %f\n", kappa);
                m_best.status = result_status::kappa_max_reached;
                return m_best;
            }

            m_end = std::chrono::steady_clock::now();
            if (is_time_limit(p.time_limit, m_begin, m_end)) {
                m_ctx->info("  - Time limit reached\n");
                m_best.status = result_status::time_limit_reached;
                return m_best;
            }
        }
    }
};

template <typename modeT, typename constraintOrderT, typename randomT>
struct optimize_functor
{
    using mode_type = modeT;
    using constraint_order_type = constraintOrderT;
    using random_generator_type = randomT;

    std::chrono::time_point<std::chrono::steady_clock> m_begin;
    std::chrono::time_point<std::chrono::steady_clock> m_end;

    std::shared_ptr<context> m_ctx;
    x_type m_best_x;
    result m_best;

    optimize_functor(std::shared_ptr<context> ctx)
      : m_ctx(ctx)
    {
    }

    result operator()(const std::vector<merged_constraint>& constraints,
                      index variables,
                      const c_type& original_costs,
                      const c_type& norm_costs,
                      int cost_constant,
                      const parameters& p,
                      randomT& rng)
    {
        m_begin = std::chrono::steady_clock::now();
        m_end = m_begin;

        long int i{ 0 };
        long int i2{ 0 };
        double kappa_old{ 0 };
        double kappa = p.kappa_min;

        long int pushed{ -1 };
        long int pushing_iteration{ 0 };

        solver<mode_type, random_generator_type> slv(
          rng, variables, norm_costs, constraints);

        constraint_order_type compute(rng);

        for (; not is_time_limit(p.time_limit, m_begin, m_end);
             m_end = std::chrono::steady_clock::now(), ++i) {

            index remaining = compute.run(slv, kappa, p.delta, p.theta);

            auto current = slv.results(original_costs, cost_constant);
            current.loop = i;
            current.remaining_constraints = remaining;

            if (store_if_better(current)) {
                m_best_x = slv.x;
                pushed = 0;
            }

            if (i2 <= p.w) {
                kappa = p.kappa_min;
                i2++;
            } else {
                i2 = 0;
                kappa = kappa_old +
                        p.kappa_step *
                          std::pow(static_cast<double>(remaining) /
                                     static_cast<double>(current.constraints),
                                   p.alpha);

                kappa_old = kappa;
            }

            if (i >= p.limit or kappa > p.kappa_max or
                pushed > p.pushes_limit) {
                if (m_best.status == result_status::success) {
                    slv.reinit(rng, m_best_x);
                } else {
                    slv.reinit(rng);
                }

                i = 0;
                i2 = 0;
                kappa_old = 0.0;
                kappa = p.kappa_min;
                pushed = -1;
                pushing_iteration = 0;

                continue;
            }

            if (pushed >= 0) {
                ++pushing_iteration;

                if (pushing_iteration >= p.pushing_iteration_limit) {
                    pushed++;
                    pushing_iteration = 0;

                    auto c_copy = slv.c;
                    for (index i{ 0 }; i != slv.n; ++i)
                        slv.c(i) += slv.c(i) * p.pushing_objective_amplifier;

                    remaining = compute.run_all(
                      slv, p.pushing_k_factor * kappa, p.delta, p.theta);

                    slv.c = c_copy;

                    if (remaining == 0) {
                        current = slv.results(original_costs, cost_constant);
                        current.loop = i;
                        if (store_if_better(current))
                            m_best_x = slv.x;
                    }
                }
            }
        }

        return m_best;
    }

private:
    bool store_if_better(const result& current) noexcept
    {
        if (current.status != result_status::success)
            return false;

        if (m_best.status != result_status::success or
            is_better_solution(current.value, m_best.value, mode_type())) {

            double t =
              std::chrono::duration_cast<std::chrono::duration<double>>(
                m_end - m_begin)
                .count();

            m_ctx->info("  - Solution found: %f (i=%ld t=%fs)\n",
                        current.value,
                        current.loop,
                        t);

            m_best = current;
            m_best.duration = t;

            return true;
        }

        return false;
    }
};

/*
 * Normalizes the cost vector, i.e. divides it by its l{1,2, +oo}norm. If the
 * input vector is too small the c is unchanged.
 */
inline c_type
normalize_costs(std::shared_ptr<context> ctx, const c_type& c)
{
    c_type ret(c);

    {
        ctx->info("  -C ompute infinity-norm\n");

        double sum{ static_cast<double>(c.maxCoeff()) };

        if (std::isnormal(sum))
            return ret /= sum;
    }

    // {
    //     ctx->info("Compute l1 norm\n");
    //     int sum{ 0 };
    //     for (long i{ 0 }, e{ c.rows() }; i != e; ++i)
    //         sum += std::abs(c(i));

    //     if (std::isnormal(sum))
    //         return ret /= sum;
    // }

    // {
    //     ctx->info("Compute l2 norm\n");
    //     int sum{ 0 };
    //     for (long i{ 0 }, e{ c.rows() }; i != e; ++i)
    //         sum += c(i) * c(i);

    //     if (std::isnormal(sum))
    //         return ret /= sum;
    // }

    return ret;
}

template <typename modeT, typename constraintOrderT, typename randomT>
inline result
solve(std::shared_ptr<context> ctx,
      problem& pb,
      const parameters& p,
      randomT& rng)
{
    ctx->info("Solver initializing\n");

    auto names = pb.vars.names;
    auto constraints{ make_merged_constraints(ctx, pb, p) };
    auto variables = lp::numeric_cast<index>(pb.vars.values.size());
    auto cost = make_objective_function(pb.objective, variables);
    auto norm_costs = normalize_costs(ctx, cost);
    auto cost_constant = pb.objective.constant;
    pb.clear();

    solver_functor<modeT, constraintOrderT, randomT> slv(ctx);

    auto result =
      slv(constraints, variables, cost, norm_costs, cost_constant, p, rng);

    result.method = "inequalities_1coeff solver";
    result.variable_name = names;

    return result;
}

template <typename modeT, typename constraintOrderT, typename randomT>
inline result
optimize(std::shared_ptr<context> ctx,
         problem& pb,
         const parameters& p,
         randomT& rng,
         long int thread)
{
    Expects(thread >= 1, "optimize: bad thread number");

    ctx->info("Optimizer initializing\n");

    auto names = pb.vars.names;
    auto constraints{ make_merged_constraints(ctx, pb, p) };
    auto variables = lp::numeric_cast<index>(pb.vars.values.size());
    auto cost = make_objective_function(pb.objective, variables);
    auto norm_costs = normalize_costs(ctx, cost);
    auto cost_constant = pb.objective.constant;
    pb.clear();

    std::vector<std::thread> pool(thread);
    pool.clear();
    std::vector<std::future<result>> results(thread);
    results.clear();

    if (thread == 1)
        ctx->info("optimizer starts with one thread\n");
    else
        ctx->info("Optimizer starts with %ld threads\n", thread);

    for (long int i{ 0 }; i != thread; ++i) {
        std::packaged_task<result()> task(
          std::bind(optimize_functor<modeT, constraintOrderT, randomT>(ctx),
                    std::ref(constraints),
                    variables,
                    std::ref(cost),
                    std::ref(norm_costs),
                    cost_constant,
                    std::ref(p),
                    std::ref(rng)));

        results.emplace_back(task.get_future());

        pool.emplace_back(std::thread(std::move(task)));
    }

    for (auto& t : pool)
        t.join();

    result best = results[0].get();
    for (long int i{ 1 }; i != thread; ++i) {
        auto current = results[i].get();
        if (current.status == lp::result_status::success) {
            if (best.status != lp::result_status::success or
                is_better_solution(current.value, best.value, modeT()))
                best = current;
        }
    }

    best.method = "inequalities_1coeff optimizer";
    best.variable_name = names;

    return best;
}

} // inequalities_1coeff

inline result
inequalities_1coeff_wedelin_solve(
  std::shared_ptr<lp::context> ctx,
  problem& pb,
  const std::map<std::string, parameter>& params)
{
    namespace ine_1 = lp::inequalities_1coeff;
    ine_1::parameters p(ctx, params);

    using random_generator_type = std::default_random_engine;

    //
    // TODO we need to add parameters to select the type of the generator to
    // use several type of PRNG.
    //

    random_generator_type::result_type seed = ine_1::get_integer(
      ctx,
      params,
      "seed",
      std::chrono::system_clock::now().time_since_epoch().count());

    random_generator_type rng(seed);

    switch (p.order) {
    case ine_1::constraint_order::none:
        if (pb.type == lp::objective_function_type::maximize)
            return ine_1::solve<ine_1::maximize_tag,
                                ine_1::compute_none<random_generator_type>,
                                random_generator_type>(ctx, pb, p, rng);
        return ine_1::solve<ine_1::minimize_tag,
                            ine_1::compute_none<random_generator_type>,
                            random_generator_type>(ctx, pb, p, rng);
    case ine_1::constraint_order::reversing:
        if (pb.type == lp::objective_function_type::maximize)
            return ine_1::solve<
              ine_1::maximize_tag,
              ine_1::compute_reversing<random_generator_type>>(
              ctx, pb, p, rng);
        return ine_1::solve<ine_1::minimize_tag,
                            ine_1::compute_reversing<random_generator_type>>(
          ctx, pb, p, rng);
    case ine_1::constraint_order::random_sorting:
        if (pb.type == lp::objective_function_type::maximize)
            return ine_1::solve<ine_1::maximize_tag,
                                ine_1::compute_random<random_generator_type>>(
              ctx, pb, p, rng);
        return ine_1::solve<ine_1::minimize_tag,
                            ine_1::compute_random<random_generator_type>>(
          ctx, pb, p, rng);
    case ine_1::constraint_order::infeasibility_decr:
        if (pb.type == lp::objective_function_type::maximize)
            return ine_1::solve<
              ine_1::maximize_tag,
              ine_1::compute_infeasibility<random_generator_type,
                                           ine_1::compute_infeasibility_decr>>(
              ctx, pb, p, rng);
        return ine_1::solve<
          ine_1::minimize_tag,
          ine_1::compute_infeasibility<random_generator_type,
                                       ine_1::compute_infeasibility_decr>>(
          ctx, pb, p, rng);
    case ine_1::constraint_order::infeasibility_incr:
        if (pb.type == lp::objective_function_type::maximize)
            return ine_1::solve<
              ine_1::maximize_tag,
              ine_1::compute_infeasibility<random_generator_type,
                                           ine_1::compute_infeasibility_incr>>(
              ctx, pb, p, rng);
        return ine_1::solve<
          ine_1::minimize_tag,
          ine_1::compute_infeasibility<random_generator_type,
                                       ine_1::compute_infeasibility_incr>>(
          ctx, pb, p, rng);
    }

    ctx->error("inequalities_1coeff_wedelin_solve: internal error");

    return {};
}

inline result
inequalities_1coeff_wedelin_optimize(
  std::shared_ptr<lp::context> ctx,
  problem& pb,
  const std::map<std::string, parameter>& params,
  long int thread)
{

    namespace ine_1 = lp::inequalities_1coeff;
    ine_1::parameters p(ctx, params);

    using random_generator_type = std::default_random_engine;

    //
    // TODO we need to add parameters to select the type of the generator to
    // use several type of PRNG.
    //

    random_generator_type::result_type seed = ine_1::get_integer(
      ctx,
      params,
      "seed",
      std::chrono::system_clock::now().time_since_epoch().count());

    random_generator_type rng(seed);

    switch (p.order) {
    case ine_1::constraint_order::none:
        if (pb.type == lp::objective_function_type::maximize)
            return ine_1::optimize<ine_1::maximize_tag,
                                   ine_1::compute_none<random_generator_type>,
                                   random_generator_type>(
              ctx, pb, p, rng, thread);
        return ine_1::optimize<ine_1::minimize_tag,
                               ine_1::compute_none<random_generator_type>,
                               random_generator_type>(ctx, pb, p, rng, thread);
    case ine_1::constraint_order::reversing:
        if (pb.type == lp::objective_function_type::maximize)
            return ine_1::optimize<
              ine_1::maximize_tag,
              ine_1::compute_reversing<random_generator_type>>(
              ctx, pb, p, rng, thread);
        return ine_1::optimize<
          ine_1::minimize_tag,
          ine_1::compute_reversing<random_generator_type>>(
          ctx, pb, p, rng, thread);
    case ine_1::constraint_order::random_sorting:
        if (pb.type == lp::objective_function_type::maximize)
            return ine_1::optimize<
              ine_1::maximize_tag,
              ine_1::compute_random<random_generator_type>>(
              ctx, pb, p, rng, thread);
        return ine_1::optimize<ine_1::minimize_tag,
                               ine_1::compute_random<random_generator_type>>(
          ctx, pb, p, rng, thread);
    case ine_1::constraint_order::infeasibility_decr:
        if (pb.type == lp::objective_function_type::maximize)
            return ine_1::optimize<
              ine_1::maximize_tag,
              ine_1::compute_infeasibility<random_generator_type,
                                           ine_1::compute_infeasibility_decr>>(
              ctx, pb, p, rng, thread);
        return ine_1::optimize<
          ine_1::minimize_tag,
          ine_1::compute_infeasibility<random_generator_type,
                                       ine_1::compute_infeasibility_decr>>(
          ctx, pb, p, rng, thread);
    case ine_1::constraint_order::infeasibility_incr:
        if (pb.type == lp::objective_function_type::maximize)
            return ine_1::optimize<
              ine_1::maximize_tag,
              ine_1::compute_infeasibility<random_generator_type,
                                           ine_1::compute_infeasibility_incr>>(
              ctx, pb, p, rng, thread);
        return ine_1::optimize<
          ine_1::minimize_tag,
          ine_1::compute_infeasibility<random_generator_type,
                                       ine_1::compute_infeasibility_incr>>(
          ctx, pb, p, rng, thread);
    }

    ctx->error("inequalities_1coeff_wedelin_solve: internal error");

    return {};
}

} // lp

#endif
