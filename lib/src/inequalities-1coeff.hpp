/* Copyright (C) 2017 INRA
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
#include <iostream>
#include <iterator>
#include <random>
#include <set>
#include <unordered_map>

namespace lp {
namespace inequalities_1coeff {

using A_type = lp::SparseArray<int>;
using b_type = Eigen::Matrix<double, 2, Eigen::Dynamic>;
using c_type = Eigen::VectorXf;
using x_type = Eigen::VectorXi;
using P_type = lp::SparseArray<double>;
using pi_type = Eigen::VectorXf;
using u_type = Eigen::VectorXi;

enum class constraint_order
{
    none,
    reversing,
    random_sorting,
    infeasibility_decr,
    infeasibility_incr,
    adaptative
};

const char*
constraint_order_to_string(constraint_order type)
{
    static const char* ret[] = { "none",
                                 "reversing",
                                 "random-sorting",
                                 "infeasibility-decr",
                                 "infeasibility-incr",
                                 "adaptative" };

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
        case inequalities_1coeff::constraint_order::adaptative:
            return ret[5];
    }

    return nullptr;
}

double
get_real(const std::map<std::string, parameter>& params,
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

    std::fprintf(
      stderr, " [FAIL] fail to convert parameter %s\n", param.c_str());

    return def;
}

long int
get_integer(const std::map<std::string, parameter>& params,
            std::string param,
            double def)
{
    auto it = params.find(param);
    if (it == params.cend())
        return def;

    if (it->second.type == parameter::tag::integer)
        return it->second.l;

    if (it->second.type == parameter::tag::real)
        return static_cast<long>(it->second.d);

    std::fprintf(
      stderr, " [FAIL] fail to convert parameter %s\n", param.c_str());

    return def;
}

inequalities_1coeff::constraint_order
get_constraint_order(const std::map<std::string, parameter>& params,
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
    if (it->second.s == "adaptative")
        return inequalities_1coeff::constraint_order::adaptative;

    return def;
}

struct parameters
{
    parameters(const std::map<std::string, parameter>& params)
      : time_limit(get_real(params, "time-limit", -1.0))
      , theta(get_real(params, "theta", 0.5))
      , delta(get_real(params, "delta", 0.01))
      , kappa_min(get_real(params, "kappa-min", 0.0))
      , kappa_step(get_real(params, "kappa-step", 1.e-3))
      , kappa_max(get_real(params, "kappa-max", 0.6))
      , alpha(get_real(params, "alpha", 1.0))
      , pushing_k_factor(get_real(params, "pushing-k-factor", 0.9))
      , pushes_limit(get_integer(params, "pushes-limit", 100l))
      , pushing_objective_amplifier(
          get_integer(params, "pushing-objective-amplifier", 5l))
      , pushing_iteration_limit(
          get_integer(params, "pushing-iteration-limit", 50l))
      , limit(get_integer(params, "limit", 1000l))
      , w(get_integer(params, "w", 20l))
      , order(get_constraint_order(params,
                                   "constraint-order",
                                   constraint_order::none))
      , serialize(get_integer(params, "serialize", 0l))
    {
    }

    void print() const noexcept
    {
        printf("* solver inequalities_1coeff_wedelin\n"
               "  - constraint-order: %s\n"
               "  - time_limit: %.10g\n"
               "  - theta: %.10g\n"
               "  - delta: %.10g\n"
               "  - limit: %ld\n"
               "  - kappa-min: %.10g\n"
               "  - kappa-step: %.10g\n"
               "  - kappa-max: %.10g\n"
               "  - alpha: %.10g\n"
               "  - w: %ld\n"
               "  - serialise: %d\n",
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
               serialize);
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
        std::stable_sort(begin, end, [](const auto& lhs, const auto& rhs) {
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
        std::stable_sort(begin, end, [](const auto& lhs, const auto& rhs) {
            return rhs.value < lhs.value;
        });

        random_shuffle_unique(begin, end, rng);
    }
}

bool
stop_iterating(double value, minimize_tag)
{
    return value > 0;
}

bool
stop_iterating(double value, maximize_tag)
{
    return value < 0;
}

double default_solution_value(minimize_tag)
{
    return std::numeric_limits<double>::infinity();
}

double default_solution_value(maximize_tag)
{
    return -std::numeric_limits<double>::infinity();
}

bool
is_better_solution(double current, double old, minimize_tag)
{
    return current < old;
}

bool
is_better_solution(double current, double old, maximize_tag)
{
    return current > old;
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
    A_type& A;
    b_type& b;
    c_type& cost;
    x_type& x;
    P_type& P;
    pi_type& pi;
    std::vector<index> I; // Stores variables with non null coefficient.
    std::vector<r_data> r;
    std::vector<index> C; // Stores variables with negative coefficient.
    index m;
    index n;

    constraint_calculator(random_generator_type& rng_,
                          index k,
                          index m_,
                          index n_,
                          A_type& A_,
                          b_type& b_,
                          c_type& c_,
                          x_type& x_,
                          P_type& P_,
                          pi_type& pi_)
      : rng(rng_)
      , A(A_)
      , b(b_)
      , cost(c_)
      , x(x_)
      , P(P_)
      , pi(pi_)
      , m(m_)
      , n(n_)
    {
        // for (index i = 0; i != n; ++i) {
        //     if (A(k, i)) {
        //         I.emplace_back(i);
        //         r.emplace_back(0.0, i);
        //     }

        //     if (A(k, i) < 0)
        //         C.emplace_back(i);
        // }

        const auto& ak{ A.row(k) };
        const auto& values{ A.values() };

        for (std::size_t i{ 0 }, e{ ak.size() }; i != e; ++i) {
            if (values[ak[i].value]) {
                I.emplace_back(ak[i].position);
                r.emplace_back(0.0, ak[i].position);
            }

            if (values[ak[i].value] < 0)
                C.emplace_back(ak[i].position);
        }
    }

    void serialize(index k, std::ostream& os) const
    {
        os << "[k: " << k << " pi(k): " << pi(k) << " P(k, i): ";
        for (auto i : I)
            os << P(k, i) << ' ';
        os << "] ### ";

        os << b(0, k) << " <= ";

        for (index i = 0, endi = numeric_cast<index>(I.size()); i != endi; ++i)
            os << A(k, I[i]) << ' ';

        os << " <= " << b(1, k) << '\n';
    }

    void update_row(index k, double kappa, double delta, double theta)
    {
        for (std::size_t i{ 0 }, endi{ I.size() }; i != endi; ++i)
            P.set(k, I[i], P(k, I[i]) * theta);

        std::uniform_real_distribution<> dst(0.0, 1e-4);
        for (std::size_t i{ 0 }, endi{ I.size() }; i != endi; ++i) {
            double sum_a_pi{ 0 };
            double sum_a_p{ 0 };

            const auto& ai{ A.column(I[i]) };
            const auto& values{ A.values() };

            for (std::size_t h{ 0 }, eh{ ai.size() }; h != eh; ++h) {
                sum_a_pi += values[ai[h].value] * pi(ai[h].position);
                sum_a_p += values[ai[h].value] * P(ai[h].position, I[i]);
            }

            // for (index h = 0; h != m; ++h) {
            //     if (A(h, I[i])) {
            //         sum_a_pi += A(h, I[i]) * pi(h);
            //         sum_a_p += A(h, I[i]) * P(h, I[i]);

            //         // We use std::abs(A(h, I[i])) that is always 1. So we
            //         use
            //         // only the pi(h) and P(h, I[i]).

            //         // sum_a_pi += pi(h);
            //         // sum_a_p += P(h, I[i]);
            //     }
            // }

            r[i].value = cost(I[i]) - sum_a_pi - sum_a_p;
            if (is_essentially_equal(r[i].value, 0.0, 1e-7))
                r[i].value += (std::signbit(r[i].value) ? -1. : 1.) * dst(rng);

            r[i].id = I[i];
        }

        //
        // Negate reduced costs and coefficients of these variables.
        //

        for (std::size_t i{ 0 }, endi{ I.size() }; i != endi; ++i) {
            if (A(k, I[i]) < 0) {
                r[i].value = -r[i].value;
                A.set(k, I[i], -A(k, I[i]));
                P.set(k, I[i], -P(k, I[i]));
            }
        }

        calculator_sort(r.begin(), r.end(), rng, mode_type());

        b(0, k) += C.size();
        b(1, k) += C.size();

        //
        // The bkmin and bkmax constraint bounds are not equal and can be
        // assigned to -infinity or +infinity. We have to scan the r vector and
        // search a value j such as b(0, k) <= Sum A(k, r[j]) < b(1, k).
        //

        index i{ 0 }, selected{ -1 }, first, second;
        const index endi{ numeric_cast<index>(r.size()) };
        int sum{ 0 };

        for (; i != endi; ++i) {
            sum += A(k, r[i].id);

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
                sum += A(k, r[i].id);

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
                P.set(k, r[j].id, P(k, r[j].id) - delta);
            }
        } else if (second >= endi) {
            for (index j{ 0 }; j < endi; ++j) {
                x(r[j].id) = 1;
                P.set(k, r[j].id, P(k, r[j].id) + delta);
            }
        } else {
            pi(k) += ((r[first].value + r[second].value) / 2.0);

            double d = delta + ((kappa / (1.0 - kappa)) *
                                (r[second].value - r[first].value));

            index j{ 0 };
            for (; j <= selected; ++j) {
                x(r[j].id) = 1;
                P.set(k, r[j].id, P(k, r[j].id) + d);
            }

            for (; j != endi; ++j) {
                x(r[j].id) = 0;
                P.set(k, r[j].id, P(k, r[j].id) - d);
            }
        }

        b(0, k) -= C.size();
        b(1, k) -= C.size();

        //
        // Clean up: correct negated costs and adjust value of negated
        // variables.
        //

        for (std::size_t i{ 0 }, e{ C.size() }; i != e; ++i) {
            A.set(k, C[i], -A(k, C[i]));
            P.set(k, C[i], -P(k, C[i]));
            x[C[i]] = 1 - x[C[i]];
        }
    }
};

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
 * remove the constraint and we remove all reference to this variable in
 * all
 * constraints.
 */
template <typename constraintsT>
long
remove_small_constraints(constraintsT& csts) noexcept
{
    for (auto& elem : csts) {
        if (elem.elements.size() == 1) {

            // if (elem.max ==)

            //     if (elem.min == elem.max max) {
            //         remove_variable(ret, elem.elements.variable_index);
            //         re
            // }
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
make_merged_constraints(const lp::problem& pb)
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
    // Merge less and greater equal constraints if function elements are
    // the
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

    printf("  - removed constraints (merged less and greater operator): %ld\n",
           numeric_cast<long>(origin_constraints_number - ret.size()));

    //
    // Remove element from constraint functions where the factor equal 0.
    //

    {
        auto removed = remove_element_with_factor_0(ret);

        printf("  - removed elements in constraints: %ld\n"
               "  - removed empty functions in constraints: %ld\n",
               std::get<0>(removed),
               std::get<1>(removed));
    }

    {
        auto removed = remove_small_constraints(ret);

        printf("  - removed small constraints: %ld\n", removed);
    }

    printf("  - constraints stored in: constraints.tmp.lp\n");
    std::ofstream ofs("constraints.tmp.lp");
    for (auto& elem : ret) {
        ofs << elem.min << " <= ";

        for (auto& f : elem.elements)
            ofs << ((f.factor < 0) ? "- " : "+ ") << f.variable_index << ' ';

        ofs << " <= " << elem.max << '\n';
    }

    //
    // We sort constraints according to the maximum presence of variables in
    // the constraints.
    // - vars counts the number of use of the variables
    // - link the link between vars

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

    {
        auto it = vars.begin(), end = vars.end();

        int moy{ 0 };
        for (it = vars.begin(); it != end; ++it)
            moy += *it;
        moy = moy / vars.size();

        // auto it = std::min_element(vars.begin(), vars.end());
        // int min{ *it };

        for (it = vars.begin(); it != end; ++it)
            *it = std::max(0, *it - moy);
    }

    // std::cout << "linkvars:\n";
    // for (std::size_t i{ 0 }, e{ linkvars.size() }; i != e; ++i) {
    //     std::cout << i << '=';
    //     for (auto elem : linkvars[i])
    //         std::cout << elem << ' ';
    //     std::cout << '\n';
    // }

    // std::cout << "linkcst:\n";
    // for (std::size_t i{ 0 }, e{ linkcst.size() }; i != e; ++i) {
    //     std::cout << i << '=';
    //     for (auto elem : linkcst[i])
    //         std::cout << elem << ' ';
    //     std::cout << '\n';
    // }

    // std::cout << "variable:\n";
    // for (auto elem : vars)
    //     std::cout << elem << ' ';
    // std::cout << '\n';

    // std::vector<std::pair<merged_constraint, long>> tosort;
    // long i{ 0 };
    // for (auto& cst : ret) {
    //     tosort.emplace_back(cst, 0);
    //     for (auto& elem : cst.elements) {
    //         for (auto& s : linkcst[elem.variable_index])
    //             tosort.back().second += linkvars[s].size();
    //     }

    //     std::cout << i++ << ' ' << tosort.back().second << '\n';
    // }

    // std::sort(
    //   tosort.begin(), tosort.end(), [](const auto& lhs, const auto& rhs) {
    //       return rhs.second < lhs.second;
    //   });

    // ret.clear();
    // for (auto& elem : tosort)
    //     ret.emplace_back(elem.first);

    // std::sort(ret.begin(),
    //           ret.end(),
    //           [linkvars, linkcst](const auto& lhs, const auto& rhs) {
    //               long sumlhs{ 1 };
    //               long sumrhs{ 1 };
    //               // std::size_t i, e;

    //               for (auto& f : lhs.elements)
    //                   sumlhs *= linkcst[f.variable_index].size();

    //               for (auto& f : rhs.elements)
    //                   sumrhs *= linkcst[f.variable_index].size();

    //               return sumrhs > sumlhs;
    //           });

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
    A_type A;
    b_type b;
    c_type c;
    x_type x;
    P_type P;
    pi_type pi;
    u_type u;
    const lp::problem& pb;

    solver(random_generator_type& rng_,
           const lp::problem& pb,
           const std::vector<merged_constraint>& csts)
      : rng(rng_)
      , m(csts.size())
      , n(pb.vars.values.size())
      , A(m, n)
      , b(b_type::Zero(2, m))
      , c(c_type::Zero(n))
      , x(x_type::Zero(n))
      , P(m, n)
      , pi(pi_type::Zero(m))
      , u(u_type::Zero(n))
      , pb(pb)
    {
        for (std::size_t i{ 0 }, e = pb.vars.values.size(); i != e; ++i)
            u(i) = pb.vars.values[i].max;

        for (std::size_t i{ 0 }, e{ csts.size() }; i != e; ++i) {
            int lower{ 0 }, upper{ 0 };

            for (const auto& cst : csts[i].elements) {
                A.set(i, cst.variable_index, cst.factor);
                P.set(i, cst.variable_index, 0.0);

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

        for (const auto& elem : pb.objective.elements) {
            c(elem.variable_index) += elem.factor;
            x(elem.variable_index) = c(elem.variable_index) <= 0 ? 1 : 0;
        }

        for (index k = 0; k != m; ++k)
            row_updaters.emplace_back(rng, k, m, n, A, b, c, x, P, pi);
    }

    void serialize(std::ostream& os) const
    {
        std::vector<index> vec;

        for (index k{ 0 }; k != m; ++k) {
            int v = 0;

            const auto& ak{ A.row(k) };
            const auto& values{ A.values() };

            for (std::size_t i{ 0 }, e{ ak.size() }; i != e; ++i)
                v += values[ak[i].value] * x(ak[i].position);

            // for (index i{ 0 }; i != n; ++i)
            //     v += A(k, i) * x(i);

            if (not(b(0, k) <= v and v <= b(1, k)))
                vec.push_back(k);
        }

        for (auto it{ vec.cbegin() }, et{ vec.cend() }; it != et; ++it)
            row_updaters[*it].serialize(*it, os);
    }

    double compute_value() const noexcept
    {
        double ret = pb.objective.constant;

        for (auto& elem : pb.objective.elements)
            ret += elem.factor * x(elem.variable_index);

        return ret;
    }

    bool is_valid_solution() const noexcept
    {
        for (index k{ 0 }; k != m; ++k) {
            int v{ 0 };
            const auto& ak{ A.row(k) };
            const auto& values{ A.values() };

            for (std::size_t i{ 0 }, e{ ak.size() }; i != e; ++i)
                v += values[ak[i].value] * x(ak[i].position);

            // for (index i{ 0 }; i != n; ++i)
            //     v += A(k, i) * x(i);

            if (not(b(0, k) <= v and v <= b(1, k)))
                return false;
        }

        return true;
    }

    result results() const
    {
        result ret;
        ret.method = "inequalities_101coeff_wedelin";
        ret.variables = n;
        ret.constraints = m;
        ret.value = compute_value();
        ret.solution_found = is_valid_solution();
        ret.variable_name.resize(n);
        ret.variable_value.resize(n, 0);

        for (index i = 0; i != n; ++i) {
            ret.variable_name[i] = pb.vars.names[i];
            ret.variable_value[i] = x(i);
        }

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

    no_cycle_avoidance(std::size_t limit_ = 0) noexcept { (void)limit_; }

    bool have_cycle(const x_type&) const noexcept { return false; }
    bool have_cycle(const std::vector<T>&) const noexcept { return false; }
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

    ~cycle_avoidance() { printf("  - Cycle: %ld\n", nb); }

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
        const auto& ak{ solver.A.row(k) };
        const auto& values{ solver.A.values() };

        for (std::size_t i{ 0 }, e{ ak.size() }; i != e; ++i)
            v += values[ak[i].value] * solver.x(ak[i].position);

        // for (index i{ 0 }; i != solver.n; ++i)
        //     v += solver.A(k, i) * solver.x(i);

        if (not(solver.b(0, k) <= v and v <= solver.b(1, k)))
            R.push_back(k);
    }

    return R.size();
}

template <typename randomT>
struct compute_reversing
{
    std::vector<index> R;
    cycle_avoidance<index> detect_infeasability_cycle;

    compute_reversing(randomT&) {}

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
    cycle_avoidance<index> detect_infeasability_cycle;

    compute_none(randomT&) {}

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

    cycle_avoidance<index> detect_infeasability_cycle;
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
    cycle_avoidance<index> detect_infeasability_cycle;

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

            const auto& ak{ solver.A.row(k) };
            const auto& values{ solver.A.values() };

            for (std::size_t i{ 0 }, e{ ak.size() }; i != e; ++i)
                v += values[ak[i].value] * solver.x(ak[i].position);

            // for (index i{ 0 }; i != solver.n; ++i)
            //     v += solver.A(k, i) * solver.x(i);

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

            const auto& ak{ solver.A.row(k) };
            const auto& values{ solver.A.values() };

            for (std::size_t i{ 0 }, e{ ak.size() }; i != e; ++i)
                v += values[ak[i].value] * solver.x(ak[i].position);

            // for (index i{ 0 }; i != solver.n; ++i)
            //     v += solver.A(k, i) * solver.x(i);

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

template <typename randomT>
struct compute_adaptative
{
    using random_generator_type = randomT;

    random_generator_type& rng;

    std::chrono::time_point<std::chrono::steady_clock> begin;
    constraint_order constraint;
    std::size_t best;
    bool have_switch;

    compute_none<random_generator_type> cpt_none;
    compute_random<random_generator_type> cpt_random;
    compute_reversing<random_generator_type> cpt_reversing;
    compute_infeasibility<random_generator_type, compute_infeasibility_incr>
      cpt_inf_incr;
    compute_infeasibility<random_generator_type, compute_infeasibility_decr>
      cpt_inf_decr;

    compute_adaptative(random_generator_type& rng_)
      : rng(rng_)
      , begin(std::chrono::steady_clock::now())
      , constraint(constraint_order::reversing)
      , best(std::numeric_limits<std::size_t>::max())
      , have_switch(false)
      , cpt_none(rng)
      , cpt_random(rng)
      , cpt_reversing(rng)
      , cpt_inf_incr(rng)
      , cpt_inf_decr(rng)
    {
    }

    template <typename solverT>
    std::size_t run_all(solverT& solver,
                        double kappa,
                        double delta,
                        double theta)
    {
        switch (constraint) {
            case constraint_order::none:
                return cpt_none.run_all(solver, kappa, delta, theta);
                break;
            case constraint_order::reversing:
                return cpt_reversing.run_all(solver, kappa, delta, theta);
                break;
            case constraint_order::random_sorting:
                return cpt_random.run_all(solver, kappa, delta, theta);
                break;
            case constraint_order::infeasibility_decr:
                return cpt_inf_incr.run_all(solver, kappa, delta, theta);
                break;
            case constraint_order::infeasibility_incr:
                return cpt_inf_decr.run_all(solver, kappa, delta, theta);
                break;
            default:
                return cpt_none.run_all(solver, kappa, delta, theta);
                break;
        }
    }

    template <typename solverT>
    std::size_t run(solverT& solver, double kappa, double delta, double theta)
    {
        std::size_t ret;

        if (not have_switch) {
            switch (constraint) {
                case constraint_order::none:
                    ret = cpt_none.run(solver, kappa, delta, theta);
                    break;
                case constraint_order::reversing:
                    ret = cpt_reversing.run(solver, kappa, delta, theta);
                    break;
                case constraint_order::random_sorting:
                    ret = cpt_random.run(solver, kappa, delta, theta);
                    break;
                case constraint_order::infeasibility_decr:
                    ret = cpt_inf_incr.run(solver, kappa, delta, theta);
                    break;
                case constraint_order::infeasibility_incr:
                    ret = cpt_inf_decr.run(solver, kappa, delta, theta);
                    break;
                default:
                    ret = cpt_none.run(solver, kappa, delta, theta);
                    break;
            }
        } else {
            have_switch = false;
            switch (constraint) {
                case constraint_order::none:
                    ret = cpt_none.run_all(solver, kappa, delta, theta);
                    break;
                case constraint_order::reversing:
                    ret = cpt_reversing.run_all(solver, kappa, delta, theta);
                    break;
                case constraint_order::random_sorting:
                    ret = cpt_random.run_all(solver, kappa, delta, theta);
                    break;
                case constraint_order::infeasibility_decr:
                    ret = cpt_inf_incr.run_all(solver, kappa, delta, theta);
                    break;
                case constraint_order::infeasibility_incr:
                    ret = cpt_inf_decr.run_all(solver, kappa, delta, theta);
                    break;
                default:
                    ret = cpt_none.run_all(solver, kappa, delta, theta);
                    break;
            }
        }

        if (ret < best) {
            begin = std::chrono::steady_clock::now();
            best = ret;
        } else {
            if (std::chrono::duration_cast<std::chrono::seconds>(
                  std::chrono::steady_clock::now() - begin)
                  .count() > 30) {
                switch (constraint) {
                    case constraint_order::none:
                        constraint = constraint_order::reversing;
                        break;
                    case constraint_order::reversing:
                        constraint = constraint_order::random_sorting;
                        break;
                    case constraint_order::random_sorting:
                        constraint = constraint_order::infeasibility_decr;
                        break;
                    case constraint_order::infeasibility_decr:
                        constraint = constraint_order::reversing;
                        break;
                    case constraint_order::infeasibility_incr:
                        constraint = constraint_order::reversing;
                        break;
                    default:
                        constraint = constraint_order::random_sorting;
                        break;
                };
                printf("     * Switch to %s\n",
                       constraint_order_to_string(constraint));
                begin = std::chrono::steady_clock::now();
                best = std::numeric_limits<std::size_t>::max();
                have_switch = true;
            }
        }

        return ret;
    }
};

template <typename modeT, typename constraintOrderT, typename randomT>
inline result
solve(const problem& pb, const parameters& p, randomT& rng)
{
    using mode_type = modeT;
    using constraint_order_type = constraintOrderT;
    using random_generator_type = randomT;

    auto begin = std::chrono::steady_clock::now();
    long int i2{ 0 };
    double kappa_old{ 0 };
    double kappa = p.kappa_min;

    solver<mode_type, random_generator_type> slv(
      rng, pb, make_merged_constraints(pb));
    result best;
    best.remaining_constraints = std::numeric_limits<index>::max();
    constraint_order_type compute(rng);

    for (long int i{ 0 }; i < p.limit; ++i) {

        index remaining = compute.run(slv, kappa, p.delta, p.theta);

        auto current = slv.results();
        current.loop = i;
        current.remaining_constraints = remaining;
        current.begin = begin;
        current.end = std::chrono::steady_clock::now();

        if (remaining < best.remaining_constraints) {
            auto t = std::chrono::duration_cast<std::chrono::duration<double>>(
              current.end - current.begin);

            printf("  - constraints remaining: %ld/%ld at %fs\n",
                   remaining,
                   current.constraints,
                   t.count());

            best = current;

            if (p.serialize) {
                std::ofstream ofs("current-solver.lp.dat");
                slv.serialize(ofs);
            }
        }

        if (current.solution_found) {
            std::cout << '\n' << current << '\n';
            return current;
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

        if (kappa > p.kappa_max) {
            std::cout << "\nFail: kappa-max reached\n";
            return best;
        }

        if (p.time_limit > 0 and
            std::chrono::duration_cast<std::chrono::duration<double>>(
              current.end - current.begin)
                .count() > p.time_limit) {
            std::cout << "\nFail: time limit\n";
            return best;
        }
    }

    std::cout << "\nFail: limit reached\n";
    return best;
}

template <typename modeT, typename constraintOrderT, typename randomT>
inline result
optimize(const problem& pb, const parameters& p, randomT& rng)
{
    using mode_type = modeT;
    using constraint_order_type = constraintOrderT;
    using random_generator_type = randomT;

    auto begin = std::chrono::steady_clock::now();
    long int i2{ 0 };
    double kappa_old{ 0 };
    double kappa = p.kappa_min;

    long int pushed{ -1 };
    long int pushing_iteration{ 0 };

    solver<mode_type, random_generator_type> slv(
      rng, pb, make_merged_constraints(pb));
    result best;
    best.remaining_constraints = std::numeric_limits<index>::max();
    best.value = default_solution_value(mode_type());
    constraint_order_type compute(rng);

    for (long int i{ 0 }; i < p.limit; ++i) {
        index remaining = compute.run(slv, kappa, p.delta, p.theta);

        auto current = slv.results();
        current.loop = i;
        current.remaining_constraints = remaining;
        current.begin = begin;
        current.end = std::chrono::steady_clock::now();
        current.value = slv.compute_value();

        if (remaining < best.remaining_constraints) {
            auto t = std::chrono::duration_cast<std::chrono::duration<double>>(
              current.end - current.begin);

            best = current;

            printf("  - constraints remaining: %ld/%ld at %fs\n",
                   remaining,
                   current.constraints,
                   t.count());

            if (best.solution_found)
                std::cout << "  - Solution value: " << current.value << '\n';

            if (p.serialize) {
                std::ofstream ofs("current-solver.lp.dat");
                slv.serialize(ofs);
            }
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

        //
        // If kappa greater than kappa_max, optimization finished.
        //

        if (kappa > p.kappa_max) {
            return best;
        }

        //
        // If time greater thant time_limit, optimization is finish.
        //

        if (p.time_limit > 0 and
            std::chrono::duration_cast<std::chrono::duration<double>>(
              current.end - current.begin)
                .count() > p.time_limit) {
            return best;
        }

        if (pushed >= 0)
            ++pushing_iteration;
        else if (current.solution_found) {
            pushed = 0;
            pushing_iteration = 0;
        }

        if (current.solution_found and
            (pushed >= 0 and pushing_iteration >= p.pushing_iteration_limit)) {

            if (is_better_solution(current.value, best.value, mode_type())) {
                std::cout << "  - Solution value: " << current.value << '\n';
                best = current;
            }

            pushed++;
            pushing_iteration = 0;

            if (pushed > p.pushes_limit)
                return best;

            auto c_copy = slv.c;
            for (index i{ 0 }; i != slv.n; ++i)
                slv.c += slv.c * p.pushing_objective_amplifier;

            remaining = compute.run_all(
              slv, p.pushing_k_factor * kappa, p.delta, p.theta);

            slv.c = c_copy;
        }
    }

    return best;
}

} // inequalities_1coeff

inline result
inequalities_1coeff_wedelin_solve(
  const problem& pb,
  const std::map<std::string, parameter>& params)
{
    using random_generator_type = std::default_random_engine;

    //
    // TODO we need to add parameters to select the type of the generator to
    // use several type of PRNG and perhaps the seed too.
    //

    random_generator_type rng(
      std::chrono::system_clock::now().time_since_epoch().count());

    namespace ine_1 = lp::inequalities_1coeff;
    ine_1::parameters p(params);
    p.print();

    switch (p.order) {
        case ine_1::constraint_order::none:
            if (pb.type == lp::objective_function_type::maximize)
                return ine_1::solve<ine_1::maximize_tag,
                                    ine_1::compute_none<random_generator_type>,
                                    random_generator_type>(pb, p, rng);
            return ine_1::solve<ine_1::minimize_tag,
                                ine_1::compute_none<random_generator_type>,
                                random_generator_type>(pb, p, rng);
        case ine_1::constraint_order::reversing:
            if (pb.type == lp::objective_function_type::maximize)
                return ine_1::solve<
                  ine_1::maximize_tag,
                  ine_1::compute_reversing<random_generator_type>>(pb, p, rng);
            return ine_1::solve<
              ine_1::minimize_tag,
              ine_1::compute_reversing<random_generator_type>>(pb, p, rng);
        case ine_1::constraint_order::random_sorting:
            if (pb.type == lp::objective_function_type::maximize)
                return ine_1::solve<
                  ine_1::maximize_tag,
                  ine_1::compute_random<random_generator_type>>(pb, p, rng);
            return ine_1::solve<ine_1::minimize_tag,
                                ine_1::compute_random<random_generator_type>>(
              pb, p, rng);
        case ine_1::constraint_order::infeasibility_decr:
            if (pb.type == lp::objective_function_type::maximize)
                return ine_1::solve<ine_1::maximize_tag,
                                    ine_1::compute_infeasibility<
                                      random_generator_type,
                                      ine_1::compute_infeasibility_decr>>(
                  pb, p, rng);
            return ine_1::solve<
              ine_1::minimize_tag,
              ine_1::compute_infeasibility<random_generator_type,
                                           ine_1::compute_infeasibility_decr>>(
              pb, p, rng);
        case ine_1::constraint_order::infeasibility_incr:
            if (pb.type == lp::objective_function_type::maximize)
                return ine_1::solve<ine_1::maximize_tag,
                                    ine_1::compute_infeasibility<
                                      random_generator_type,
                                      ine_1::compute_infeasibility_incr>>(
                  pb, p, rng);
            return ine_1::solve<
              ine_1::minimize_tag,
              ine_1::compute_infeasibility<random_generator_type,
                                           ine_1::compute_infeasibility_incr>>(
              pb, p, rng);
        case ine_1::constraint_order::adaptative:
            if (pb.type == lp::objective_function_type::maximize)
                return ine_1::solve<
                  ine_1::maximize_tag,
                  ine_1::compute_adaptative<random_generator_type>>(
                  pb, p, rng);
            return ine_1::solve<
              ine_1::minimize_tag,
              ine_1::compute_adaptative<random_generator_type>>(pb, p, rng);
    }

    throw "TODO internal error";
}

inline result
inequalities_1coeff_wedelin_optimize(
  const problem& pb,
  const std::map<std::string, parameter>& params)
{
    using random_generator_type = std::default_random_engine;

    //
    // TODO we need to add parameters to select the type of the generator to
    // use several type of PRNG and perhaps the seed too.
    //

    random_generator_type rng(
      std::chrono::system_clock::now().time_since_epoch().count());

    namespace ine_1 = lp::inequalities_1coeff;
    ine_1::parameters p(params);
    p.print();

    switch (p.order) {
        case ine_1::constraint_order::none:
            if (pb.type == lp::objective_function_type::maximize)
                return ine_1::optimize<
                  ine_1::maximize_tag,
                  ine_1::compute_none<random_generator_type>,
                  random_generator_type>(pb, p, rng);
            return ine_1::optimize<ine_1::minimize_tag,
                                   ine_1::compute_none<random_generator_type>,
                                   random_generator_type>(pb, p, rng);
        case ine_1::constraint_order::reversing:
            if (pb.type == lp::objective_function_type::maximize)
                return ine_1::optimize<
                  ine_1::maximize_tag,
                  ine_1::compute_reversing<random_generator_type>>(pb, p, rng);
            return ine_1::optimize<
              ine_1::minimize_tag,
              ine_1::compute_reversing<random_generator_type>>(pb, p, rng);
        case ine_1::constraint_order::random_sorting:
            if (pb.type == lp::objective_function_type::maximize)
                return ine_1::optimize<
                  ine_1::maximize_tag,
                  ine_1::compute_random<random_generator_type>>(pb, p, rng);
            return ine_1::optimize<
              ine_1::minimize_tag,
              ine_1::compute_random<random_generator_type>>(pb, p, rng);
        case ine_1::constraint_order::infeasibility_decr:
            if (pb.type == lp::objective_function_type::maximize)
                return ine_1::optimize<ine_1::maximize_tag,
                                       ine_1::compute_infeasibility<
                                         random_generator_type,
                                         ine_1::compute_infeasibility_decr>>(
                  pb, p, rng);
            return ine_1::optimize<
              ine_1::minimize_tag,
              ine_1::compute_infeasibility<random_generator_type,
                                           ine_1::compute_infeasibility_decr>>(
              pb, p, rng);
        case ine_1::constraint_order::infeasibility_incr:
            if (pb.type == lp::objective_function_type::maximize)
                return ine_1::optimize<ine_1::maximize_tag,
                                       ine_1::compute_infeasibility<
                                         random_generator_type,
                                         ine_1::compute_infeasibility_incr>>(
                  pb, p, rng);
            return ine_1::optimize<
              ine_1::minimize_tag,
              ine_1::compute_infeasibility<random_generator_type,
                                           ine_1::compute_infeasibility_incr>>(
              pb, p, rng);
        case ine_1::constraint_order::adaptative:
            if (pb.type == lp::objective_function_type::maximize)
                return ine_1::optimize<
                  ine_1::maximize_tag,
                  ine_1::compute_adaptative<random_generator_type>>(
                  pb, p, rng);
            return ine_1::optimize<
              ine_1::minimize_tag,
              ine_1::compute_adaptative<random_generator_type>>(pb, p, rng);
    }

    throw "TODO internal error";
}

} // lp

#endif
