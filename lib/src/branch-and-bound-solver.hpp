/* Copyright (C) 2016-2019 INRA
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

#ifndef ORG_VLEPROJECT_BARYONYX_BRANCH_AND_BOUND_HPP
#define ORG_VLEPROJECT_BARYONYX_BRANCH_AND_BOUND_HPP

#include <limits>
#include <tuple>
#include <vector>

#include "itm-common.hpp"
#include "private.hpp"
#include "utils.hpp"

#include <cstdint>

namespace baryonyx {
namespace itm {

/**
 * @brief A classic @c std::vector shared between accessors.
 *
 * @details This class only used by the branch and bound algorithm avoid the
 *     use of dynamical memory for each node. Indeed, a global @c std::vector
 *     is used to store all boolean values and each node store only an
 *     accessors to this global vector. Moreover, we use a free list to reuse
 *     accessors from a removed node.
 */
class shared_subvector
{
public:
    using container_type = std::vector<std::int8_t>;
    using value_type = container_type::value_type;
    using size_type = std::uint32_t;
    using iterator = container_type::iterator;
    using const_iterator = container_type::const_iterator;

private:
    container_type m_container;
    size_type m_element_size;
    size_type m_free_list_head;

    static inline const size_type npos = std::numeric_limits<size_type>::max();

public:
    shared_subvector()
      : m_element_size(0)
      , m_free_list_head(npos)
    {}

    void reserve(std::size_t reserve_size)
    {
        m_container.reserve(reserve_size);
    }

    [[nodiscard]] bool init(size_type element_size) noexcept
    {
        bx_assert(element_size > 0);

        m_container.resize(0);
        m_element_size = element_size;
        m_free_list_head = npos;

        if (sizeof(size_type) > (element_size * sizeof(value_type)))
            return false;

        return true;
    }

    [[nodiscard]] size_type emplace() {
        size_type new_position = 0;

        if (m_free_list_head != npos) {
            new_position = m_free_list_head;

            m_free_list_head =
              *(reinterpret_cast<size_type*>(&m_container[new_position]));

        } else {
            new_position = static_cast<size_type>(m_container.size());
            m_container.insert(std::end(m_container), m_element_size, 0);
        }

        check_index(new_position);

        return new_position;
    }

      [[nodiscard]] value_type
      operator[](size_type i) const noexcept
    {
        check_access(i);

        return m_container[i];
    }

    [[nodiscard]] value_type& operator[](size_type i) noexcept
    {
        check_access(i);

        return m_container[i];
    }

    [[nodiscard]] auto element(size_type i) const noexcept
      -> std::tuple<size_type, size_type>
    {
        check_index(i);

        return std::make_tuple(i, i + m_element_size);
    }

    [[nodiscard]] const std::int8_t* element_data(size_type i) const noexcept
    {
        check_index(i);

        return &m_container[i];
    }

    [[nodiscard]] std::int8_t* element_data(size_type i) noexcept
    {
        check_index(i);

        return &m_container[i];
    }

    void remove(size_type i)
    {
        check_index(i);

        *(reinterpret_cast<size_type*>(&m_container[i])) = m_free_list_head;
        m_free_list_head = i;
    }

    size_type size() const noexcept
    {
        return static_cast<size_type>(m_container.size());
    }

    size_type element_size() const noexcept
    {
        return m_element_size;
    }

    void fill(size_type i, std::int8_t value)
    {
        std::fill_n(&m_container[i], m_element_size, value);
    }

    void copy(size_type source, size_type destination)
    {
        std::copy_n(
          &m_container[source], m_element_size, &m_container[destination]);
    }

    void invert(size_type i, int element)
    {
        m_container[i + element] = m_container[i + element] ? 0 : 1;
    }

private:
    void check_index([[maybe_unused]] size_type i) const noexcept
    {
        bx_assert(i + m_element_size <= m_container.size());
        bx_assert(i % m_element_size == 0);
    }

    void check_access([[maybe_unused]] size_type i) const noexcept
    {
        bx_assert(i < m_container.size());
    }
};

template<typename Mode, typename Float>
struct branch_and_bound_solver
{
    struct item
    {
        Float r = 0.0;
        Float sum_z = 0.0;
        int factor = 0;
        int original_factor = 0;
        int variable = 0;
        int result = 0;

        item() noexcept = default;

        constexpr item(Float r_, int factor_, int variable_) noexcept
          : r(r_)
          , factor(factor_)
          , variable(variable_)
        {}

        constexpr bool operator<(const item& other) const noexcept
        {
            if constexpr (std::is_same_v<Mode, minimize_tag>)
                return (r / static_cast<Float>(factor)) <
                       (other.r / static_cast<Float>(factor));
            else
                return (other.r / static_cast<Float>(factor)) <
                       (r / static_cast<Float>(factor));
        }
    };

    struct node
    {
        node() noexcept = default;

        constexpr node(Float _sumr, int _sumfactor, unsigned _level) noexcept
          : z(init_z())
          , sumr(_sumr)
          , sumfactor(_sumfactor)
          , level(_level)
        {}

        unsigned int variables = 0u;
        Float z = 0.0;
        Float sumr = 0.0;
        int sumfactor = 0;
        unsigned level = 0u;
    };

    shared_subvector subvector;
    std::vector<item> items;
    std::vector<node> nodes;
    Float upper_bound;
    Float lower_bound;
    unsigned int solution;
    int b_min;
    int b_max;

    branch_and_bound_solver() = default;

    void reserve(std::size_t max_variables)
    {
        items.reserve(max_variables);
        subvector.reserve(max_variables * 1024u);
    }

    static Float init_z() noexcept
    {
        if constexpr (std::is_same_v<Mode, minimize_tag>)
            return +std::numeric_limits<Float>::infinity();
        else
            return -std::numeric_limits<Float>::infinity();
    }

    static Float compute_lower_bound() noexcept
    {
        return init_z();
    }

    Float compute_upper_bound() noexcept
    {
        Float ret = 0;

        if constexpr (std::is_same_v<Mode, minimize_tag>) {
            for (size_t i = 0, e = items.size(); i != e; ++i)
                if (items[i].r < 0)
                    ret += items[i].r;
        } else {
            for (size_t i = 0, e = items.size(); i != e; ++i)
                if (items[i].r > 0)
                    ret += items[i].r;
        }

        return ret;
    }

    bool is_best_solution(Float node_z) const noexcept
    {
        if constexpr (std::is_same_v<Mode, minimize_tag>)
            return node_z < lower_bound;
        else
            return node_z > lower_bound;
    }

    void make_init_node_0()
    {
        nodes.emplace_back(0.0, 0, 0);
        nodes.back().variables = subvector.emplace();
        nodes.back().z = items[1].sum_z;

        subvector.fill(nodes.back().variables, 0);
    }

    void make_init_node_1()
    {
        nodes.emplace_back(items[0].r, items[0].factor, 0);
        nodes.back().variables = subvector.emplace();
        nodes.back().z = items[0].sum_z;

        subvector.fill(nodes.back().variables, 0);
        subvector[nodes.back().variables] = 1;
    }

    void make_next_equality_node()
    {
        node& node_0 = nodes.back();

        const unsigned level = node_0.level + 1;
        if (level >= subvector.element_size()) {
            subvector.remove(node_0.variables);
            nodes.pop_back();
            return;
        }

        node_0.level = level;
        node_0.z = node_0.sumr;

        if (level + 1 < subvector.element_size())
            node_0.z += items[level + 1].sum_z;

        const Float sumr = node_0.sumr;
        const int sumfactor = node_0.sumfactor;
        const unsigned int variables = node_0.variables;

        node& node_1 = nodes.emplace_back(
          sumr + items[level].r, sumfactor + items[level].factor, level);
        node_1.variables = subvector.emplace();
        subvector.copy(variables, node_1.variables);
        subvector[node_1.variables + level] = 1;
        if (level + 1 < subvector.element_size())
            node_1.z = node_1.sumr + items[level + 1].sum_z;

        if (node_1.sumfactor == b_min) {
            if (is_best_solution(node_1.sumr)) {
                subvector.copy(variables, solution);
                subvector.invert(solution, level);
                lower_bound = nodes.back().sumr;
            }
            subvector.remove(node_1.variables);
            nodes.pop_back();
        }
    }

    void make_next_inequality_node()
    {
        node& node_0 = nodes.back();

        const unsigned level = node_0.level + 1;
        if (level >= subvector.element_size()) {
            subvector.remove(node_0.variables);
            nodes.pop_back();
            return;
        }

        node_0.level = level;
        node_0.z = node_0.sumr;

        if (level + 1 < subvector.element_size())
            node_0.z += items[level + 1].sum_z;

        const Float sumr = node_0.sumr;
        const int sumfactor = node_0.sumfactor;
        const unsigned int variables = node_0.variables;

        node& node_1 = nodes.emplace_back(
          sumr + items[level].r, sumfactor + items[level].factor, level);
        node_1.variables = subvector.emplace();
        subvector.copy(variables, node_1.variables);
        subvector[node_1.variables + level] = 1;
        if (level + 1 < subvector.element_size())
            node_1.z = node_1.sumr + items[level + 1].sum_z;

        if (b_min <= node_1.sumfactor && node_1.sumfactor <= b_max) {
            if (is_best_solution(node_1.sumr)) {
                subvector.copy(variables, solution);
                subvector.invert(solution, level);
                lower_bound = nodes.back().sumr;
            }
        } else if (node_1.sumfactor > b_max) {
            subvector.remove(node_1.variables);
            nodes.pop_back();
        }
    }

    void solve_equality()
    {
        lower_bound = compute_lower_bound();
        upper_bound = compute_upper_bound();

        make_init_node_0();
        if (nodes.back().sumfactor == b_min) {
            upper_bound = nodes.back().sumr;
            subvector.copy(nodes.back().variables, solution);
            subvector.remove(nodes.back().variables);
            return;
        }

        make_init_node_1();
        if (nodes.back().sumfactor == b_min) {
            upper_bound = nodes.back().sumr;
            subvector.copy(nodes.back().variables, solution);
            subvector.remove(nodes.back().variables);
            nodes.pop_back(); // Remove the branch
        }

        while (!nodes.empty())
            make_next_inequality_node();
    }

    void solve_inequality()
    {
        lower_bound = compute_lower_bound();
        upper_bound = compute_upper_bound();

        make_init_node_0();
        if (b_min <= nodes.back().sumfactor &&
            nodes.back().sumfactor <= b_max) {
            upper_bound = nodes.back().sumr;
            subvector.copy(nodes.back().variables, solution);
        } else if (nodes.back().sumfactor > b_max) {
            subvector.remove(nodes.back().variables);
            nodes.pop_back(); // Remove the branch
        }

        make_init_node_1();
        if (b_min <= nodes.back().sumfactor &&
            nodes.back().sumfactor <= b_max) {
            upper_bound = nodes.back().sumr;
            subvector.copy(nodes.back().variables, solution);
        } else if (nodes.back().sumfactor > b_max) {
            subvector.remove(nodes.back().variables);
            nodes.pop_back(); // Remove the branch
        }

        while (!nodes.empty())
            make_next_inequality_node();
    }

    template<typename R>
    int solve(R& reduced_cost, int r_size, int bk_min, int bk_max)
    {
        bx_assert(r_size >= 4);
        bx_assert(bk_min <= bk_max);

        b_min = bk_min;
        b_max = bk_max;

        items.resize(r_size);
        if (!subvector.init(r_size))
            throw solver_failure(solver_error_tag::no_solver_available);

        solution = subvector.emplace();

        bx_assert(solution == 0u);

        nodes.clear();

        for (int i = 0; i != r_size; ++i) {
            items[i].r = reduced_cost[i].value;
            items[i].variable = reduced_cost[i].id;
            items[i].original_factor = reduced_cost[i].f;
            items[i].factor = std::abs(reduced_cost[i].f);
            items[i].result = 0;
        }

        std::fill_n(
          subvector.element_data(solution), subvector.element_size(), 0);

        std::sort(std::begin(items), std::end(items));

        for (std::size_t i = 0, e = items.size(); i != e; ++i) {
            if (items[i].original_factor < 0) {
                b_min += items[i].factor;
                b_max += items[i].factor;
            }
        }

        // Compute the inversion of the sum z
        items.back().sum_z = items.back().r;
        for (int i = length(items) - 2; i >= 0; --i)
            items[i].sum_z = items[i + 1].sum_z + items[i].r;

        if (b_max == b_min)
            solve_equality();
        else
            solve_inequality();

        for (int i = 0, e = length(items); i != e; ++i) {
            if (items[i].original_factor > 0)
                items[i].result = subvector[solution + i] == 0 ? 0 : 1;
            else
                items[i].result = subvector[solution + i] == 0 ? 1 : 0;
        }

        std::sort(std::begin(items),
                  std::end(items),
                  [](const auto& lhs, const auto& rhs) {
                      if (lhs.result == rhs.result) {
                          if constexpr (std::is_same_v<Mode, minimize_tag>)
                              return lhs.r < rhs.r;
                          else
                              return lhs.r > rhs.r;
                      } else
                          return lhs.result > rhs.result;
                  });

        auto middle =
          std::find_if(std::begin(items),
                       std::end(items),
                       [](const auto& item) { return item.result == 0; });

        for (std::size_t i = 0, e = items.size(); i != e; ++i) {
            reduced_cost[i].value = items[i].r;
            reduced_cost[i].id = items[i].variable;
            reduced_cost[i].f = items[i].original_factor;
        }

        if (middle == std::end(items))
            return items[0].result == 0 ? -1 : r_size;

        return static_cast<int>(std::distance(std::begin(items), middle) - 1);
    }
};

} // namespace itm
} // namespace baryonyx

#endif // ORG_VLEPROJECT_BARYONYX_BRANCH_AND_BOUND_HPP
