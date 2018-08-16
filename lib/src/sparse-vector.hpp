/* Copyright (C) 2018 INRA
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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_SPARSE_VECTOR_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_SPARSE_VECTOR_HPP

#include <baryonyx/core>

#include <algorithm>
#include <memory>
#include <string>

namespace baryonyx {

template<typename T>
class sparse_vector
{
public:
    using value_type = T;
    using reference = value_type&;
    using const_reference = const value_type&;
    using iterator = value_type*;
    using const_iterator = const value_type*;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    using size_type = unsigned int;
    using difference_type = std::ptrdiff_t;

protected:
    using accessors = std::unique_ptr<int[]>;
    using values = std::unique_ptr<T[]>;

    accessors m_access;
    values m_values;

    int m_access_length;
    int m_values_length;

public:
    template<typename Constraint>
    explicit sparse_vector(const std::vector<Constraint>& csts)
      : m_access(std::make_unique<int[]>(csts.size() + 1))
      , m_access_length(static_cast<int>(csts.size()))
      , m_values_length(0)
    {
        //
        // Fill the accessors vector with a cumulative number of element for
        // each element of the vector.
        //
        // m_access[0] = 0
        // m_access[1] = m_access[0] + number_of_negative_element_in[0]
        // m_access[2] = m_access[1] + number_of_negative_element_in[2]
        // etc.
        //

        for (int i = 0; i != m_access_length; ++i) {
            for (const auto& elem : csts[i].elements) {
                if (elem.factor < 0) {
                    ++m_access[i + 1];
                    ++m_values_length;
                }
            }
        }

        for (int i = 1, e = m_access_length + 1; i != e; ++i)
            m_access[i] += m_access[i - 1];

        bx_assert(m_values_length > 0);
        m_values = std::make_unique<T[]>(m_values_length);

        //
        // Fill the values vector with the index of the absolute position of
        // the negative coefficient in the constraints.
        //

        for (int id = 0, i = 0; i != m_access_length; ++i) {
            int id_in_r = 0;
            for (const auto& elem : csts[i].elements) {
                if (elem.factor < 0)
                    m_values[id++].id_r = id_in_r;

                ++id_in_r;
            }
        }
    }

    bool empty(int i) const noexcept
    {
        bx_expects(i >= 0 && i < m_access_length);

        return m_access[i] == m_access[i + 1];
    }

    std::tuple<iterator, iterator> range(int i) noexcept
    {
        bx_expects(i >= 0 && i < m_access_length);

        iterator begin = m_values.get() + m_access[i];
        iterator end = m_values.get() + m_access[i + 1];

        return std::make_tuple(begin, end);
    }

    std::tuple<const_iterator, const_iterator> range(int i) const noexcept
    {
        bx_expects(i >= 0 && i < m_access_length);

        iterator begin = m_values.get() + m_access[i];
        iterator end = m_values.get() + m_access[i + 1];

        return std::make_tuple(begin, end);
    }

    int values_size() const noexcept
    {
        return m_values_length;
    }
};
}

#endif
