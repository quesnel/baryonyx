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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_SPARSE_MATRIX_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_SPARSE_MATRIX_HPP

#include <baryonyx/core>

#include <algorithm>
#include <memory>
#include <string>

#include "debug.hpp"
#include "fixed-array.hpp"

namespace baryonyx {

template<typename Index>
class sparse_matrix
{
public:
    using size_type = std::size_t;
    using index_type = Index;

    struct col_value
    {
        col_value() = default;

        col_value(int value_, int row_)
          : value(value_)
          , row(row_)
        {}

        int value = -1;
        int row = -1;
    };

    struct row_value
    {
        row_value() = default;

        row_value(int value_, int column_)
          : value(value_)
          , column(column_)
        {}

        int value = -1;
        int column = -1;
    };

protected:
    using vector_access = fixed_array<int>;
    using row_vector_access = fixed_array<row_value>;
    using col_vector_access = fixed_array<col_value>;

    vector_access m_rows_access;
    vector_access m_cols_access;
    row_vector_access m_rows_values;
    col_vector_access m_cols_values;

public:
    using row_iterator = typename row_vector_access::iterator;
    using col_iterator = typename col_vector_access::iterator;
    using const_row_iterator = typename row_vector_access::const_iterator;
    using const_col_iterator = typename col_vector_access::const_iterator;

    template<typename Constraint>
    explicit sparse_matrix(const std::vector<Constraint>& csts,
                           const int rows,
                           const int cols)
      : m_rows_access(rows + 1)
      , m_cols_access(cols + 1)
    {
        int elem = 0;

        for (int i = 0, e = rows; i != e; ++i)
            elem += static_cast<int>(csts[i].elements.size());

        row_vector_access(elem).swap(m_rows_values);
        col_vector_access(elem).swap(m_cols_values);

        //
        // First, we build a vector of access to all elements in the matrix
        // different to 0. Then, we sort accessors according to the column id
        // then we assign a unique identifier for each accessors and then we
        // sort accessors according to the row id. This code enable the
        // construction of the r_rows_values with pointers to the
        // r_cols_values.
        //

        struct access
        {
            access() = default;

            access(int row_, int col_, function_element fct_)
              : row(row_)
              , col(col_)
              , fct(fct_)
            {}

            int row;
            int col;
            int id;
            function_element fct;
        };

        std::vector<access> accessors;
        accessors.reserve(elem);
        for (int i = 0; i != rows; ++i)
            for (const auto& fct_elem : csts[i].elements)
                accessors.emplace_back(i, fct_elem.variable_index, fct_elem);

        fixed_array<int> rinit(rows, 0);
        int row = accessors[0].row;
        int nb_row = 0;
        for (int i = 0; i != elem; ++i) {
            if (row != accessors[i].row) {
                rinit[row] = nb_row;
                nb_row = 0;
                row = accessors[i].row;
            }
            ++nb_row;

            accessors[i].id = i;
            m_rows_values[i] = { i, accessors[i].col };
        }

        rinit[row] = nb_row; /* Be sure to affect the good number of elements
                                for this constaints. */

        for (int i = 0; i != rows; ++i)
            bx_assert(rinit[i] == static_cast<int>(csts[i].elements.size()));

        // We sort the accessors to improve memory friendly access.

        std::stable_sort(
          accessors.begin(),
          accessors.end(),
          [](const auto& lhs, const auto& rhs) { return lhs.col < rhs.col; });

        fixed_array<int> cinit(cols, 0);
        int col = accessors[0].col;
        int nb_col = 0;

        for (int i = 0; i != elem; ++i) {
            if (col != accessors[i].col) {
                cinit[col] = nb_col;
                nb_col = 0;
                col = accessors[i].col;
            }
            ++nb_col;

            m_cols_values[i] = { m_rows_values[accessors[i].id].value,
                                 accessors[i].row };
        }

        cinit[col] = nb_col; /* Be sure to affect the good number of elements
                                for this constaints. */

        m_rows_access[0] = 0;
        for (int i = 0; i != rows; ++i)
            m_rows_access[i + 1] = m_rows_access[i] + rinit[i];

        m_cols_access[0] = 0;
        for (int i = 0; i != cols; ++i)
            m_cols_access[i + 1] = m_cols_access[i] + cinit[i];
    }

    std::tuple<row_iterator, row_iterator> row(int row) noexcept
    {
        m_check_index(row, 0);

        row_iterator begin = m_rows_values.begin() + m_rows_access[row];
        row_iterator end = m_rows_values.begin() + m_rows_access[row + 1];

        return std::make_tuple(begin, end);
    }

    std::tuple<col_iterator, col_iterator> column(int col) noexcept
    {
        m_check_index(0, col);

        col_iterator begin = m_cols_values.begin() + m_cols_access[col];
        col_iterator end = m_cols_values.begin() + m_cols_access[col + 1];

        return std::make_tuple(begin, end);
    }

    std::tuple<const_row_iterator, const_row_iterator> row(int row) const
      noexcept
    {
        m_check_index(row, 0);

        const_row_iterator begin = m_rows_values.begin() + m_rows_access[row];
        const_row_iterator end =
          m_rows_values.begin() + m_rows_access[row + 1];

        return std::make_tuple(begin, end);
    }

    std::tuple<const_col_iterator, const_col_iterator> column(int col) const
      noexcept
    {
        m_check_index(0, col);

        const_col_iterator begin = m_cols_values.begin() + m_cols_access[col];
        const_col_iterator end =
          m_cols_values.begin() + m_cols_access[col + 1];

        return std::make_tuple(begin, end);
    }

    size_type size() const noexcept
    {
        return m_rows_values.size();
    }

    int length() const noexcept
    {
        return static_cast<int>(size());
    }

private:
    void m_check_index(index_type row, index_type col) const noexcept
    {
        bx_expects(col >= 0 &&
                   static_cast<size_type>(col) < (m_cols_access.size() - 1));

        bx_expects(row >= 0 &&
                   static_cast<size_type>(row) < (m_rows_access.size() - 1));
    }
};
}

#endif
