/* Copyright (C) 2018-2021 INRAE
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

#include <tuple>
#include <type_traits>
#include <vector>

#include "debug.hpp"

namespace baryonyx {

namespace details {

template<typename Factor>
struct sparse_col_value
{
    int value = -1;
    int row = -1;
    Factor factor{};
};

template<>
struct sparse_col_value<void>
{
    int value = -1;
    int row = -1;
};

} // namespace details

/**
 * Compressed row and column views over the same non-zero cells.
 *
 * When `Factor` is not `void`, the column view also stores the coefficient of
 * each cell so that column scans in the solvers avoid an indirect load.
 */
template<typename Factor = void>
class sparse_matrix
{
public:
    using size_type = std::size_t;
    using col_value = details::sparse_col_value<Factor>;

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
    using vector_access = std::vector<int>;
    using row_vector_access = std::vector<row_value>;
    using col_vector_access = std::vector<col_value>;

    vector_access m_rows_access;
    vector_access m_cols_access;
    row_vector_access m_rows_values;
    col_vector_access m_cols_values;

public:
    using row_iterator = typename row_vector_access::iterator;
    using col_iterator = typename col_vector_access::iterator;
    using const_row_iterator = typename row_vector_access::const_iterator;
    using const_col_iterator = typename col_vector_access::const_iterator;

    /**
     * Builds both the row-major and the column-major views in O(nnz + rows +
     * cols). Inside a column, cells are ordered by increasing row then by
     * position in the row.
     */
    template<typename Constraint>
    explicit sparse_matrix(const std::vector<Constraint>& csts,
                           const int rows,
                           const int cols)
      : m_rows_access(static_cast<std::size_t>(rows) + 1u, 0)
      , m_cols_access(static_cast<std::size_t>(cols) + 1u, 0)
    {
        bx_expects(rows >= 0 && cols >= 0);
        bx_expects(static_cast<size_type>(rows) <= csts.size());

        for (int r = 0; r != rows; ++r) {
            const auto& elements = csts[r].elements;
            m_rows_access[r + 1] =
              m_rows_access[r] + static_cast<int>(elements.size());

            for (const auto& elem : elements) {
                bx_expects(elem.variable_index >= 0 &&
                           elem.variable_index < cols);
                ++m_cols_access[elem.variable_index + 1];
            }
        }

        for (int c = 0; c != cols; ++c)
            m_cols_access[c + 1] += m_cols_access[c];

        const auto nnz = static_cast<std::size_t>(m_rows_access.back());
        m_rows_values.reserve(nnz);
        m_cols_values.resize(nnz);

        std::vector<int> col_cursor(m_cols_access.begin(),
                                    m_cols_access.end() - 1);

        int id = 0;
        for (int r = 0; r != rows; ++r) {
            for (const auto& elem : csts[r].elements) {
                m_rows_values.emplace_back(id, elem.variable_index);
                auto& cell = m_cols_values[col_cursor[elem.variable_index]++];
                cell.value = id;
                cell.row = r;
                if constexpr (!std::is_void_v<Factor>)
                    cell.factor = static_cast<Factor>(elem.factor);
                ++id;
            }
        }
    }

    std::tuple<row_iterator, row_iterator> row(int row) noexcept
    {
        bx_expects(m_is_valid_row(row));

        row_iterator begin = m_rows_values.begin() + m_rows_access[row];
        row_iterator end = m_rows_values.begin() + m_rows_access[row + 1];

        return std::make_tuple(begin, end);
    }

    std::tuple<col_iterator, col_iterator> column(int col) noexcept
    {
        bx_expects(m_is_valid_column(col));

        col_iterator begin = m_cols_values.begin() + m_cols_access[col];
        col_iterator end = m_cols_values.begin() + m_cols_access[col + 1];

        return std::make_tuple(begin, end);
    }

    std::tuple<const_row_iterator, const_row_iterator> row(int row) const
      noexcept
    {
        bx_expects(m_is_valid_row(row));

        const_row_iterator begin = m_rows_values.begin() + m_rows_access[row];
        const_row_iterator end =
          m_rows_values.begin() + m_rows_access[row + 1];

        return std::make_tuple(begin, end);
    }

    std::tuple<const_col_iterator, const_col_iterator> column(int col) const
      noexcept
    {
        bx_expects(m_is_valid_column(col));

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
    bool m_is_valid_row(int row) const noexcept
    {
        return row >= 0 &&
               static_cast<size_type>(row) + 1u < m_rows_access.size();
    }

    bool m_is_valid_column(int col) const noexcept
    {
        return col >= 0 &&
               static_cast<size_type>(col) + 1u < m_cols_access.size();
    }
};
}

#endif
