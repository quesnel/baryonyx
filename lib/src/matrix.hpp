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

#ifndef ORG_VLEPROJECT_LP_MATRIX_HPP
#define ORG_VLEPROJECT_LP_MATRIX_HPP

#include "utils.hpp"
#include <algorithm>
#include <stdexcept>
#include <vector>

namespace lp {

/**
 * An \e SparseArray defined a two-dimensional template array. Informations are
 * stored into a \e std::vector<T> and accessor to the sparse array use \e
 * std::vector<acess>.
 *
 * \tparam T Type of element
 */
template <typename T>
class SparseArray
{
public:
    using value_type = T;
    using index_type = long int;
    using container_type = std::vector<value_type>;
    using reference = typename container_type::reference;
    using const_reference = typename container_type::const_reference;
    using iterator = typename container_type::iterator;
    using const_iterator = typename container_type::const_iterator;
    using reverse_iterator = typename container_type::reverse_iterator;
    using const_reverse_iterator =
      typename container_type::const_reverse_iterator;
    using size_type = typename container_type::size_type;

    struct access
    {
        using value_type = T;

        access() {}
        access(std::size_t position_, std::size_t value_)
          : position(position_)
          , value(value_)
        {
        }

        index_type position;
        index_type value;
    };

    using accessor_type = std::vector<access>;

protected:
    std::vector<accessor_type> m_rows;
    std::vector<accessor_type> m_cols;
    container_type m_values;

public:
    SparseArray();
    explicit SparseArray(index_type rows, index_type cols);

    ~SparseArray() = default;

    SparseArray(const SparseArray& q) = default;
    SparseArray(SparseArray&& q) = default;

    SparseArray& operator=(const SparseArray& q) = default;
    SparseArray& operator=(SparseArray&& q) = default;

    void resize(index_type rows, index_type cols);

    bool empty() const noexcept;
    size_type size() const noexcept;

    size_type rows() const noexcept;
    size_type columns() const noexcept;

    const accessor_type& row(index_type row) const noexcept;
    const accessor_type& column(index_type col) const noexcept;

    void set(index_type row, index_type col, const value_type& x);
    void set(index_type row, index_type col, value_type&& x);

    template <class... Args>
    void emplace(index_type row, index_type col, Args&&... args);

    value_type operator()(index_type row, index_type col) const;

    container_type& values();
    const container_type& values() const;

    void swap(SparseArray& c) noexcept(noexcept(swap(m_values, c.m_values)) &&
                                       noexcept(swap(m_cols, c.m_cols)) &&
                                       noexcept(swap(m_rows, c.m_rows)));

private:
    void m_check_index(index_type row, index_type col) const noexcept;
};

template <typename T>
SparseArray<T>::SparseArray()
{
}

template <typename T>
SparseArray<T>::SparseArray(index_type rows, index_type cols)
  : m_rows(rows)
  , m_cols(cols)
{
}

template <typename T>
void
SparseArray<T>::resize(index_type rows, index_type cols)
{
    m_cols.clear();
    m_cols.resize(cols);

    m_rows.clear();
    m_rows.resize(rows);

    m_values.clear();
}

template <typename T>
bool
SparseArray<T>::empty() const noexcept
{
    return m_values.empty();
}

template <typename T>
typename SparseArray<T>::size_type
SparseArray<T>::size() const noexcept
{
    return m_values.size();
}

template <typename T>
typename SparseArray<T>::size_type
SparseArray<T>::rows() const noexcept
{
    return m_rows.size();
}

template <typename T>
typename SparseArray<T>::size_type
SparseArray<T>::columns() const noexcept
{
    return m_cols.size();
}

template <typename T>
const typename SparseArray<T>::accessor_type&
SparseArray<T>::row(index_type row) const noexcept
{
    m_check_index(row, 0);

    return m_rows[row];
}

template <typename T>
const typename SparseArray<T>::accessor_type&
SparseArray<T>::column(index_type col) const noexcept
{
    m_check_index(0, col);

    return m_cols[col];
}

template <typename T>
void
SparseArray<T>::set(index_type row, index_type col, const value_type& x)
{
    m_check_index(row, col);

    auto it = std::find_if(
      m_rows[row].begin(), m_rows[row].end(), [col](const auto& access) {
          return access.position == col;
      });

    if (it == m_rows[row].end()) {
        m_rows[row].emplace_back(col, m_values.size());
        m_cols[col].emplace_back(row, m_values.size());
        m_values.emplace_back(x);
    } else {
        m_values[it->value] = x;
    }
}

template <typename T>
void
SparseArray<T>::set(index_type row, index_type col, value_type&& x)
{
    m_check_index(row, col);

    auto it = std::find_if(
      m_rows[row].begin(), m_rows[row].end(), [col](const auto& access) {
          return access.position == col;
      });

    if (it == m_rows[row].end()) {
        m_rows[row].emplace_back(col, m_values.size());
        m_cols[col].emplace_back(row, m_values.size());
        m_values.emplace_back(x);
    } else {
        m_values[it->value] = std::move(x);
    }
}

template <typename T>
template <class... Args>
void
SparseArray<T>::emplace(index_type row, index_type col, Args&&... args)
{
    m_check_index(row, col);

    auto it = std::find_if(
        m_rows[row].begin(), m_rows[row].end(), [col](const auto& access) {
          return access.position == col;
      });

    if (it == m_rows[row].end()) {
        m_rows[row].emplace_back(col, m_values.size());
        m_cols[col].emplace_back(row, m_values.size());
        m_values.emplace_back(std::forward<Args>(args)...);
    } else {
        m_values[it->value] = value_type(std::forward<Args>(args)...);
    }
}

template <typename T>
typename SparseArray<T>::value_type
SparseArray<T>::operator()(index_type row, index_type col) const
{
    m_check_index(row, col);

    auto it = std::find_if(
      m_rows[row].begin(), m_rows[row].end(), [col](const auto& access) {
          return access.position == col;
      });

    if (it == m_rows[row].end())
        throw std::out_of_range("SparseArray::operator()");

    return m_values[it->value];
}

template <typename T>
typename SparseArray<T>::container_type&
SparseArray<T>::values()
{
    return m_values;
}

template <typename T>
const typename SparseArray<T>::container_type&
SparseArray<T>::values() const
{
    return m_values;
}

template <typename T>
void
SparseArray<T>::swap(SparseArray& c) noexcept(
  noexcept(swap(m_values, c.m_values)) && noexcept(swap(m_cols, c.m_cols)) &&
  noexcept(swap(m_rows, c.m_rows)))
{
    std::swap(m_values, c.m_values);
    std::swap(m_cols, c.m_cols);
    std::swap(m_rows, c.m_rows);
}

#ifndef NDEBUG
template <typename T>
void
SparseArray<T>::m_check_index(index_type row, index_type col) const noexcept
{
    Expects(col >= 0 and numeric_cast<std::size_t>(col) < m_cols.size(),
            "SparseArray: bad column access");

    Expects(row >= 0 and numeric_cast<std::size_t>(row) < m_rows.size(),
            "SparseArray: bad row access");
}
#else
template <typename T>
void SparseArray<T>::m_check_index(index_type, index_type) const noexcept
{
}
#endif

} // namespace lp

#endif
