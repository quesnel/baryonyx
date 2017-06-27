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
 * stored into a \e std::vector<T> and accessors to the sparse array use \e
 * std::vector<access>.
 *
 * \tparam A_T Type of element of matrix A.
 * \tparam P_T Type of element of matrix P.
 */
template <typename A_T, typename P_T>
class SparseArray
{
public:
    static_assert(std::is_integral<A_T>::value, "A_T Integer required.");
    static_assert(std::is_floating_point<P_T>::value, "P_T Real required.");

    using a_type = A_T;
    using p_type = P_T;
    using index_type = int;

    struct access
    {
        access() = default;

        access(index_type position_, index_type value_)
          : position(position_)
          , value(value_)
        {
        }

        index_type position;
        index_type value;
    };

    struct access_compare
    {
        bool operator()(const access& s, index_type i)
        {
            return s.position < i;
        }

        bool operator()(index_type i, const access& s)
        {
            return i < s.position;
        }
    };

    using accessor_type = std::vector<access>;
    using size_type = typename accessor_type::size_type;

protected:
    std::vector<accessor_type> m_rows;
    std::vector<accessor_type> m_cols;
    std::vector<a_type> m_a;
    std::vector<p_type> m_p;

    // container_type m_values;

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

    void set(index_type row, index_type col, a_type x, p_type y);
    void set_p(index_type row, index_type col, p_type y);

    void add_a(index_type row, index_type col, a_type x);
    void add_p(index_type row, index_type col, p_type y);

    void invert(index_type row, index_type col);
    void mult_p(index_type row, index_type col, p_type y);

    void mult_row_p(index_type row, p_type y);

    void sort() noexcept;

    a_type A(index_type row, index_type col) const;
    p_type P(index_type row, index_type col) const;

    std::vector<a_type>& A();
    const std::vector<a_type>& A() const;
    std::vector<p_type>& P();
    const std::vector<p_type>& P() const;

    void swap(SparseArray& c) noexcept(noexcept(swap(m_a, c.m_a)) &&
                                       noexcept(swap(m_p, c.m_p)) &&
                                       noexcept(swap(m_cols, c.m_cols)) &&
                                       noexcept(swap(m_rows, c.m_rows)));

private:
    void m_check_index(index_type row, index_type col) const noexcept;

    inline index_type binary_find(index_type row, index_type col) const;
};

template <typename A_T, typename P_T>
SparseArray<A_T, P_T>::SparseArray()
{
}

template <typename A_T, typename P_T>
SparseArray<A_T, P_T>::SparseArray(index_type rows, index_type cols)
  : m_rows(rows)
  , m_cols(cols)
{
}

template <typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::resize(index_type rows, index_type cols)
{
    m_cols.clear();
    m_cols.resize(cols);

    m_rows.clear();
    m_rows.resize(rows);

    m_a.clear();
    m_p.clear();
}

template <typename A_T, typename P_T>
bool
SparseArray<A_T, P_T>::empty() const noexcept
{
    return m_a.empty();
}

template <typename A_T, typename P_T>
typename SparseArray<A_T, P_T>::size_type
SparseArray<A_T, P_T>::size() const noexcept
{
    return m_a.size();
}

template <typename A_T, typename P_T>
typename SparseArray<A_T, P_T>::size_type
SparseArray<A_T, P_T>::rows() const noexcept
{
    return m_rows.size();
}

template <typename A_T, typename P_T>
typename SparseArray<A_T, P_T>::size_type
SparseArray<A_T, P_T>::columns() const noexcept
{
    return m_cols.size();
}

template <typename A_T, typename P_T>
const typename SparseArray<A_T, P_T>::accessor_type&
SparseArray<A_T, P_T>::row(index_type row) const noexcept
{
    m_check_index(row, 0);

    return m_rows[row];
}

template <typename A_T, typename P_T>
const typename SparseArray<A_T, P_T>::accessor_type&
SparseArray<A_T, P_T>::column(index_type col) const noexcept
{
    m_check_index(0, col);

    return m_cols[col];
}

template <typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::set(index_type row, index_type col, a_type x, p_type y)
{
    m_check_index(row, col);

    auto it = std::find_if(
      m_rows[row].begin(), m_rows[row].end(), [col](const auto& access) {
          return access.position == col;
      });

    if (it == m_rows[row].end()) {
        m_rows[row].emplace_back(col, m_a.size());
        m_cols[col].emplace_back(row, m_a.size());
        m_a.emplace_back(x);
        m_p.emplace_back(y);
    } else {
        m_a[it->value] = x;
        m_p[it->value] = y;
    }
}

template <typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::set_p(index_type row, index_type col, p_type y)
{
    m_check_index(row, col);

    auto it = binary_find(row, col);

    m_p[it] = y;
}

template <typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::add_a(index_type row, index_type col, a_type x)
{
    m_check_index(row, col);

    auto it = binary_find(row, col);

    m_a[it] += x;
}

template <typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::add_p(index_type row, index_type col, p_type y)
{
    m_check_index(row, col);

    auto it = binary_find(row, col);

    m_p[it] += y;
}

template <typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::invert(index_type row, index_type col)
{
    m_check_index(row, col);

    auto it = binary_find(row, col);

    m_a[it] = -m_a[it];
    m_p[it] = -m_p[it];
}

template <typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::mult_p(index_type row, index_type col, p_type y)
{
    m_check_index(row, col);

    auto it = binary_find(row, col);

    m_p[it] *= y;
}

template <typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::mult_row_p(index_type row, p_type y)
{
    auto it = m_rows[row].cbegin();
    auto et = m_rows[row].cend();

    for (; it != et; ++it)
        m_p[it->value] *= y;
}

template <typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::sort() noexcept
{
    for (auto& v : m_rows) {
        std::sort(v.begin(), v.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.position < rhs.position;
        });
    }

    for (auto& v : m_cols) {
        std::sort(v.begin(), v.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.position < rhs.position;
        });
    }
}

template <typename A_T, typename P_T>
typename SparseArray<A_T, P_T>::a_type
SparseArray<A_T, P_T>::A(index_type row, index_type col) const
{
    m_check_index(row, col);

    auto it = binary_find(row, col);

    return m_a[it];
}

template <typename A_T, typename P_T>
typename SparseArray<A_T, P_T>::p_type
SparseArray<A_T, P_T>::P(index_type row, index_type col) const
{
    m_check_index(row, col);

    auto it = binary_find(row, col);

    return m_p[it];
}

template <typename A_T, typename P_T>
typename std::vector<typename SparseArray<A_T, P_T>::a_type>&
SparseArray<A_T, P_T>::A()
{
    return m_a;
}

template <typename A_T, typename P_T>
const typename std::vector<typename SparseArray<A_T, P_T>::a_type>&
SparseArray<A_T, P_T>::A() const
{
    return m_a;
}

template <typename A_T, typename P_T>
typename std::vector<typename SparseArray<A_T, P_T>::p_type>&
SparseArray<A_T, P_T>::P()
{
    return m_p;
}

template <typename A_T, typename P_T>
const typename std::vector<typename SparseArray<A_T, P_T>::p_type>&
SparseArray<A_T, P_T>::P() const
{
    return m_p;
}

template <typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::swap(SparseArray& c) noexcept(
  noexcept(swap(m_a, c.m_a)) && noexcept(swap(m_p, c.m_p)) &&
  noexcept(swap(m_cols, c.m_cols)) &&
  noexcept(swap(m_rows, c.m_rows)))
{
    std::swap(m_a, c.m_a);
    std::swap(m_p, c.m_p);
    std::swap(m_cols, c.m_cols);
    std::swap(m_rows, c.m_rows);
}

#ifndef LP_FULL_OPTIMIZATION
template <typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::m_check_index(index_type row, index_type col) const
  noexcept
{
    Expects(col >= 0 and numeric_cast<std::size_t>(col) < m_cols.size(),
            "SparseArray: bad column access");

    Expects(row >= 0 and numeric_cast<std::size_t>(row) < m_rows.size(),
            "SparseArray: bad row access");
}
#else
template <typename A_T, typename P_T>
void SparseArray<A_T, P_T>::m_check_index(index_type, index_type) const
  noexcept
{
}
#endif

template <typename A_T, typename P_T>
typename SparseArray<A_T, P_T>::index_type
SparseArray<A_T, P_T>::binary_find(
  typename SparseArray<A_T, P_T>::index_type row,
  typename SparseArray<A_T, P_T>::index_type col) const
{
    if (m_rows[row].size() < m_cols[col].size()) {
        auto it = std::lower_bound(
          m_rows[row].begin(), m_rows[row].end(), col, access_compare());

#ifndef LP_FULL_OPTIMIZATION
        if (!(it != m_rows[row].end() && it->position == col))
            throw std::out_of_range("SparseArray::binary_find");
#endif

        return it->value;
    } else {
        auto it = std::lower_bound(
          m_cols[col].begin(), m_cols[col].end(), row, access_compare());

#ifndef LP_FULL_OPTIMIZATION
        if (!(it != m_cols[col].end() && it->position == row))
            throw std::out_of_range("SparseArray::binary_find");
#endif

        return it->value;
    }
}

} // namespace lp

#endif
