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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_MATRIX_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_MATRIX_HPP

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <vector>

#include "fixed_array.hpp"
#include "utils.hpp"

#include <iomanip>

#include <cassert>

namespace baryonyx {

template<typename A_T, typename P_T>
class SparseArray;

template<typename A_T, typename P_T>
std::ostream&
operator<<(std::ostream& os, const SparseArray<A_T, P_T>& m);

/**
 * @brief A specfic class to store A and P sparse matrices.
 *
 * The @c SparseArray defined a two-dimensional template array to access A
 *     and P matrices. A and P data are stored into a @e std::vector<T>
 *     and accessors to the sparse array use @e std::vector<access>.
 *
 * @tparam A_T Type of element of matrix A.
 * @tparam P_T Type of element of matrix P.
 */
template<typename A_T, typename P_T>
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

        index_type position{ -1 };
        index_type value{ -1 };
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

    using size_type = std::size_t;
    using iterator = typename fixed_array<access>::iterator;
    using const_iterator = typename fixed_array<access>::const_iterator;

protected:
    fixed_array<int> m_rows_access;
    fixed_array<int> m_cols_access;
    fixed_array<access> m_rows;
    fixed_array<access> m_cols;

    std::vector<a_type> m_a;
    std::vector<p_type> m_p;

public:
    explicit SparseArray(index_type rows, index_type cols);

    SparseArray() = delete;

    ~SparseArray() = default;

    SparseArray(const SparseArray& q) = default;
    SparseArray(SparseArray&& q) = default;

    SparseArray& operator=(const SparseArray& q) = default;
    SparseArray& operator=(SparseArray&& q) = default;

    bool empty() const noexcept;
    size_type size() const noexcept;

    std::tuple<iterator, iterator> row(index_type row) noexcept;
    std::tuple<iterator, iterator> column(index_type col) noexcept;
    std::tuple<const_iterator, const_iterator> row(index_type row) const
      noexcept;
    std::tuple<const_iterator, const_iterator> column(index_type col) const
      noexcept;

    void set(index_type row, index_type col, a_type x, p_type y);
    void set_p(index_type row, index_type col, p_type y);

    void add_a(index_type row, index_type col, a_type x);
    void add_p(index_type row, index_type col, p_type y);

    void invert_p(index_type row, index_type col);
    void mult_p(index_type row, index_type col, p_type y);

    void mult_row_p(index_type row, p_type y);

    template<typename InputIterator>
    void reserve(index_type elem,
                 InputIterator row_begin,
                 InputIterator row_end,
                 InputIterator col_begin,
                 InputIterator col_end);

    void sort() noexcept;

    a_type A(index_type row, index_type col) const;
    p_type P(index_type row, index_type col) const;

    std::vector<a_type>& A();
    const std::vector<a_type>& A() const;
    std::vector<p_type>& P();
    const std::vector<p_type>& P() const;

    void swap(SparseArray& c) noexcept(
      noexcept(swap(m_a, c.m_a)) && noexcept(swap(m_p, c.m_p)) &&
      noexcept(swap(m_cols_access, c.m_cols_acces)) &&
      noexcept(swap(m_rows_access, c.m_rows_acces)));

    friend std::ostream& operator<<<>(std::ostream& os, const SparseArray& m);

private:
    void m_check_index(index_type row, index_type col) const noexcept;

    inline index_type binary_find(index_type row, index_type col) const;
};

template<typename A_T, typename P_T>
std::ostream&
operator<<(std::ostream& os, const SparseArray<A_T, P_T>& m)
{
    const std::size_t rows = m.m_rows_access.size();
    const std::size_t cols = m.m_cols_access.size();
    std::vector<P_T> to_show(cols);

    for (std::size_t i{ 0 }; i != rows; ++i) {
        std::fill(std::begin(to_show),
                  std::end(to_show),
                  std::numeric_limits<double>::infinity());

        auto its = m.row(i);

        for (; std::get<0>(its) != std::get<1>(its); ++std::get<0>(its))
            to_show[std::get<0>(its)->position] =
              m.m_p[std::get<0>(its)->value];

        for (auto elem : to_show)
            if (elem == std::numeric_limits<double>::infinity())
                os << std::setw(8) << std::left << "        ";
            else
                os << std::setw(8) << std::right << elem;

        os << '\n';
    }

    return os;
}

template<typename A_T, typename P_T>
SparseArray<A_T, P_T>::SparseArray(index_type rows, index_type cols)
  : m_rows_access(rows, -1)
  , m_cols_access(cols, -1)
{
}

template<typename A_T, typename P_T>
bool
SparseArray<A_T, P_T>::empty() const noexcept
{
    return m_a.empty();
}

template<typename A_T, typename P_T>
typename SparseArray<A_T, P_T>::size_type
SparseArray<A_T, P_T>::size() const noexcept
{
    return m_a.size();
}

template<typename A_T, typename P_T>
std::tuple<typename SparseArray<A_T, P_T>::iterator,
           typename SparseArray<A_T, P_T>::iterator>
SparseArray<A_T, P_T>::row(index_type row) noexcept
{
    m_check_index(row, 0);

    iterator begin = m_rows.begin() + m_rows_access[row];

    if (static_cast<std::size_t>(row + 1) < m_rows_access.size())
        return std::make_tuple(
          begin, begin + m_rows_access[row + 1] - m_rows_access[row]);
    else
        return std::make_tuple(begin, m_rows.end());
}

template<typename A_T, typename P_T>
std::tuple<typename SparseArray<A_T, P_T>::iterator,
           typename SparseArray<A_T, P_T>::iterator>
SparseArray<A_T, P_T>::column(index_type col) noexcept
{
    m_check_index(0, col);

    iterator begin = m_cols.begin() + m_cols_access[col];

    if (static_cast<std::size_t>(col + 1) < m_cols_access.size())
        return std::make_tuple(
          begin, begin + m_cols_access[col + 1] - m_cols_access[col]);
    else
        return std::make_tuple(begin, m_cols.end());
}

template<typename A_T, typename P_T>
std::tuple<typename SparseArray<A_T, P_T>::const_iterator,
           typename SparseArray<A_T, P_T>::const_iterator>
SparseArray<A_T, P_T>::row(index_type row) const noexcept
{
    m_check_index(row, 0);

    const_iterator begin = std::next(m_rows.begin(), m_rows_access[row]);

    if (static_cast<std::size_t>(row + 1) < m_rows_access.size())
        return std::make_tuple(
          begin, begin + m_rows_access[row + 1] - m_rows_access[row]);
    else
        return std::make_tuple(begin, m_rows.end());
}

template<typename A_T, typename P_T>
std::tuple<typename SparseArray<A_T, P_T>::const_iterator,
           typename SparseArray<A_T, P_T>::const_iterator>
SparseArray<A_T, P_T>::column(index_type col) const noexcept
{
    m_check_index(0, col);

    const_iterator begin = m_cols.begin() + m_cols_access[col];

    if (static_cast<std::size_t>(col + 1) < m_cols_access.size())
        return std::make_tuple(
          begin, begin + m_cols_access[col + 1] - m_cols_access[col]);
    else
        return std::make_tuple(begin, m_cols.end());
}

template<typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::set(index_type row, index_type col, a_type x, p_type y)
{
    m_check_index(row, col);

    auto id = m_a.size();

    for (auto i = SparseArray<A_T, P_T>::row(row);
         std::get<0>(i) != std::get<1>(i);
         ++std::get<0>(i)) {
        if (std::get<0>(i)->position == -1) {
            std::get<0>(i)->position = col;
            std::get<0>(i)->value = id;
            m_a.emplace_back(x);
            m_p.emplace_back(y);
            break;
        }
    }

    for (auto i = SparseArray<A_T, P_T>::column(col);
         std::get<0>(i) != std::get<1>(i);
         ++std::get<0>(i)) {
        if (std::get<0>(i)->position == -1) {
            std::get<0>(i)->position = row;
            std::get<0>(i)->value = id;
            break;
        }
    }
}

template<typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::set_p(index_type row, index_type col, p_type y)
{
    m_check_index(row, col);

    auto it = binary_find(row, col);

    m_p[it] = y;
}

template<typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::add_a(index_type row, index_type col, a_type x)
{
    m_check_index(row, col);

    auto it = binary_find(row, col);

    m_a[it] += x;
}

template<typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::add_p(index_type row, index_type col, p_type y)
{
    m_check_index(row, col);

    auto it = binary_find(row, col);

    m_p[it] += y;
}

template<typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::invert_p(index_type row, index_type col)
{
    m_check_index(row, col);

    auto it = binary_find(row, col);

    m_p[it] = -m_p[it];
}

template<typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::mult_p(index_type row, index_type col, p_type y)
{
    m_check_index(row, col);

    auto it = binary_find(row, col);

    m_p[it] *= y;
}

template<typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::mult_row_p(index_type row, p_type y)
{
    auto i = SparseArray<A_T, P_T>::row(row);

    for (; std::get<0>(i) != std::get<1>(i); ++std::get<0>(i))
        m_p[std::get<0>(i)->value] *= y;
}

template<typename A_T, typename P_T>
template<typename InputIterator>
void
SparseArray<A_T, P_T>::reserve(index_type elem,
                               InputIterator row_begin,
                               InputIterator row_end,
                               InputIterator col_begin,
                               InputIterator col_end)
{
    m_a.reserve(elem);
    m_p.reserve(elem);

    fixed_array<access>(elem).swap(m_rows);
    fixed_array<access>(elem).swap(m_cols);

    for (index i{ 0 }, current{ 0 }; row_begin != row_end; ++row_begin, ++i) {
        m_rows_access[i] = current;
        current += *row_begin;
    }

    for (index i{ 0 }, current{ 0 }; col_begin != col_end; ++col_begin, ++i) {
        m_cols_access[i] = current;
        current += *col_begin;
    }
}

template<typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::sort() noexcept
{
    if (m_rows_access.size() <= 1)
        return;

    for (index i{ 0 }, e{ static_cast<index>(m_rows_access.size() - 1) };
         i < e;
         ++i) {
        std::sort(&m_rows[m_rows_access[i]],
                  &m_rows[m_rows_access[i + 1]],
                  [](const auto& lhs, const auto& rhs) {
                      return lhs.position < rhs.position;
                  });
    }

    if (m_cols_access.size() <= 1)
        return;

    for (std::size_t i{ 0 }, e{ m_cols_access.size() - 1 }; i != e; ++i) {
        std::sort(&m_cols[m_cols_access[i]],
                  &m_cols[m_cols_access[i + 1]],
                  [](const auto& lhs, const auto& rhs) {
                      return lhs.position < rhs.position;
                  });
    }
}

template<typename A_T, typename P_T>
typename SparseArray<A_T, P_T>::a_type
SparseArray<A_T, P_T>::A(index_type row, index_type col) const
{
    m_check_index(row, col);

    auto it = binary_find(row, col);

    return m_a[it];
}

template<typename A_T, typename P_T>
typename SparseArray<A_T, P_T>::p_type
SparseArray<A_T, P_T>::P(index_type row, index_type col) const
{
    m_check_index(row, col);

    auto it = binary_find(row, col);

    return m_p[it];
}

template<typename A_T, typename P_T>
typename std::vector<typename SparseArray<A_T, P_T>::a_type>&
SparseArray<A_T, P_T>::A()
{
    return m_a;
}

template<typename A_T, typename P_T>
const typename std::vector<typename SparseArray<A_T, P_T>::a_type>&
SparseArray<A_T, P_T>::A() const
{
    return m_a;
}

template<typename A_T, typename P_T>
typename std::vector<typename SparseArray<A_T, P_T>::p_type>&
SparseArray<A_T, P_T>::P()
{
    return m_p;
}

template<typename A_T, typename P_T>
const typename std::vector<typename SparseArray<A_T, P_T>::p_type>&
SparseArray<A_T, P_T>::P() const
{
    return m_p;
}

template<typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::swap(SparseArray& c) noexcept(
  noexcept(swap(m_a, c.m_a)) && noexcept(swap(m_p, c.m_p)) &&
  noexcept(swap(m_cols_access, c.m_cols_acces)) &&
  noexcept(swap(m_rows_access, c.m_rows_acces)))
{
    std::swap(m_a, c.m_a);
    std::swap(m_p, c.m_p);
    std::swap(m_cols_access, c.m_cols_access);
    std::swap(m_rows_access, c.m_rows_access);
}

template<typename A_T, typename P_T>
void
SparseArray<A_T, P_T>::m_check_index(index_type row, index_type col) const
  noexcept
{
    (void)row;
    (void)col;
    assert(col >= 0 && "SparseArray: bad column access");
    assert(row >= 0 && "SparseArray: bad row access");
}

template<typename A_T, typename P_T>
typename SparseArray<A_T, P_T>::index_type
SparseArray<A_T, P_T>::binary_find(
  typename SparseArray<A_T, P_T>::index_type row,
  typename SparseArray<A_T, P_T>::index_type col) const
{

    auto i = SparseArray<A_T, P_T>::row(row);
    auto ret =
      std::lower_bound(std::get<0>(i), std::get<1>(i), col, access_compare());

#ifndef BARYONYX_FULL_OPTIMIZATION
    if (!(ret != std::get<1>(i) && ret->position == col))
        throw std::out_of_range("SparseArray::binary_find");
#endif

    return ret->value;
}

} // namespace baryonyx

#endif
