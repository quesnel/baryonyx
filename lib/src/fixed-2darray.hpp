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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_FIXED_2DARRAY_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_FIXED_2DARRAY_HPP

#include <algorithm>
#include <iterator>
#include <memory>

#include <cassert>
#include <cstddef>

namespace baryonyx {

/**
 * @brief A two-dimensional template array.
 *
 * @details @c fixed_2array stores a smart pointer to a dynamically allocated
 *     array (with @c std::make_unique<T[]>). @c fixed_2darray garantes to
 *     delete the pointer on destruction of the @c fixed_2darray. The rows and
 *     columns are stored into this container so, @c fixed_2darray provides a
 *     std::vector's like API.
 *
 * @tparam T Type of element.
 */
template<typename T>
class fixed_2darray
{
public:
    using value_type = T;
    using reference = T&;
    using const_reference = const T&;
    using iterator = T*;
    using const_iterator = const T*;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

protected:
    std::size_t m_rows;
    std::size_t m_columns;
    std::unique_ptr<T[]> m_buffer;

public:
    /**
     * @brief Constructs and empty container without any elements.
     */
    fixed_2darray() noexcept;

    explicit fixed_2darray(size_type rows, size_type columns);
    explicit fixed_2darray(size_type rows,
                           size_type columns,
                           const value_type& def);

    ~fixed_2darray() = default;

    fixed_2darray(const fixed_2darray& q);
    fixed_2darray& operator=(const fixed_2darray& q);

    fixed_2darray(fixed_2darray&& q) = default;
    fixed_2darray& operator=(fixed_2darray&& q) = default;

    template<class InputIterator>
    void assign(InputIterator first, InputIterator last);

    explicit operator bool() const noexcept;
    bool empty() const noexcept;
    size_type size() const noexcept;
    size_type rows() const noexcept;
    size_type columns() const noexcept;

    iterator begin() noexcept;
    const_iterator begin() const noexcept;
    const_iterator cbegin() const noexcept;
    iterator end() noexcept;
    const_iterator end() const noexcept;
    const_iterator cend() const noexcept;

    reverse_iterator rbegin() noexcept;
    const_reverse_iterator rbegin() const noexcept;
    const_reverse_iterator crbegin() const noexcept;
    reverse_iterator rend() noexcept;
    const_reverse_iterator rend() const noexcept;
    const_reverse_iterator crend() const noexcept;

    T& front() noexcept;
    T const& front() const noexcept;
    T& back() noexcept;
    const T& back() const noexcept;
    T* data() noexcept;
    const T* data() const noexcept;

    void set(size_type row, size_type col, const value_type& x);
    void set(size_type row, size_type col, value_type&& x);

    template<class... Args>
    void emplace(size_type row, size_type col, Args&&... args);

    const_reference operator()(size_type row, size_type col) const;
    reference operator()(size_type row, size_type col);

    void swap(fixed_2darray& c) noexcept(noexcept(m_buffer.swap(c.m_buffer)));

private:
    bool m_is_valid(size_type row, size_type col) const noexcept;
};

template<typename T>
fixed_2darray<T>::fixed_2darray() noexcept
  : m_rows{ 0 }
  , m_columns{ 0 }
{}

template<typename T>
fixed_2darray<T>::fixed_2darray(size_type rows, size_type columns)
  : m_rows{ rows }
  , m_columns{ columns }
  , m_buffer{ std::make_unique<T[]>(columns * rows) }
{}

template<typename T>
fixed_2darray<T>::fixed_2darray(size_type rows,
                                size_type columns,
                                const value_type& def)
  : m_rows{ rows }
  , m_columns{ columns }
  , m_buffer{ std::make_unique<T[]>(columns * rows) }
{
    std::fill(begin(), end(), def);
}

template<typename T>
fixed_2darray<T>::fixed_2darray(const fixed_2darray& o)
  : m_rows{ o.m_rows }
  , m_columns{ o.m_columns }
  , m_buffer{ std::make_unique<T[]>(m_columns * m_rows) }
{
    std::copy(o.m_buffer.begin(), o.m_buffer.end(), begin());
}

template<typename T>
fixed_2darray<T>&
fixed_2darray<T>::operator=(const fixed_2darray& o)
{
    auto tmp = std::make_unique<T[]>(o.m_columns * o.m_columns);

    std::copy(o.begin(), o.end(), tmp.get());

    m_columns = o.m_columns;
    m_rows = o.m_rows;
    m_buffer = std::move(tmp);

    return *this;
}

template<typename T>
template<class InputIterator>
void
fixed_2darray<T>::assign(InputIterator first, InputIterator last)
{
    for (size_type r = 0; r != m_rows; ++r) {
        for (size_type c = 0; c != m_columns; ++c) {
            set(r, c, *first++);
            if (first == last)
                return;
        }
    }
}

template<class T>
fixed_2darray<T>::operator bool() const noexcept
{
    return m_buffer.get() != nullptr;
}

template<typename T>
bool
fixed_2darray<T>::empty() const noexcept
{
    return m_columns * m_rows == 0;
}

template<typename T>
typename fixed_2darray<T>::size_type
fixed_2darray<T>::size() const noexcept
{
    return m_columns * m_rows;
}

template<typename T>
typename fixed_2darray<T>::size_type
fixed_2darray<T>::rows() const noexcept
{
    return m_rows;
}

template<typename T>
typename fixed_2darray<T>::size_type
fixed_2darray<T>::columns() const noexcept
{
    return m_columns;
}

template<typename T>
typename fixed_2darray<T>::iterator
fixed_2darray<T>::begin() noexcept
{
    return fixed_2darray<T>::iterator(data());
}

template<typename T>
typename fixed_2darray<T>::const_iterator
fixed_2darray<T>::begin() const noexcept
{
    return fixed_2darray<T>::const_iterator(data());
}

template<typename T>
typename fixed_2darray<T>::const_iterator
fixed_2darray<T>::cbegin() const noexcept
{
    return fixed_2darray<T>::const_iterator(data());
}

template<typename T>
typename fixed_2darray<T>::iterator
fixed_2darray<T>::end() noexcept
{
    return fixed_2darray<T>::iterator(data() + size());
}

template<typename T>
typename fixed_2darray<T>::const_iterator
fixed_2darray<T>::end() const noexcept
{
    return fixed_2darray<T>::const_iterator(data() + size());
}

template<typename T>
typename fixed_2darray<T>::const_iterator
fixed_2darray<T>::cend() const noexcept
{
    return fixed_2darray<T>::const_iterator(data() + size());
}

template<typename T>
typename fixed_2darray<T>::reverse_iterator
fixed_2darray<T>::rbegin() noexcept
{
    return fixed_2darray<T>::reverse_iterator(end());
}

template<typename T>
typename fixed_2darray<T>::const_reverse_iterator
fixed_2darray<T>::rbegin() const noexcept
{
    return fixed_2darray<T>::const_reverse_iterator(end());
}

template<typename T>
typename fixed_2darray<T>::const_reverse_iterator
fixed_2darray<T>::crbegin() const noexcept
{
    return fixed_2darray<T>::const_reverse_iterator(end());
}

template<typename T>
typename fixed_2darray<T>::reverse_iterator
fixed_2darray<T>::rend() noexcept
{
    return fixed_2darray<T>::reverse_iterator(begin());
}

template<typename T>
typename fixed_2darray<T>::const_reverse_iterator
fixed_2darray<T>::rend() const noexcept
{
    return fixed_2darray<T>::const_reverse_iterator(begin());
}

template<typename T>
typename fixed_2darray<T>::const_reverse_iterator
fixed_2darray<T>::crend() const noexcept
{
    return fixed_2darray<T>::const_reverse_iterator(begin());
}

template<typename T>
T&
fixed_2darray<T>::front() noexcept
{
    assert(not empty());

    return *begin();
}

template<typename T>
const T&
fixed_2darray<T>::front() const noexcept
{
    assert(not empty());

    return *begin();
}

template<typename T>
T&
fixed_2darray<T>::back() noexcept
{
    assert(not empty());

    return *(begin() + size() - 1);
}

template<typename T>
const T&
fixed_2darray<T>::back() const noexcept
{
    assert(not empty());

    return *(begin() + size() - 1);
}

template<typename T>
T*
fixed_2darray<T>::data() noexcept
{
    return m_buffer.get();
}

template<typename T>
const T*
fixed_2darray<T>::data() const noexcept
{
    return m_buffer.get();
}

template<typename T>
void
fixed_2darray<T>::set(size_type row, size_type column, const value_type& x)
{
    assert(m_is_valid(row, column));

    m_buffer[row * m_columns + column] = x;
}

template<typename T>
void
fixed_2darray<T>::set(size_type row, size_type column, value_type&& x)
{
    assert(m_is_valid(row, column));

    m_buffer.emplace(std::begin(m_buffer) + (row * m_columns + column),
                     std::move(x));
}

template<typename T>
template<class... Args>
void
fixed_2darray<T>::emplace(size_type row, size_type column, Args&&... args)
{
    assert(m_is_valid(row, column));

    m_buffer.emplace(std::begin(m_buffer) + (row * m_columns + column),
                     std::forward<Args>(args)...);
}

template<typename T>
typename fixed_2darray<T>::const_reference
fixed_2darray<T>::operator()(size_type row, size_type column) const
{
    assert(m_is_valid(row, column));

    return m_buffer[row * m_columns + column];
}

template<typename T>
typename fixed_2darray<T>::reference
fixed_2darray<T>::operator()(size_type row, size_type column)
{
    assert(m_is_valid(row, column));

    return m_buffer[row * m_columns + column];
}

template<typename T>
void
fixed_2darray<T>::swap(fixed_2darray& c) noexcept(
  noexcept(m_buffer.swap(c.m_buffer)))
{
    std::swap(m_buffer, c.m_buffer);
    std::swap(m_columns, c.m_columns);
    std::swap(m_rows, c.m_rows);
}

template<typename T>
bool
fixed_2darray<T>::m_is_valid(size_type row, size_type column) const noexcept
{
    return column < m_columns && row < m_rows;
}
}

#endif
