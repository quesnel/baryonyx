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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_DYNARRAY_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_DYNARRAY_HPP

#include <vector>

namespace baryonyx {

template<typename T>
class dynamic_array_view;

template<typename T>
class dynamic_array
{
public:
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
    using value_type = typename std::vector<T>::value_type;
    using size_type = typename std::vector<T>::size_type;
    using index_type = typename std::vector<T>::difference_type;

private:
    constexpr static size_type size_max{ 512 };

    std::vector<T> m_data;
    size_type m_line_size = { 0 };
    size_type m_capacity = { 0 };
    size_type m_size = { 0 };

    constexpr size_type compute_capacity(size_type rows)
    {
        return (1 + (rows / size_max)) * size_max;
    }

public:
    dynamic_array() = default;
    ~dynamic_array() = default;

    iterator begin() noexcept;
    const_iterator begin() const noexcept;
    const_iterator cbegin() const noexcept;
    iterator end() noexcept;
    const_iterator end() const noexcept;
    const_iterator cend() const noexcept;

    void init(size_type cols)
    {
        m_line_size = cols;
        m_capacity = size_max;
        m_size = 0;

        m_data.resize(m_capacity * m_line_size);
    }

    void init(size_type rows, size_type cols)
    {
        m_capacity = compute_capacity(rows);
        m_line_size = cols;
        m_size = rows;

        m_data.resize(m_capacity * m_line_size);
    }

    void push_line()
    {
        ++m_size;

        if (m_size > m_capacity) {
            m_capacity *= 3;
            m_capacity /= 2;
            m_data.resize(m_capacity * m_line_size);
        }
    }

    void pop_line()
    {
        if (m_size)
            --m_size;
    }

    T operator()(index_type row, index_type col) const
    {
        return m_data[row * m_line_size + col];
    }

    T& operator()(index_type row, index_type col)
    {
        return m_data[row * m_line_size + col];
    }

    T at(index_type row, index_type col) const
    {
        return m_data.at(row * m_line_size + col);
    }

    T& at(index_type row, index_type col)
    {
        return m_data.at(row * m_line_size + col);
    }

    void swap(dynamic_array<T>& other)
    {
        using std::swap;

        swap(m_data, other.m_data);
        swap(m_line_size, other.m_line_size);
        swap(m_capacity, other.m_capacity);
        swap(m_size, other.m_size);
    }

    dynamic_array_view<T> row(index_type row);
    dynamic_array_view<T> row(index_type row) const;

    size_t rows() const
    {
        return m_size;
    }

    size_t cols() const
    {
        return m_line_size;
    }
};

template<typename T>
class dynamic_array_view
{
public:
    using iterator = typename dynamic_array<T>::iterator;
    using const_iterator = typename dynamic_array<T>::const_iterator;
    using value_type = typename dynamic_array<T>::value_type;
    using size_type = typename dynamic_array<T>::size_type;
    using index_type = typename dynamic_array<T>::difference_type;

private:
    const_iterator m_first;
    const_iterator m_last;

public:
    dynamic_array_view(const_iterator first, const_iterator last)
      : m_first(first)
      , m_last(last)
    {}

    template<typename Integer>
    value_type operator[](Integer i) const
    {
        return *(m_first + i);
    }

    template<typename Integer>
    value_type operator[](Integer i)
    {
        return *(m_first + i);
    }
};

template<typename T>
typename dynamic_array<T>::iterator
dynamic_array<T>::begin() noexcept
{
    return m_data.begin();
}

template<typename T>
typename dynamic_array<T>::const_iterator
dynamic_array<T>::begin() const noexcept
{
    return m_data.begin();
}

template<typename T>
typename dynamic_array<T>::const_iterator
dynamic_array<T>::cbegin() const noexcept
{
    return m_data.cbegin();
}

template<typename T>
typename dynamic_array<T>::iterator
dynamic_array<T>::end() noexcept
{
    return m_data.end();
}

template<typename T>
typename dynamic_array<T>::const_iterator
dynamic_array<T>::end() const noexcept
{
    return m_data.end();
}

template<typename T>
typename dynamic_array<T>::const_iterator
dynamic_array<T>::cend() const noexcept
{
    return m_data.cend();
}

template<typename T>
dynamic_array_view<T>
dynamic_array<T>::row(index_type row)
{
    return dynamic_array_view<T>(m_data.begin() + (row * m_line_size),
                                 m_data.begin() + ((row + 1) * m_line_size));
}

template<typename T>
dynamic_array_view<T>
dynamic_array<T>::row(index_type row) const
{
    return dynamic_array_view<T>(m_data.begin() + (row * m_line_size),
                                 m_data.begin() + ((row + 1) * m_line_size));
}

} // namespace efyj

#endif
