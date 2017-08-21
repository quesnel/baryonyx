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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_FIXED_ARRAY_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_FIXED_ARRAY_HPP

#include <cstddef>

#include <iterator>
#include <memory>

namespace baryonyx {

/**
 * @brief A @c std::vector<T> array without resize support.
 *
 * @details @c fixed_array stores a smart pointer to a dynamically
 *     allocated array (with @c std::make_unique<T[]>). @c fixed_array
 *     garantes to delete the pointer on destruction of the @c
 *     fixed_array. The length is stored into this container so, @c
 *     fixed_array provides a std::vector's like API.
 */
template<typename T>
class fixed_array
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

private:
    std::size_t m_length;
    std::unique_ptr<T[]> m_buffer;

public:
    /**
     * Constructs an empty container with no elements.
     */
    fixed_array() noexcept;

    /**
     * Constructs a container with @c n elements.
     *
     * @param n Container size (i.e. the number of elements in the
     * container at construction.
     *
     * @exception @c std::bad_alloc.
     */
    explicit fixed_array(std::size_t n);

    /**
     * Constructs a container with @c n elements initialized to @c def.
     *
     * @param n Container size (i.e. the number of elements in the
     * container at construction).
     * @param def Default value to assign for each elements.
     *
     * @exception @c std::bad_alloc.
     */
    fixed_array(std::size_t n, const value_type& def);

    fixed_array(fixed_array&&) = default;
    fixed_array& operator=(fixed_array&&) = default;
    fixed_array(const fixed_array& o);
    fixed_array& operator=(fixed_array const& o);

    ~fixed_array() noexcept = default;

    explicit operator bool() const noexcept;

    std::size_t size() const noexcept;

    bool empty() const noexcept;

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

    T& operator[](std::size_t i) noexcept;
    const T& operator[](std::size_t i) const noexcept;
    T& operator()(std::size_t i) noexcept;
    const T& operator()(std::size_t i) const noexcept;

    void swap(fixed_array& other) noexcept;
};

template<typename T>
fixed_array<T>::fixed_array() noexcept
  : m_length{ 0 }
{
}

template<typename T>
fixed_array<T>::fixed_array(std::size_t n)
  : m_length{ n }
  , m_buffer{ std::make_unique<T[]>(n) }
{
}

template<typename T>
fixed_array<T>::fixed_array(std::size_t n, const value_type& def)
  : m_length{ n }
  , m_buffer{ std::make_unique<T[]>(m_length) }
{
    std::fill(begin(), end(), def);
}

template<typename T>
fixed_array<T>::fixed_array(const fixed_array& o)
  : m_length{ o.m_length }
  , m_buffer{ std::make_unique<T[]>(o.m_length) }
{
    std::copy(o.begin(), o.end(), begin());
}

template<typename T>
fixed_array<T>&
fixed_array<T>::operator=(const fixed_array& o)
{
    auto tmp = std::make_unique<T[]>(o.m_length);

    std::copy(o.begin(), o.end(), tmp.get());

    m_length = o.m_length;
    m_buffer = std::move(tmp);

    return *this;
}

template<class T>
fixed_array<T>::operator bool() const noexcept
{
    return m_buffer.get() != nullptr;
}

template<typename T>
std::size_t
fixed_array<T>::size() const noexcept
{
    return m_length;
}

template<typename T>
bool
fixed_array<T>::empty() const noexcept
{
    return m_length == 0;
}

template<typename T>
typename fixed_array<T>::iterator
fixed_array<T>::begin() noexcept
{
    return data();
}

template<typename T>
typename fixed_array<T>::const_iterator
fixed_array<T>::begin() const noexcept
{
    return data();
}

template<typename T>
typename fixed_array<T>::const_iterator
fixed_array<T>::cbegin() const noexcept
{
    return data();
}

template<typename T>
typename fixed_array<T>::iterator
fixed_array<T>::end() noexcept
{
    return data() + size();
}

template<typename T>
typename fixed_array<T>::const_iterator
fixed_array<T>::end() const noexcept
{
    return data() + size();
}

template<typename T>
typename fixed_array<T>::const_iterator
fixed_array<T>::cend() const noexcept
{
    return data() + size();
}

template<typename T>
typename fixed_array<T>::reverse_iterator
fixed_array<T>::rbegin() noexcept
{
    return { end() };
}

template<typename T>
typename fixed_array<T>::const_reverse_iterator
fixed_array<T>::rbegin() const noexcept
{
    return { end() };
}

template<typename T>
typename fixed_array<T>::const_reverse_iterator
fixed_array<T>::crbegin() const noexcept
{
    return { end() };
}

template<typename T>
typename fixed_array<T>::reverse_iterator
fixed_array<T>::rend() noexcept
{
    return { begin() };
}

template<typename T>
typename fixed_array<T>::const_reverse_iterator
fixed_array<T>::rend() const noexcept
{
    return { begin() };
}

template<typename T>
typename fixed_array<T>::const_reverse_iterator
fixed_array<T>::crend() const noexcept
{
    return { begin() };
}

template<typename T>
T&
fixed_array<T>::front() noexcept
{
    return *begin();
}

template<typename T>
const T&
fixed_array<T>::front() const noexcept
{
    return *begin();
}

template<typename T>
T&
fixed_array<T>::back() noexcept
{
    return *(begin() + size() - 1);
}

template<typename T>
const T&
fixed_array<T>::back() const noexcept
{
    return *(begin() + size() - 1);
}

template<typename T>
T*
fixed_array<T>::data() noexcept
{
    return m_buffer.get();
}

template<typename T>
const T*
fixed_array<T>::data() const noexcept
{
    return m_buffer.get();
}

template<typename T>
T& fixed_array<T>::operator[](std::size_t i) noexcept
{
    return data()[i];
}

template<typename T>
const T& fixed_array<T>::operator[](std::size_t i) const noexcept
{
    return data()[i];
}

template<typename T>
T&
fixed_array<T>::operator()(std::size_t i) noexcept
{
    return data()[i];
}

template<typename T>
const T&
fixed_array<T>::operator()(std::size_t i) const noexcept
{
    return data()[i];
}

template<class T>
void
fixed_array<T>::swap(fixed_array& other) noexcept
{
    std::swap(m_buffer, other.m_buffer);
    std::swap(m_length, other.m_length);
}

template<class T>
constexpr inline void
swap(fixed_array<T>& lhs, fixed_array<T>& rhs) noexcept
{
    lhs.swap(rhs);
}

} // namespace baryonyx

#endif
