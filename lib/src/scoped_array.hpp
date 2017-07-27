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

#ifndef ORG_VLEPROJECT_LP_SCOPED_ARRAY_HPP
#define ORG_VLEPROJECT_LP_SCOPED_ARRAY_HPP

#include <memory>

namespace lp {

/**
 * @brief Manage a single pointer to an array of T.
 *
 * @details @c scoped_array stores a smart pointer to a dynamically
 *      allocated array (with @c std::make_unique<T[]>().  @c scoped_array
 *      garantes to delete the pointer on destruction of the @c
 *      scoped_array
 *
 * @note The length is not stored into this containter.
 */
template<typename T>
class scoped_array
{
public:
    using value_type = T;

private:
    std::unique_ptr<T[]> m_buffer;

public:
    /**
     * Constructs a container with @c n elements.
     *
     * @param n Container size (i.e the number of elemens in the container at
     * construction).
     *
     * @exception @c std::bad_alloc.
     */
    explicit scoped_array(std::size_t n);

    /**
     * Constructs a container with @c n elements initialized to @c def.
     *
     * @param n Container size (i.e the number of elemens in the container at
     * construction).
     * @param def Default value to assign for each elements.
     *
     * @exception @c std::bad_alloc.
     */
    scoped_array(std::size_t n, const value_type& def);

    scoped_array(scoped_array&& other) = default;
    scoped_array& operator=(scoped_array&& other) = default;
    scoped_array(const scoped_array&) = delete;
    scoped_array& operator=(const scoped_array&) = delete;

    ~scoped_array() noexcept = default;

    explicit operator bool() const noexcept;

    T* data() noexcept;
    const T* data() const noexcept;

    T& operator[](std::ptrdiff_t i) noexcept;
    const T& operator[](std::ptrdiff_t i) const noexcept;

    T& operator()(std::ptrdiff_t i) noexcept;
    const T& operator()(std::ptrdiff_t i) const noexcept;

    void swap(scoped_array& other) noexcept;
};

template<class T>
scoped_array<T>::scoped_array(std::size_t n)
  : m_buffer{ std::make_unique<T[]>(n) }
{
}

template<class T>
scoped_array<T>::scoped_array(std::size_t n, const value_type& def)
  : m_buffer{ std::make_unique<T[]>(n) }
{
    std::fill(m_buffer.get(), m_buffer.get() + n, def);
}

template<class T>
T& scoped_array<T>::operator[](std::ptrdiff_t i) noexcept
{
    return m_buffer[i];
}

template<class T>
const T& scoped_array<T>::operator[](std::ptrdiff_t i) const noexcept
{
    return m_buffer[i];
}

template<class T>
T&
scoped_array<T>::operator()(std::ptrdiff_t i) noexcept
{
    return m_buffer[i];
}

template<class T>
const T&
scoped_array<T>::operator()(std::ptrdiff_t i) const noexcept
{
    return m_buffer[i];
}

template<class T>
T*
scoped_array<T>::data() noexcept
{
    return m_buffer.get();
}

template<class T>
const T*
scoped_array<T>::data() const noexcept
{
    return m_buffer.get();
}

template<class T>
scoped_array<T>::operator bool() const noexcept
{
    return m_buffer.get() != nullptr;
}

template<class T>
void
scoped_array<T>::swap(scoped_array& other) noexcept
{
    auto tmp = other.m_buffer;
    other.m_buffer = tmp;
    m_buffer = tmp;
}

template<class T>
constexpr inline void
swap(scoped_array<T>& lhs, scoped_array<T>& rhs) noexcept
{
    lhs.swap(rhs);
}

} // namespace lp

#endif
