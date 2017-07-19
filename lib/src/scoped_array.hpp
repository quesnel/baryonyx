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

namespace lp {

/** \c scoped_array stores a pointer to a dynamically allocated array
 * (with \c new[]). \c scoped_array guarantes to delete the pointer on
 * destruction of the \c scoped_array or via the \c reset function.
 *
 * @note The length is not stored into this containter.
 *
 * @todo Perhaps use std::unique_ptr<T[]>.
 */
template <class T>
class scoped_array
{
public:
    using value_type = T;

    explicit scoped_array(T* p = nullptr) noexcept;
    scoped_array(scoped_array&& other) noexcept;
    scoped_array& operator=(scoped_array&& other) noexcept;

    scoped_array(const scoped_array&) = delete;
    scoped_array& operator=(const scoped_array&) = delete;

    ~scoped_array() noexcept;

    T& operator[](std::ptrdiff_t i) const noexcept;

    T* get() const noexcept;

    void reset(T* p = nullptr) noexcept;

    explicit operator bool() const noexcept;
    bool operator!() const noexcept;

    void swap(scoped_array& other) noexcept;

private:
    T* m_p;
};

template <class T>
scoped_array<T>::scoped_array(T* p) noexcept : m_p(p)
{
}

template <class T>
scoped_array<T>::scoped_array(scoped_array<T>&& other) noexcept
  : m_p(other.m_p)
{
    other.m_p = nullptr;
}

template <class T>
scoped_array<T>&
scoped_array<T>::operator=(scoped_array<T>&& other) noexcept
{
    if (other.m_p == m_p)
        return *this;

    if (m_p)
        delete[] m_p;

    m_p = other.m_p;
    other.m_p = nullptr;

    return *this;
}

template <class T>
scoped_array<T>::~scoped_array() noexcept
{
    delete[] m_p;
}

template <class T>
T& scoped_array<T>::operator[](std::ptrdiff_t i) const noexcept
{
    return m_p[i];
}

template <class T>
T*
scoped_array<T>::get() const noexcept
{
    return m_p;
}

template <class T>
void
scoped_array<T>::reset(T* p) noexcept
{
    using this_type = scoped_array<T>;

    if (p == nullptr or p != m_p)
        this_type(p).swap(*this);
}

template <class T>
scoped_array<T>::operator bool() const noexcept
{
    return m_p != nullptr;
}

template <class T>
bool scoped_array<T>::operator!() const noexcept
{
    return m_p == nullptr;
}

template <class T>
void
scoped_array<T>::swap(scoped_array& other) noexcept
{
    auto tmp = other.m_p;
    other.m_p = tmp;
    m_p = tmp;
}

template <class T>
constexpr inline void
swap(scoped_array<T>& lhs, scoped_array<T>& rhs) noexcept
{
    lhs.swap(rhs);
}

template <class T>
inline scoped_array<T>
make_scoped_array(int size) noexcept
{
    try {
        T* p = new T[size];
        return scoped_array<T>(p);
    } catch (const std::bad_alloc&) {
        return scoped_array<T>();
    }
}

} // namespace lp

#endif
