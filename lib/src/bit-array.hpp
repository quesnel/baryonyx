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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_BIT_ARRAY_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_BIT_ARRAY_HPP

#include <algorithm>
#include <limits>
#include <memory>

#include <cstdint>

namespace baryonyx {

/**
 * @brief Simple and efficient structure to store bit array.
 *
 * @details @c bit_array stores a smart pointer to a dynamically allocated
 *     array of unsigned integer (with @c std::make_unique<std::uint32_t[]>).
 *     @c bit_array guarantees to delete the pointer on destruction of the @c
 *     bit_array.
 */
class bit_array_impl
{
public:
    using underlying_type = uintptr_t;

    constexpr static inline int bit_per_block =
      std::numeric_limits<underlying_type>::digits;

    bit_array_impl()
      : m_size(0)
      , m_block_size(0)
    {}

    bit_array_impl(int size)
      : m_size(size)
      , m_block_size((m_size / bit_per_block) + 1)
      , m_data(std::make_unique<underlying_type[]>(m_block_size))
    {}

    bit_array_impl(const bit_array_impl& other)
      : m_size(other.m_size)
      , m_block_size(other.m_block_size)
      , m_data(std::make_unique<underlying_type[]>(m_block_size))
    {
        std::copy_n(other.m_data.get(), m_block_size, m_data.get());
    }

    bit_array_impl(bit_array_impl&& other) noexcept
      : m_size(other.m_size)
      , m_block_size(other.m_block_size)
      , m_data(std::move(other.m_data))
    {}

    bit_array_impl& operator=(const bit_array_impl& other)
    {
        if (this != &other) {
            if (!(m_size == other.m_size &&
                  m_block_size == other.m_block_size)) {
                m_size = other.m_size;
                m_block_size = other.m_block_size;
                m_data = std::make_unique<underlying_type[]>(m_block_size);
            }
            std::copy_n(other.m_data.get(), m_block_size, m_data.get());
        }

        return *this;
    }

    bit_array_impl& operator=(bit_array_impl&& other)
    {
        if (this != &other) {
            m_size = other.m_size;
            m_block_size = other.m_block_size;
            m_data = std::move(other.m_data);
        }

        return *this;
    }

    void swap(bit_array_impl& other) noexcept
    {
        if (this != &other) {
            std::swap(m_size, other.m_size);
            std::swap(m_block_size, other.m_block_size);
            m_data.swap(other.m_data);
        }
    }

    ~bit_array_impl() noexcept = default;

    /**
     * @brief Affect 1 to the bits at @c index.
     *
     * @param index
     */
    void set(int index) noexcept
    {
        m_data[b_index(index)] |= UINTMAX_C(1) << b_offset(index);
    }

    /**
     * @brief Affect 0 to the bits at @c index.
     *
     * @param index
     */
    void unset(int index) noexcept
    {
        m_data[b_index(index)] &= ~(UINTMAX_C(1) << b_offset(index));
    }

    void invert(int index) noexcept
    {
        m_data[b_index(index)] ^= UINTMAX_C(1) << b_offset(index);
    }

    int get(int index) const noexcept
    {
        return !!(m_data[b_index(index)] & (UINTMAX_C(1) << b_offset(index)));
    }

    int operator[](int index) const noexcept
    {
        return !!(m_data[b_index(index)] & (UINTMAX_C(1) << b_offset(index)));
    }

    template<typename Iterator>
    void assign(Iterator first, Iterator last) noexcept
    {
        int i = 0;
        for (; first != last; ++first)
            set(i++, *first);
    }

    void ones() noexcept
    {
        std::fill_n(m_data.get(), m_block_size, UINTMAX_MAX);
    }

    void zeros() noexcept
    {
        std::fill_n(m_data.get(), m_block_size, UINTMAX_C(0));
    }

    /**
     * @brief Return the number of bit in the @e bit_array.
     */
    int size() const noexcept
    {
        return m_size;
    }

    /**
     * @brief Returns the number of @c underlying_type (std::uint64_t) used to
     * store bits in the @c bit_array.
     */
    int block_size() const noexcept
    {
        return m_block_size;
    }

    bool operator==(const bit_array_impl& other) const noexcept
    {
        if (this == &other)
            return true;

        if (size() != other.size())
            return false;

        for (int i = 0; i != m_block_size; ++i)
            if (other.m_data[i] != m_data[i])
                return false;

        return true;
    }

    bool operator!=(const bit_array_impl& other) const noexcept
    {
        if (size() != other.size())
            return true;

        for (int i = 0; i != m_block_size; ++i)
            if (other.m_data[i] != m_data[i])
                return true;

        return false;
    }

protected:
    int m_size;
    int m_block_size;
    std::unique_ptr<underlying_type[]> m_data;

    constexpr static int b_index(int b) noexcept
    {
        return b / bit_per_block;
    }

    constexpr static int b_offset(int b) noexcept
    {
        return b % bit_per_block;
    }
};

class bit_array : public bit_array_impl
{
public:
    bit_array() = default;

    bit_array(int size)
      : bit_array_impl(size)
    {}

    bit_array(const bit_array& other) = default;
    bit_array(bit_array&& other) = default;
    bit_array& operator=(const bit_array& other) = default;
    bit_array& operator=(bit_array&& other) = default;

    ~bit_array() noexcept = default;

    int get(int index) const noexcept
    {
        return !!(m_data[b_index(index)] & (UINTMAX_C(1) << b_offset(index)));
    }

    int operator[](int index) const noexcept
    {
        return !!(m_data[b_index(index)] & (UINTMAX_C(1) << b_offset(index)));
    }
};

template<int Value1 = 0, int Value2 = 1>
class value_bit_array : public bit_array_impl
{
public:
    constexpr static inline int return_value[] = { Value1, Value2 };

    value_bit_array() = default;

    value_bit_array(int size)
      : bit_array_impl(size)
    {}

    value_bit_array(const value_bit_array& other) = default;
    value_bit_array(value_bit_array&& other) = default;
    value_bit_array& operator=(const value_bit_array& other) = default;
    value_bit_array& operator=(value_bit_array&& other) = default;

    ~value_bit_array() noexcept = default;

    int get(int index) const noexcept
    {
        return return_value[!!(m_data[b_index(index)] &
                               (UINTMAX_C(1) << b_offset(index)))];
    }

    int operator[](int index) const noexcept
    {
        return return_value[!!(m_data[b_index(index)] &
                               (UINTMAX_C(1) << b_offset(index)))];
    }
};

} // namespace baryonyx

#endif
