/* Copyright (C) 2016-2021 INRAE
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

#include "debug.hpp"

#include <algorithm>
#include <limits>
#include <memory>

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
    using underlying_type = std::size_t;

    static inline constexpr size_t k_one = underlying_type{ 1 };
    static inline constexpr size_t k_ones =
      std::numeric_limits<underlying_type>::max();
    static inline constexpr size_t k_zeros = underlying_type{ 0 };

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
    {
        std::fill_n(m_data.get(), m_block_size, k_zeros);
    }

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

    underlying_type block(int index) const noexcept
    {
        bx_assert(index >= 0 && index < m_block_size);

        return m_data[index];
    }

    void set_block(int index, underlying_type value) noexcept
    {
        m_data[index] = value;
    }

    /**
     * @brief Affect 1 to the bits at @c index.
     *
     * @param index
     */
    void set(int index) noexcept
    {
        bx_assert(index >= 0 && index < m_size);

        m_data[b_index(index)] |= k_one << b_offset(index);
    }

    void set(int index, bool val) noexcept
    {
        auto mask = (k_one << b_offset(index));
        auto& ref = m_data[b_index(index)];
        ref ^= (underlying_type(val) ^ ref) & mask;
    }

    /**
     * @brief Affect 0 to the bits at @c index.
     *
     * @param index
     */
    void unset(int index) noexcept
    {
        bx_assert(index >= 0 && index < m_size);

        m_data[b_index(index)] &= ~(k_one << b_offset(index));
    }

    void invert(int index) noexcept
    {
        bx_assert(index >= 0 && index < m_size);

        m_data[b_index(index)] ^= k_one << b_offset(index);
    }

    int get(int index) const noexcept
    {
        bx_assert(index >= 0 && index < m_size);

        return !!(m_data[b_index(index)] & (k_one << b_offset(index)));
    }

    int operator[](int index) const noexcept
    {
        bx_assert(index >= 0 && index < m_size);

        return !!(m_data[b_index(index)] & (k_one << b_offset(index)));
    }

    void assign(const bit_array_impl& other, int from, int to)
    {
        bx_assert(other.size() == size());
        bx_assert(from >= 0 && from < to && from < size() && to <= size());

        int loop_min, loop_max;

        {
            const auto index = b_index(from);
            const auto offset = b_offset(from);
            loop_min = index;

            if (offset > 0) {
                ++loop_min;

                underlying_type value_orig = m_data[index];
                value_orig = (value_orig >> (bit_per_block - offset));
                value_orig = (value_orig << (bit_per_block - offset));

                underlying_type value_new = other.m_data[index];
                value_new = (value_new << offset);
                value_new = (value_new >> offset);

                m_data[index] = value_new | value_orig;
            }
        }

        {
            const auto index = b_index(to);
            const auto offset = b_offset(to);
            loop_max = index;

            if (offset > 0) {
                --loop_max;

                underlying_type value_orig = m_data[index];
                value_orig = (value_orig << offset);
                value_orig = (value_orig >> offset);

                underlying_type value_new = other.m_data[index];
                value_new = (value_new >> (bit_per_block - offset));
                value_new = (value_new << (bit_per_block - offset));

                m_data[index] = value_new | value_orig;
            }
        }

        if (loop_max > loop_min)
            std::copy_n(other.m_data.get() + loop_min,
                        loop_max - loop_min,
                        m_data.get() + loop_min);
    }

    void ones() noexcept
    {
        std::fill_n(m_data.get(), m_block_size, k_ones);
    }

    void zeros() noexcept
    {

        std::fill_n(m_data.get(), m_block_size, k_zeros);
    }

    bool empty() const noexcept
    {
        return m_data.get() == nullptr;
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
        return !!(m_data[b_index(index)] & (k_one << b_offset(index)));
    }

    int operator[](int index) const noexcept
    {
        return !!(m_data[b_index(index)] & (k_one << b_offset(index)));
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
                               (k_one << b_offset(index)))];
    }

    int operator[](int index) const noexcept
    {
        return return_value[!!(m_data[b_index(index)] &
                               (k_one << b_offset(index)))];
    }
};

//
// Swap functions
//

inline void
swap(bit_array& lhs, bit_array& rhs) noexcept
{
    lhs.swap(rhs);
}

template<int Value1, int Value2>
void
swap(value_bit_array<Value1, Value2>& lhs,
     value_bit_array<Value1, Value2>& rhs) noexcept
{
    lhs.swap(rhs);
}

//
// Hash functions
//

template<int Value1, int Value2>
struct value_bit_array_hash
{
    std::size_t operator()(const value_bit_array<Value1, Value2>& array) const
      noexcept
    {
        std::size_t ret{ 0 };

        for (int i = 0, e = array.block_size(); i != e; ++i)
            ret ^=
              std::hash<bit_array_impl::underlying_type>()(array.block(i)) +
              0x9e3779b9 + (ret << 6) + (ret >> 2);

        return ret;
    }
};

struct bit_array_hash
{
    std::size_t operator()(const bit_array& array) const noexcept
    {
        std::size_t ret{ 0 };

        for (int i = 0, e = array.block_size(); i != e; ++i)
            ret ^=
              std::hash<bit_array_impl::underlying_type>()(array.block(i)) +
              0x9e3779b9 + (ret << 6) + (ret >> 2);

        return ret;
    }
};

} // namespace baryonyx

#endif
