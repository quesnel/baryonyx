/* Copyright (C) 2016-2026 INRAE
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
#include <cstdint>
#include <limits>
#include <memory>
#include <type_traits>
#include <utility>

namespace baryonyx {

/**
 * Fixed-size array of bits stored in 64-bit blocks. The bits past `size()`
 * in the last block are always zero so that `operator==` and the hash can
 * compare whole blocks.
 */
class bit_array_impl
{
public:
    using underlying_type = std::uint64_t;

    static constexpr underlying_type k_one = 1;
    static constexpr underlying_type k_ones =
      std::numeric_limits<underlying_type>::max();
    static constexpr underlying_type k_zeros = 0;

    static constexpr int bit_per_block =
      std::numeric_limits<underlying_type>::digits;

    bit_array_impl()
      : m_size(0)
      , m_block_size(0)
    {
    }

    explicit bit_array_impl(int size)
      : m_size(size)
      , m_block_size((size + bit_per_block - 1) / bit_per_block)
      , m_data(std::make_unique<underlying_type[]>(m_block_size))
    {
        bx_assert(size >= 0);
    }

    bit_array_impl(const bit_array_impl& other)
      : m_size(other.m_size)
      , m_block_size(other.m_block_size)
      , m_data(std::make_unique<underlying_type[]>(m_block_size))
    {
        std::copy_n(other.m_data.get(), m_block_size, m_data.get());
    }

    bit_array_impl(bit_array_impl&& other) noexcept
      : m_size(std::exchange(other.m_size, 0))
      , m_block_size(std::exchange(other.m_block_size, 0))
      , m_data(std::move(other.m_data))
    {
    }

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

    bit_array_impl& operator=(bit_array_impl&& other) noexcept
    {
        if (this != &other) {
            m_size = std::exchange(other.m_size, 0);
            m_block_size = std::exchange(other.m_block_size, 0);
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
        bx_assert(index >= 0 && index < m_block_size);

        m_data[index] = index + 1 == m_block_size ? value & last_block_mask()
                                                  : value;
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
        bx_assert(index >= 0 && index < m_size);

        // Negating val spreads it over the whole word before masking.
        const auto mask = k_one << b_offset(index);
        auto& ref = m_data[b_index(index)];
        ref ^= (-static_cast<underlying_type>(val) ^ ref) & mask;
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

    void ones() noexcept
    {
        std::fill_n(m_data.get(), m_block_size, k_ones);

        if (m_block_size > 0)
            m_data[m_block_size - 1] &= last_block_mask();
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

        return std::equal(
          m_data.get(), m_data.get() + m_block_size, other.m_data.get());
    }

protected:
    int m_size;
    int m_block_size;
    std::unique_ptr<underlying_type[]> m_data;

    // Indices are non-negative: unsigned arithmetic turns the division and
    // the modulo into a shift and a mask.
    constexpr static unsigned b_index(int b) noexcept
    {
        return static_cast<unsigned>(b) / bit_per_block;
    }

    constexpr static unsigned b_offset(int b) noexcept
    {
        return static_cast<unsigned>(b) % bit_per_block;
    }

    underlying_type last_block_mask() const noexcept
    {
        const auto used = b_offset(m_size);
        return used == 0 ? k_ones : (k_one << used) - 1;
    }
};

class bit_array : public bit_array_impl
{
public:
    bit_array() = default;

    explicit bit_array(int size)
      : bit_array_impl(size)
    {
    }

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
    static constexpr int return_value[] = { Value1, Value2 };

    value_bit_array() = default;

    explicit value_bit_array(int size)
      : bit_array_impl(size)
    {
    }

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

inline std::size_t
hash_bit_array(const bit_array_impl& array) noexcept
{
    static_assert(
      std::is_same_v<bit_array_impl::underlying_type, std::uint64_t>);

    // Mixing steps of the MurmurHash3 64 bits finalizer.
    std::size_t seed = static_cast<std::size_t>(array.block_size());

    for (int i = 0, e = array.block_size(); i != e; ++i) {
        seed ^= array.block(i);
        seed ^= (seed >> 33);
        seed *= 0xff51afd7ed558ccd;
        seed ^= (seed >> 33);
        seed *= 0xc4ceb9fe1a85ec53;
        seed ^= (seed >> 33);
    }

    return seed;
}

template<int Value1, int Value2>
struct value_bit_array_hash
{
    std::size_t operator()(
      const value_bit_array<Value1, Value2>& array) const noexcept
    {
        return hash_bit_array(array);
    }
};

struct bit_array_hash
{
    std::size_t operator()(const bit_array& array) const noexcept
    {
        return hash_bit_array(array);
    }
};

} // namespace baryonyx

#endif
