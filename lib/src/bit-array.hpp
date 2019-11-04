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

#include <cstddef>
#include <cstdint>

#include <algorithm>
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
class bit_array
{
public:
    // Define the number of bits available in a block
    constexpr static inline int block_t =
      static_cast<int>(sizeof(std::uint32_t) * 8);

    bit_array()
      : m_size(0)
      , m_block_size(0)
    {}

    bit_array(int size)
      : m_size(size)
      , m_block_size((m_size / block_t) + 1)
      , m_data(std::make_unique<std::uint32_t[]>(m_block_size))
    {}

    bit_array(const bit_array& other)
      : m_size(other.m_size)
      , m_block_size(other.m_block_size)
      , m_data(std::make_unique<std::uint32_t[]>(m_block_size))
    {
        std::copy_n(other.m_data.get(), m_block_size, m_data.get());
    }

    bit_array(bit_array&& other) noexcept
      : m_size(other.m_size)
      , m_block_size(other.m_block_size)
      , m_data(std::move(other.m_data))
    {}

    bit_array& operator=(const bit_array& other)
    {
        if (this != &other) {
            if (!(m_size == other.m_size &&
                  m_block_size == other.m_block_size)) {
                m_size = other.m_size;
                m_block_size = other.m_block_size;
                m_data = std::make_unique<std::uint32_t[]>(m_block_size);
            }
            std::copy_n(other.m_data.get(), m_block_size, m_data.get());
        }

        return *this;
    }

    bit_array& operator=(bit_array&& other)
    {
        if (this != &other) {
            m_size = other.m_size;
            m_block_size = other.m_block_size;
            m_data = std::move(other.m_data);
        }

        return *this;
    }

    void swap(bit_array& other) noexcept
    {
        if (this != &other) {
            std::swap(m_size, other.m_size);
            std::swap(m_block_size, other.m_block_size);
            m_data.swap(other.m_data);
        }
    }

    ~bit_array() noexcept = default;

    /**
     * @brief Affect 1 to the bits at @c index.
     *
     * @param index
     */
    void set(int index) noexcept
    {
        m_data[b_index(index)] |= 1u << b_offset(index);
    }

    /**
     * @brief Affect 0 to the bits at @c index.
     *
     * @param index
     */
    void unset(int index) noexcept
    {
        m_data[b_index(index)] &= ~(1u << b_offset(index));
    }

    void invert(int index) noexcept
    {
        m_data[b_index(index)] ^= 1u << b_offset(index);
    }

    int get(int index) const noexcept
    {
        return !!(m_data[b_index(index)] & (1u << b_offset(index)));
    }

    int operator[](int index) const noexcept
    {
        return !!(m_data[b_index(index)] & (1u << b_offset(index)));
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
        std::fill_n(m_data.get(), m_block_size, 0xffffffff);
    }

    void zeros() noexcept
    {
        std::fill_n(m_data.get(), m_block_size, 0u);
    }

    int size() const noexcept
    {
        return m_size;
    }

    int block_size() const noexcept
    {
        return m_block_size;
    }

private:
    int m_size;
    int m_block_size;
    std::unique_ptr<std::uint32_t[]> m_data;

    constexpr static int b_index(int b) noexcept
    {
        return b / block_t;
    }

    constexpr static int b_offset(int b) noexcept
    {
        return b % block_t;
    }
};

template<int Value1 = 0, int Value2 = 1>
class value_bit_array
{
public:
    // Define the number of bits available in a block
    constexpr static inline int block_t =
      static_cast<int>(sizeof(std::uint32_t) * 8);

    constexpr static inline int return_value[] = { Value1, Value2 };

    value_bit_array(int size)
      : m_size(size)
      , m_block_size((m_size / block_t) + 1)
      , m_data(std::make_unique<std::uint32_t[]>(m_block_size))
    {}

    value_bit_array(const value_bit_array& other)
      : m_size(other.m_size)
      , m_block_size(other.m_block_size)
      , m_data(std::make_unique<std::uint32_t[]>(m_block_size))
    {
        std::copy_n(other.m_data.get(), m_block_size, m_data.get());
    }

    value_bit_array(value_bit_array&& other) noexcept
      : m_size(other.m_size)
      , m_block_size(other.m_block_size)
      , m_data(std::move(other.m_data))
    {}

    value_bit_array& operator=(const value_bit_array& other)
    {
        if (this != &other) {
            if (!(m_size == other.m_size &&
                  m_block_size == other.m_block_size)) {
                m_size = other.m_size;
                m_block_size = other.m_block_size;
                m_data = std::make_unique<std::uint32_t[]>(m_block_size);
            }
            std::copy_n(other.m_data.get(), m_block_size, m_data.get());
        }

        return *this;
    }

    value_bit_array& operator=(value_bit_array&& other)
    {
        if (this != &other) {
            m_size = other.m_size;
            m_block_size = other.m_block_size;
            m_data = std::move(other.m_data);
        }

        return *this;
    }

    void swap(value_bit_array& other) noexcept
    {
        if (this != &other) {
            std::swap(m_size, other.m_size);
            std::swap(m_block_size, other.m_block_size);
            m_data.swap(other.m_data);
        }
    }

    ~value_bit_array() noexcept = default;

    /**
     * @brief Affect 1 to the bits at @c index.
     *
     * @param index
     */
    void set(int index) noexcept
    {
        m_data[b_index(index)] |= 1u << b_offset(index);
    }

    /**
     * @brief Affect 0 to the bits at @c index.
     *
     * @param index
     */
    void unset(int index) noexcept
    {
        m_data[b_index(index)] &= ~(1u << b_offset(index));
    }

    void invert(int index) noexcept
    {
        m_data[b_index(index)] ^= 1u << b_offset(index);
    }

    int get(int index) const noexcept
    {
        return return_value[!!(m_data[b_index(index)] &
                               (1u << b_offset(index)))];
    }

    int operator[](int index) const noexcept
    {
        return return_value[!!(m_data[b_index(index)] &
                               (1u << b_offset(index)))];
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
        std::fill_n(m_data.get(), m_block_size, 0xffffffff);
    }

    void zeros() noexcept
    {
        std::fill_n(m_data.get(), m_block_size, 0u);
    }

    int size() const noexcept
    {
        return m_size;
    }

    int block_size() const noexcept
    {
        return m_block_size;
    }

private:
    int m_size;
    int m_block_size;
    std::unique_ptr<std::uint32_t[]> m_data;

    constexpr static int b_index(int b) noexcept
    {
        return b / block_t;
    }

    constexpr static int b_offset(int b) noexcept
    {
        return b % block_t;
    }
};

} // namespace baryonyx

#endif
