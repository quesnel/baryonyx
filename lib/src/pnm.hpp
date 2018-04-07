/* Copyright (C) 2018 INRA
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

#ifndef ORG_VLEPROJECT_BARYONYX_PNM_HPP
#define ORG_VLEPROJECT_BARYONYX_PNM_HPP

#include <fstream>
#include <memory>

#include <cassert>
#include <cstdint>

namespace baryonyx {

template<typename T>
struct colormap_2
{
    std::array<std::uint8_t, 3> lower_color{ { 000, 000, 255 } };
    std::array<std::uint8_t, 3> middle_color{ { 000, 000, 000 } };
    std::array<std::uint8_t, 3> upper_color{ { 255, 000, 000 } };
    T m_lower, m_mid, m_upper;

    colormap_2(T lower, T middle, T upper)
      : m_lower(lower)
      , m_mid(middle)
      , m_upper(upper)
    {
        assert(lower < middle and middle < upper);
    }

    std::array<std::uint8_t, 3> operator()(T v) noexcept
    {
        if (v < m_lower)
            return lower_color;

        if (m_upper < v)
            return upper_color;

        std::array<std::uint8_t, 3> ret;

        if (v < m_mid)
            for (int i = 0; i != 3; ++i)
                ret[i] = static_cast<std::uint8_t>(
                  lower_color[i] + ((v - m_mid) / (m_upper - m_mid)) *
                                     (middle_color[i] - lower_color[i]));
        else
            for (int i = 0; i != 3; ++i)
                ret[i] = static_cast<std::uint8_t>(
                  middle_color[i] + ((v - m_lower) / (m_mid - m_lower)) *
                                      (upper_color[i] - middle_color[i]));

        return ret;
    }
};

template<typename T>
struct colormap
{
    std::array<std::uint8_t, 3> lower_color{ { 000, 0, 255 } };
    std::array<std::uint8_t, 3> upper_color{ { 255, 0, 000 } };
    T m_lower, m_upper;

    colormap(T lower, T upper)
      : m_lower(lower)
      , m_upper(upper)
    {
        assert(lower < upper);
    }

    std::array<std::uint8_t, 3> operator()(T v) noexcept
    {
        if (v < m_lower)
            return lower_color;

        if (m_upper < v)
            return upper_color;

        std::array<std::uint8_t, 3> ret;

        for (int i = 0; i != 3; ++i)
            ret[i] = static_cast<std::uint8_t>(
              lower_color[i] + ((v - m_lower) / (m_upper - m_lower)) *
                                 (upper_color[i] - lower_color[i]));

        return ret;
    }
};

class pnm_iterator;

class pnm_iterator_entry
{
private:
    std::uint8_t* m_red;
    std::uint8_t* m_green;
    std::uint8_t* m_blue;

public:
    pnm_iterator_entry()
      : m_red(nullptr)
      , m_green(nullptr)
      , m_blue(nullptr)
    {}

    pnm_iterator_entry(std::uint8_t* pointer)
      : m_red(pointer)
      , m_green(pointer + 1)
      , m_blue(pointer + 2)
    {}

    pnm_iterator_entry& operator=(const std::array<uint8_t, 3>& value)
    {
        *m_red = value[0];
        *m_green = value[1];
        *m_blue = value[2];

        return *this;
    }

    pnm_iterator_entry(const pnm_iterator_entry& other) = default;
    ~pnm_iterator_entry() = default;

    std::uint8_t red() const noexcept
    {
        return *m_red;
    }

    std::uint8_t green() const noexcept
    {
        return *m_green;
    }

    std::uint8_t blue() const noexcept
    {
        return *m_blue;
    }

    bool operator==(const pnm_iterator_entry& other) const noexcept
    {
        return m_red == other.m_red;
    }

    bool operator<(const pnm_iterator_entry& other) const noexcept
    {
        return m_red < other.m_red;
    }

    friend pnm_iterator;
};

class pnm_iterator
{
private:
    pnm_iterator_entry m_entry;

public:
    pnm_iterator()
      : m_entry()
    {}

    pnm_iterator(std::uint8_t* ptr)
      : m_entry(ptr)
    {}

    pnm_iterator(const pnm_iterator& other)
      : m_entry(other.m_entry)
    {}

    pnm_iterator operator++(int)
    {
        if (m_entry.m_red == nullptr)
            return pnm_iterator();

        return pnm_iterator(m_entry.m_red + 3);
    }

    pnm_iterator& operator++()
    {
        if (m_entry.m_red == nullptr)
            return *this;

        m_entry.m_red += 3;
        m_entry.m_green += 3;
        m_entry.m_blue += 3;

        return *this;
    }

    const pnm_iterator_entry* operator->() const
    {
        assert(m_entry.m_red != nullptr && "Dereference pointer");

        return &m_entry;
    }

    pnm_iterator_entry operator*() const
    {
        assert(m_entry.m_red != nullptr && "Dereference pointer");

        pnm_iterator_entry ret(m_entry);
        return ret;
    }

    friend void swap(pnm_iterator& lhs, pnm_iterator& rhs)
    {
        auto tmp = lhs.m_entry;
        lhs.m_entry = rhs.m_entry;
        rhs.m_entry = tmp;
    }

    friend bool operator==(const pnm_iterator& lhs, const pnm_iterator& rhs)
    {
        return lhs.m_entry == rhs.m_entry;
    }

    friend bool operator!=(const pnm_iterator& lhs, const pnm_iterator& rhs)
    {
        return !(lhs.m_entry == rhs.m_entry);
    }

    friend bool operator<(const pnm_iterator& lhs, const pnm_iterator& rhs)
    {
        return lhs.m_entry < rhs.m_entry;
    }
};

class pnm_array
{
public:
    using value_type = std::uint8_t;
    using pointer_type = std::uint8_t*;
    using iterator = pnm_iterator;
    using reference = value_type&;
    using const_reference = const value_type&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

private:
    std::unique_ptr<std::uint8_t[]> m_buffer;
    unsigned int m_heigth;
    unsigned int m_width;

public:
    pnm_array(unsigned int h, unsigned int w)
      : m_heigth(h)
      , m_width(w)
    {
        if (h * w > 0)
            m_buffer = std::make_unique<std::uint8_t[]>(size());
    }

    explicit operator bool() const noexcept
    {
        return m_buffer.get() != nullptr;
    }

    size_type size() const noexcept
    {
        return 3 * m_heigth * m_width;
    }

    void operator()(std::string filename) const noexcept
    {
        std::ofstream ofs{ filename, std::ios::binary };
        auto header{ fmt::format("P6 {} {} {} ", m_heigth, m_width, 255) };

        ofs << header;
        ofs.write(reinterpret_cast<char*>(m_buffer.get()), size());
    }

    iterator begin() noexcept
    {
        return m_buffer.get();
    }

    iterator end() noexcept
    {
        return m_buffer.get() + size();
    }

    iterator begin(size_type row) noexcept
    {
        assert(row < m_heigth);

        return m_buffer.get() + row * m_width * 3;
    }

    iterator end(size_type row) noexcept
    {
        assert(row < m_heigth);

        return m_buffer.get() + (row + 1) * m_width * 3;
    }

    pointer_type operator()(size_type row, size_type col) noexcept
    {
        assert(row < m_heigth);
        assert(col < m_width);

        return m_buffer.get() + (row * m_width + col) * 3;
    }

    void clear() noexcept
    {
        std::fill(m_buffer.get(),
                  m_buffer.get() + size(),
                  static_cast<std::uint8_t>(0));
    }
};

class pnm_vector
{
private:
    std::unique_ptr<std::uint8_t[]> m_buffer;
    std::ofstream m_ofs;
    unsigned int m_heigth;
    unsigned int m_width;

public:
    using value_type = std::uint8_t;
    using pointer_type = std::uint8_t*;
    using iterator = pnm_iterator;
    using reference = value_type&;
    using const_reference = const value_type&;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    pnm_vector(std::string filename, unsigned int h, unsigned int w)
      : m_ofs(filename, std::ios::binary)
      , m_heigth(h)
      , m_width(w)
    {
        if (not m_ofs.is_open())
            return;

        if (h * w > 0) {
            m_buffer = std::make_unique<std::uint8_t[]>(size());

            auto header{ fmt::format("P6 {} {} {} ", m_heigth, m_width, 255) };
            m_ofs << header;
        }
    }

    explicit operator bool() const noexcept
    {
        return m_buffer.get() != nullptr and m_ofs.is_open();
    }

    size_type size() const noexcept
    {
        return 3 * m_width;
    }

    iterator begin() noexcept
    {
        return m_buffer.get();
    }

    iterator end() noexcept
    {
        return m_buffer.get() + size();
    }

    void clear() noexcept
    {
        std::fill(m_buffer.get(),
                  m_buffer.get() + size(),
                  static_cast<std::uint8_t>(0));
    }

    void flush() noexcept
    {
        if (m_heigth == 0)
            return;

        m_ofs.write(reinterpret_cast<char*>(m_buffer.get()), size());
        --m_heigth;
    }
};

} // namespace baryonyx

#endif
