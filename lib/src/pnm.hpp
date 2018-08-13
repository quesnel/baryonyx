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

#include <array>
#include <fstream>
#include <memory>

#include <cassert>
#include <cstdint>

namespace baryonyx {

struct colormap
{
    float m_lower, m_upper;

    struct rgb
    {
        rgb(std::uint8_t red_, std::uint8_t green_, std::uint8_t blue_)
          : red(red_)
          , green(green_)
          , blue(blue_)
        {}

        std::uint8_t red;
        std::uint8_t green;
        std::uint8_t blue;
    };

    template<typename T>
    colormap(T lower, T upper)
      : m_lower(static_cast<float>(lower))
      , m_upper(static_cast<float>(upper))
    {
        assert(lower < upper);
    }

    template<typename T>
    rgb operator()(T v_) noexcept
    {
        float v = static_cast<float>(v_);

        if (v < m_lower)
            v = m_lower;

        if (v > m_upper)
            v = m_upper;

        float dv = m_upper - m_lower;
        float r = 1.f, g = 1.f, b = 1.f;

        if (v < (m_lower + 0.25f * dv)) {
            r = 0.f;
            g = 4.f * (v - m_lower) / dv;
        } else if (v < (m_lower + 0.5f * dv)) {
            r = 0.f;
            b = 1.f + 4.f * (m_lower + 0.25f * dv - v) / dv;
        } else if (v < (m_lower + 0.75f * dv)) {
            r = 4.f * (v - m_lower - 0.5f * dv) / dv;
            b = 0.f;
        } else {
            g = 1.f + 4.f * (m_lower + 0.75f * dv - v) / dv;
            b = 0.f;
        }

        return rgb(static_cast<std::uint8_t>(r * 255.f),
                   static_cast<std::uint8_t>(g * 255.f),
                   static_cast<std::uint8_t>(b * 255.f));
    }
};

class pnm_iterator;

class pnm_iterator_entry
{
private:
    std::uint8_t* m_colors{ nullptr };

public:
    pnm_iterator_entry() = default;

    pnm_iterator_entry(std::uint8_t* pointer)
      : m_colors(pointer)
    {}

    pnm_iterator_entry& operator=(const colormap::rgb& value)
    {
        m_colors[0] = value.red;
        m_colors[1] = value.green;
        m_colors[2] = value.blue;

        return *this;
    }

    pnm_iterator_entry(const pnm_iterator_entry& other) = default;
    ~pnm_iterator_entry() = default;

    void advance() noexcept
    {
        if (m_colors == nullptr)
            return;

        m_colors += 3;
    }

    std::uint8_t red() const noexcept
    {
        return m_colors[0];
    }

    std::uint8_t green() const noexcept
    {
        return m_colors[1];
    }

    std::uint8_t blue() const noexcept
    {
        return m_colors[2];
    }

    bool operator==(const pnm_iterator_entry& other) const noexcept
    {
        return m_colors == other.m_colors;
    }

    bool operator<(const pnm_iterator_entry& other) const noexcept
    {
        return m_colors < other.m_colors;
    }

    friend pnm_iterator;
};

class pnm_iterator
  : public std::iterator<std::forward_iterator_tag, pnm_iterator_entry>
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

    pnm_iterator(const pnm_iterator& other) = default;

    pnm_iterator operator++(int)
    {
        if (m_entry.m_colors == nullptr)
            return pnm_iterator();

        pnm_iterator copy(*this);
        copy.m_entry.advance();
        return copy;
    }

    pnm_iterator& operator++()
    {
        if (m_entry.m_colors == nullptr)
            return *this;

        m_entry.advance();

        return *this;
    }

    const pnm_iterator_entry* operator->() const
    {
        assert(m_entry.m_colors && "Dereference pointer");

        return &m_entry;
    }

    pnm_iterator_entry operator*() const
    {
        assert(m_entry.m_colors && "Dereference pointer");

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
    unsigned int m_width;
    unsigned int m_heigth;

public:
    pnm_array(unsigned int m, unsigned int n)
      : m_width(n)
      , m_heigth(m)
    {
        if (m_width * m_heigth > 0)
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
        auto header{ fmt::format("P6 {} {} {} ", m_width, m_heigth, 255) };

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

    pnm_vector(std::string filename, unsigned int m, unsigned int loop)
      : m_ofs(filename, std::ios::binary)
      , m_heigth(loop)
      , m_width(m)
    {
        if (!m_ofs.is_open())
            return;

        if (m_heigth * m_width > 0) {
            m_buffer = std::make_unique<std::uint8_t[]>(size());

            auto header{ fmt::format("P6 {} {} {} ", m_width, m_heigth, 255) };
            m_ofs << header;
        }
    }

    ~pnm_vector()
    {
        for (; m_heigth > 0; --m_heigth)
            m_ofs.write(reinterpret_cast<char*>(m_buffer.get()), size());
    }

    explicit operator bool() const noexcept
    {
        return m_buffer.get() != nullptr && m_ofs.is_open();
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
