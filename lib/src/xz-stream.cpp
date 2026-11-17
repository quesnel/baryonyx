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

#include "xz-stream.hpp"
#include "private.hpp"
#include "problem.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <streambuf>
#include <vector>

#ifdef BARYONYX_HAVE_LZMA
#include <lzma.h>
#endif

namespace baryonyx {

bool
is_xz_stream(std::istream& is)
{
    static constexpr std::array<unsigned char, 6> magic{
        0xFD, '7', 'z', 'X', 'Z', 0x00
    };

    std::array<char, 6> header{};
    const auto start = is.tellg();
    is.read(header.data(), static_cast<std::streamsize>(header.size()));
    const bool found =
      is.gcount() == static_cast<std::streamsize>(header.size()) &&
      std::ranges::equal(header, magic, [](char a, unsigned char b) {
          return static_cast<unsigned char>(a) == b;
      });

    is.clear();
    is.seekg(start);
    return found;
}

#ifdef BARYONYX_HAVE_LZMA

namespace {

/**
 * Read-only stream buffer that decodes an xz input on the fly, so that the
 * parser reads a compressed file without a decompressed copy in memory.
 */
class xz_streambuf : public std::streambuf
{
public:
    explicit xz_streambuf(std::istream& in)
      : m_in(in)
      , m_input(buffer_size)
      , m_output(buffer_size)
    {
        // Concatenated streams, as produced by `cat a.xz b.xz` or pxz.
        m_initialized = lzma_stream_decoder(&m_stream,
                                            UINT64_MAX,
                                            LZMA_CONCATENATED) == LZMA_OK;
        setg(m_output.data(), m_output.data(), m_output.data());
    }

    xz_streambuf(const xz_streambuf&) = delete;
    xz_streambuf& operator=(const xz_streambuf&) = delete;

    ~xz_streambuf() override
    {
        lzma_end(&m_stream);
    }

    bool initialized() const noexcept
    {
        return m_initialized;
    }

    bool failed() const noexcept
    {
        return m_failed;
    }

protected:
    int_type underflow() override
    {
        if (gptr() < egptr())
            return traits_type::to_int_type(*gptr());

        if (!m_initialized || m_finished || m_failed)
            return traits_type::eof();

        m_stream.next_out = reinterpret_cast<std::uint8_t*>(m_output.data());
        m_stream.avail_out = m_output.size();

        while (m_stream.avail_out == m_output.size()) {
            if (m_stream.avail_in == 0 && !m_input_exhausted) {
                m_in.read(m_input.data(),
                          static_cast<std::streamsize>(m_input.size()));
                m_stream.next_in =
                  reinterpret_cast<const std::uint8_t*>(m_input.data());
                m_stream.avail_in = static_cast<std::size_t>(m_in.gcount());
                m_input_exhausted = !m_in;
            }

            const auto ret = lzma_code(
              &m_stream, m_input_exhausted ? LZMA_FINISH : LZMA_RUN);

            if (ret == LZMA_STREAM_END) {
                m_finished = true;
                break;
            }

            if (ret != LZMA_OK) {
                m_failed = true;
                break;
            }
        }

        const auto produced = m_output.size() - m_stream.avail_out;
        setg(m_output.data(),
             m_output.data(),
             m_output.data() + static_cast<std::ptrdiff_t>(produced));

        return produced ? traits_type::to_int_type(*gptr())
                        : traits_type::eof();
    }

private:
    static constexpr std::size_t buffer_size = 1u << 16;

    std::istream& m_in;
    std::vector<char> m_input;
    std::vector<char> m_output;
    lzma_stream m_stream{};
    bool m_initialized = false;
    bool m_input_exhausted = false;
    bool m_finished = false;
    bool m_failed = false;
};

} // anonymous namespace

raw_problem
make_problem_from_xz(const context& ctx, std::istream& compressed)
{
    xz_streambuf buffer(compressed);
    if (!buffer.initialized()) {
        error(ctx, "xz decoder initialization failed\n");
        return raw_problem(file_format_error_tag::bad_compression);
    }

    std::istream is(&buffer);
    auto pb = make_problem(ctx, is);

    if (buffer.failed()) {
        error(ctx, "corrupted or truncated xz file\n");
        return raw_problem(file_format_error_tag::bad_compression);
    }

    return pb;
}

bool
has_xz_support() noexcept
{
    return true;
}

#else

raw_problem
make_problem_from_xz(const context& ctx, std::istream& /*compressed*/)
{
    error(ctx,
          "xz compressed file: baryonyx was built without liblzma, "
          "decompress the file or rebuild with liblzma\n");
    return raw_problem(file_format_error_tag::unsupported_compression);
}

bool
has_xz_support() noexcept
{
    return false;
}

#endif

} // namespace baryonyx
