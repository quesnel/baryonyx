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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_RANDOM_ENGINE_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_RANDOM_ENGINE_HPP

#include <array>
#include <cstdint>

namespace baryonyx {

/**
 * xoshiro256++ by D. Blackman and S. Vigna, seeded with splitmix64 as its
 * authors recommend.
 *
 * Unlike `std::default_random_engine`, whose algorithm depends on the
 * standard library, it draws the same sequence on every platform, and its
 * 32-byte state is cheap to create for each row of the buffered solver.
 */
class xoshiro256pp
{
public:
    using result_type = std::uint64_t;

    explicit xoshiro256pp(result_type seed_value = 1) noexcept
    {
        seed(seed_value);
    }

    void seed(result_type value) noexcept
    {
        for (auto& word : m_state) {
            value += 0x9e3779b97f4a7c15ULL;
            auto z = value;
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            word = z ^ (z >> 31);
        }
    }

    static constexpr result_type min() noexcept
    {
        return 0;
    }

    static constexpr result_type max() noexcept
    {
        return ~result_type{ 0 };
    }

    result_type operator()() noexcept
    {
        const auto result = rotl(m_state[0] + m_state[3], 23) + m_state[0];
        const auto t = m_state[1] << 17;

        m_state[2] ^= m_state[0];
        m_state[3] ^= m_state[1];
        m_state[1] ^= m_state[2];
        m_state[0] ^= m_state[3];
        m_state[2] ^= t;
        m_state[3] = rotl(m_state[3], 45);

        return result;
    }

private:
    static constexpr result_type rotl(result_type x, int k) noexcept
    {
        return (x << k) | (x >> (64 - k));
    }

    std::array<result_type, 4> m_state;
};

} // namespace baryonyx

#endif
