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

#ifndef ORG_VLEPROJECT_BARYONYX_SOLVER_RANDOM_DISTRIBUTIONS_HPP
#define ORG_VLEPROJECT_BARYONYX_SOLVER_RANDOM_DISTRIBUTIONS_HPP

#include <cmath>
#include <concepts>
#include <cstdint>
#include <iterator>
#include <limits>
#include <utility>

/**
 * Replacements for the standard distributions, whose algorithms depend on
 * the standard library: with these, a seed gives the same draws on every
 * platform. They require an engine producing uniform 64-bit words.
 */

namespace baryonyx {

template<typename Engine>
concept uniform_64_engine = requires(Engine& rng) {
    { rng() } -> std::same_as<std::uint64_t>;
} && Engine::min() == 0 && Engine::max() == ~std::uint64_t{ 0 };

/** Uniform integer in `[0, bound)`, `bound > 0`, without modulo bias. */
template<uniform_64_engine Engine>
std::uint64_t
uniform_below(Engine& rng, const std::uint64_t bound) noexcept
{
    // Draws below `threshold` would make the low residues more likely.
    const std::uint64_t threshold = (0 - bound) % bound;

    for (;;)
        if (const auto r = rng(); r >= threshold)
            return r % bound;
}

/** Uniform double in `[0, 1)` built from the 53 high bits of a draw. */
template<uniform_64_engine Engine>
double
uniform_unit(Engine& rng) noexcept
{
    return static_cast<double>(rng() >> 11) * 0x1.0p-53;
}

/** Fisher-Yates shuffle; call it qualified to avoid `std::shuffle`. */
template<std::random_access_iterator Iterator, uniform_64_engine Engine>
void
shuffle(Iterator first, Iterator last, Engine& rng) noexcept
{
    const auto size = static_cast<std::uint64_t>(last - first);

    for (std::uint64_t i = size; i > 1; --i) {
        const auto j = uniform_below(rng, i);
        using std::swap;
        swap(first[static_cast<std::ptrdiff_t>(i - 1)],
             first[static_cast<std::ptrdiff_t>(j)]);
    }
}

class bernoulli_distribution
{
public:
    explicit bernoulli_distribution(double p = 0.5) noexcept
      : m_p(p)
    {}

    template<uniform_64_engine Engine>
    bool operator()(Engine& rng) const noexcept
    {
        return uniform_unit(rng) < m_p;
    }

private:
    double m_p;
};

/** Uniform integer in `[a, b]`. */
template<std::integral T>
class uniform_int_distribution
{
public:
    explicit uniform_int_distribution(
      T a = 0,
      T b = std::numeric_limits<T>::max()) noexcept
      : m_a(a)
      , m_b(b)
    {}

    template<uniform_64_engine Engine>
    T operator()(Engine& rng) const noexcept
    {
        const auto range = static_cast<std::uint64_t>(m_b) -
                           static_cast<std::uint64_t>(m_a) + 1u;

        // A zero range means the whole 64-bit span.
        const auto offset = range == 0 ? rng() : uniform_below(rng, range);

        return static_cast<T>(static_cast<std::uint64_t>(m_a) + offset);
    }

private:
    T m_a;
    T m_b;
};

/** Uniform real in `[a, b)`. */
template<std::floating_point T>
class uniform_real_distribution
{
public:
    uniform_real_distribution(T a, T b) noexcept
      : m_a(a)
      , m_b(b)
    {}

    template<uniform_64_engine Engine>
    T operator()(Engine& rng) const noexcept
    {
        return static_cast<T>(m_a + (m_b - m_a) * uniform_unit(rng));
    }

private:
    T m_a;
    T m_b;
};

/** Normal law with Marsaglia's polar method, which yields values in pairs. */
template<std::floating_point T>
class normal_distribution
{
public:
    normal_distribution(T mean, T stddev) noexcept
      : m_mean(mean)
      , m_stddev(stddev)
    {}

    T mean() const noexcept
    {
        return m_mean;
    }

    T stddev() const noexcept
    {
        return m_stddev;
    }

    template<uniform_64_engine Engine>
    T operator()(Engine& rng) noexcept
    {
        if (m_has_spare) {
            m_has_spare = false;
            return static_cast<T>(m_mean + m_stddev * m_spare);
        }

        double u, v, s;
        do {
            u = 2.0 * uniform_unit(rng) - 1.0;
            v = 2.0 * uniform_unit(rng) - 1.0;
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);

        const double factor = std::sqrt(-2.0 * std::log(s) / s);
        m_spare = v * factor;
        m_has_spare = true;

        return static_cast<T>(m_mean + m_stddev * u * factor);
    }

private:
    double m_spare = 0.0;
    T m_mean;
    T m_stddev;
    bool m_has_spare = false;
};

} // namespace baryonyx

#endif
