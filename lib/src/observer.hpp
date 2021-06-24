/* Copyright (C) 2018-2021 INRAE
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

#ifndef ORG_VLEPROJECT_BARYONYX_PNM_OBSERVER_HPP
#define ORG_VLEPROJECT_BARYONYX_PNM_OBSERVER_HPP

#include "pnm.hpp"
#include "sparse-matrix.hpp"

#include <algorithm>
#include <fstream>
#include <memory>
#include <utility>

namespace baryonyx {

namespace details {

class pi_pnm_observer
{
private:
    const int constraints;
    pnm_vector m_pnm;

public:
    pi_pnm_observer(const std::string_view filename,
                    int m,
                    int /*n*/,
                    long int loop)
      : constraints(m)
      , m_pnm(fmt::format("{}-pi.pnm", filename), constraints, loop)
    {}

    template<typename Float>
    void make_observation(const sparse_matrix<int>& /*ap*/,
                          const Float* /*P*/,
                          const Float* pi)
    {
        static_assert(std::is_floating_point<Float>::value);

        auto it = m_pnm.begin();
        colormap fct(-1.f, +1.f);

        for (int i = 0; i != constraints; ++i)
            *it = fct(pi[i]);

        m_pnm.flush();
    }
};

class ap_pnm_observer
{
private:
    const int constraints;
    const int variables;
    int m_frame;
    std::string m_basename;

public:
    ap_pnm_observer(const std::string_view filename,
                    int m,
                    int n,
                    long int /*loop*/)
      : constraints(m)
      , variables(n)
      , m_frame(0)
      , m_basename(std::move(filename))
    {}

    template<typename Float>
    void make_observation(const sparse_matrix<int>& ap,
                          const Float* P,
                          const Float* /*pi*/)
    {
        static_assert(std::is_floating_point<Float>::value);

        colormap cm(-1.0f, 1.0f);
        pnm_array pnm(constraints, variables);
        if (!pnm)
            return;

        pnm.clear();
        for (int k = 0; k != constraints; ++k) {
            auto [it, et] = ap.row(k);

            for (; it != et; ++it) {
                std::uint8_t* pointer = pnm(k, it->column);
                auto color_rgb = cm(P[it->value]);

                pointer[0] = color_rgb.red;
                pointer[1] = color_rgb.green;
                pointer[2] = color_rgb.blue;
            }
        }

        pnm(fmt::format("{}-P-{}.pnm", m_basename, m_frame++));
    }
};

class pi_file_observer
{
private:
    const int constraints;
    int m_len;
    std::ofstream m_ofs;

public:
    pi_file_observer(const std::string_view filename,
                     int m,
                     int /*n*/,
                     long int /*loop*/)
      : constraints(m)
      , m_ofs(fmt::format("{}-pi.txt", filename))
    {}

    template<typename Float>
    void make_observation(const sparse_matrix<int>& /*ap*/,
                          const Float* /*P*/,
                          const Float* pi)
    {
        static_assert(std::is_floating_point<Float>::value);

        if (m_ofs) {
            if (m_len > 1) {
                m_ofs << pi[0];

                for (int i = 1; i != constraints; ++i)
                    m_ofs << ' ' << pi[i];

            } else {
                m_ofs << pi[0];
            }

            m_ofs << '\n';
        }
    }
};

class ap_file_observer
{
private:
    const int constraints;
    const int variables;
    int m_frame;
    std::string m_basename;
    std::vector<float> m_value;

public:
    ap_file_observer(const std::string_view filename,
                     int m,
                     int n,
                     long int /*loop*/)
      : constraints(m)
      , variables(n)
      , m_frame(0)
      , m_basename(std::move(filename))
      , m_value(variables)
    {}

    template<typename Float>
    void make_observation(const sparse_matrix<int>& ap,
                          const Float* P,
                          const Float* /*pi*/)
    {
        static_assert(std::is_floating_point<Float>::value);

        std::ofstream ofs(fmt::format("{}-P-{}.txt", m_basename, m_frame++));

        for (int k = 0; k != constraints; ++k) {
            std::fill_n(std::begin(m_value), variables, 0.0f);

            auto [it, et] = ap.row(k);
            for (; it != et; ++it)
                m_value[it->column] = static_cast<float>(P[it->value]);

            std::copy(m_value.begin(),
                      m_value.end(),
                      std::ostream_iterator<Float>(ofs, " "));

            ofs << '\n';
        }
    }
};
}

class pnm_observer
{
    details::pi_pnm_observer m_pi_obs;
    details::ap_pnm_observer m_ap_obs;

public:
    pnm_observer(const std::string_view filename, int m, int n, long int loop)
      : m_pi_obs(filename, m, n, loop)
      , m_ap_obs(filename, m, n, loop)
    {}

    template<typename Float>
    void make_observation(const sparse_matrix<int>& ap,
                          const Float* P,
                          const Float* pi)
    {
        static_assert(std::is_floating_point<Float>::value);

        m_pi_obs.make_observation(ap, P, pi);
        m_ap_obs.make_observation(ap, P, pi);
    }
};

class file_observer
{
    details::pi_file_observer m_pi_obs;
    details::ap_file_observer m_ap_obs;

public:
    file_observer(const std::string_view filename, int m, int n, long int loop)
      : m_pi_obs(filename, m, n, loop)
      , m_ap_obs(filename, m, n, loop)
    {}

    template<typename Float>
    void make_observation(const sparse_matrix<int>& ap,
                          const Float* P,
                          const Float* pi)
    {
        static_assert(std::is_floating_point<Float>::value);

        m_pi_obs.make_observation(ap, P, pi);
        m_ap_obs.make_observation(ap, P, pi);
    }
};

class none_observer
{
public:
    none_observer(const std::string_view /*filename*/,
                  int /*m*/,
                  int /*n*/,
                  long int /*loop*/)
    {}

    template<typename Float>
    void make_observation(const sparse_matrix<int>& /*ap*/,
                          const Float* /*P*/,
                          const Float* /*pi*/)
    {
        static_assert(std::is_floating_point<Float>::value);
    }
};

} // namespace baryonyx

#endif
