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

template<typename floatingpointT>
class pi_pnm_observer
{
private:
    const std::unique_ptr<floatingpointT[]>& m_pi;
    int m_m;
    int m_loop;
    pnm_vector m_pnm;

public:
    pi_pnm_observer(std::string filename,
                    const std::unique_ptr<floatingpointT[]>& pi,
                    int m,
                    int loop)
      : m_pi(pi)
      , m_m(m)
      , m_loop(loop)
      , m_pnm(fmt::format("{}-pi.pnm", filename), m, loop)
    {}

    void make_observation()
    {
        std::transform(
          m_pi.get(), m_pi.get() + m_m, m_pnm.begin(), colormap(-1.0f, +1.0f));

        m_pnm.flush();
    }
};

template<typename floatingpointT>
class ap_pnm_observer
{
private:
    std::string m_basename;
    const sparse_matrix<int> m_ap;
    const std::unique_ptr<floatingpointT[]>& m_P;
    int m_m;
    int m_n;
    int m_frame;

public:
    ap_pnm_observer(std::string basename,
                    sparse_matrix<int> ap,
                    const std::unique_ptr<floatingpointT[]>& P,
                    int m,
                    int n)
      : m_basename(std::move(basename))
      , m_ap(std::move(ap))
      , m_P(P)
      , m_m(m)
      , m_n(n)
      , m_frame(0)
    {}

    void make_observation()
    {
        colormap cm(-1.0f, 1.0f);
        pnm_array pnm(m_m, m_n);
        if (!pnm)
            return;

        // floatingpointT min_val = std::numeric_limits<floatingpointT>::max();
        // floatingpointT max_val = std::numeric_limits<floatingpointT>::min();

        pnm.clear();
        for (int k = 0; k != m_m; ++k) {
            sparse_matrix<int>::const_row_iterator it, et;
            std::tie(it, et) = m_ap.row(k);

            for (; it != et; ++it) {
                std::uint8_t* pointer = pnm(k, it->column);

                // min_val = std::min(min_val, m_P[it->value]);
                // max_val = std::max(max_val, m_P[it->value]);

                auto color_rgb = cm(m_P[it->value]);

                pointer[0] = color_rgb.red;
                pointer[1] = color_rgb.green;
                pointer[2] = color_rgb.blue;
            }
        }

        // fmt::print("range [{} {}]\n", min_val, max_val);

        pnm(fmt::format("{}-P-{}.pnm", m_basename, m_frame++));
    }
};

template<typename floatingpointT>
class pi_file_observer
{
private:
    const std::unique_ptr<floatingpointT[]>& m_pi;
    int m_len;
    int m_loop;
    std::ofstream m_ofs;

public:
    pi_file_observer(std::string filename,
                     const std::unique_ptr<floatingpointT[]>& pi,
                     int len,
                     int loop)
      : m_pi(pi)
      , m_len(len)
      , m_loop(loop)
      , m_ofs(fmt::format("{}-pi.txt", filename))
    {}

    void make_observation()
    {
        if (m_ofs) {
            std::copy(m_pi.get(),
                      m_pi.get() + m_len,
                      std::ostream_iterator<floatingpointT>(m_ofs, " "));

            m_ofs << '\n';
        }
    }
};

template<typename floatingpointT>
class ap_file_observer
{
private:
    std::string m_basename;
    const sparse_matrix<int> m_ap;
    const std::unique_ptr<floatingpointT[]>& m_P;
    int m_m;
    int m_n;
    int m_frame;

public:
    ap_file_observer(std::string basename,
                     sparse_matrix<int> ap,
                     const std::unique_ptr<floatingpointT[]>& P,
                     int m,
                     int n)
      : m_basename(std::move(basename))
      , m_ap(std::move(ap))
      , m_P(P)
      , m_m(m)
      , m_n(n)
      , m_frame(0)
    {}

    void make_observation()
    {
        std::vector<floatingpointT> val(m_n);
        std::ofstream ofs(fmt::format("{}-P-{}.txt", m_basename, m_frame++));

        for (int k = 0; k != m_m; ++k) {
            val.resize(m_n, static_cast<floatingpointT>(0));

            sparse_matrix<int>::const_row_iterator it, et;
            std::tie(it, et) = m_ap.row(k);

            for (; it != et; ++it)
                val[it->column] = m_P[it->value];

            std::copy(val.begin(),
                      val.end(),
                      std::ostream_iterator<floatingpointT>(ofs, " "));

            ofs << '\n';
        }
    }
};
}

template<typename solverT, typename floatingpointT>
class pnm_observer
{
    details::pi_pnm_observer<floatingpointT> m_pi_obs;
    details::ap_pnm_observer<floatingpointT> m_ap_obs;

public:
    pnm_observer(const solverT& slv, std::string basename, int loop)
      : m_pi_obs(basename, slv.pi, slv.m, loop)
      , m_ap_obs(basename, slv.ap, slv.P, slv.m, slv.n)
    {}

    void make_observation()
    {
        m_pi_obs.make_observation();
        m_ap_obs.make_observation();
    }
};

template<typename solverT, typename floatingpointT>
class file_observer
{
    details::pi_file_observer<floatingpointT> m_pi_obs;
    details::ap_file_observer<floatingpointT> m_ap_obs;

public:
    file_observer(const solverT& slv, std::string basename, int loop)
      : m_pi_obs(basename, slv.pi, slv.m, loop)
      , m_ap_obs(basename, slv.ap, slv.P, slv.m, slv.n)
    {}

    void make_observation()
    {
        m_pi_obs.make_observation();
        m_ap_obs.make_observation();
    }
};

template<typename solverT, typename floatingpointT>
class none_observer
{
public:
    none_observer(const solverT& /*slv*/,
                  std::string /*basename*/,
                  int /*loop*/)
    {}

    void make_observation()
    {}
};

} // namespace baryonyx

#endif
