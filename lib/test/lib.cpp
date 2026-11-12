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

#include "bit-array.hpp"
#include "branch-and-bound-solver.hpp"
#include "exhaustive-solver.hpp"
#include "itm-common.hpp"
#include "itm-polish.hpp"
#include "itm-repair.hpp"
#include "itm-restart-optimizer.hpp"
#include "memory.hpp"
#include "pnm.hpp"
#include "sparse-matrix.hpp"
#include "parallel-for.hpp"
#include "random-distributions.hpp"
#include "random-engine.hpp"

#include <boost/ut.hpp>

#include <baryonyx/core>

#include <functional>
#include <cmath>
#include <numeric>

template<typename R>
std::tuple<baryonyx::real, int>
sum(const R& r, int r_size, int result)
{
    baryonyx::real sumr{ 0 };
    int sumf{ 0 };

    for (int i = 0, e = std::min(r_size - 1, result); i <= e; ++i) {
        sumr += r[i].value;
        sumf += r[i].f;
    }

    return std::make_tuple(sumr, sumf);
}

namespace {

/** Costs of the repair tests: all variables cost 1. */
struct unit_costs
{
    baryonyx::real operator()(int, const baryonyx::bit_array&) const noexcept
    {
        return 1;
    }
};

/** Linear costs of the polish tests. */
struct vector_costs
{
    std::vector<baryonyx::real> costs;

    baryonyx::real operator()(int j, const baryonyx::bit_array&) const noexcept
    {
        return costs[static_cast<std::size_t>(j)];
    }
};

bool
satisfies(const baryonyx::bit_array& x,
          const std::vector<baryonyx::itm::merged_constraint>& csts)
{
    for (const auto& cst : csts) {
        int value = 0;
        for (const auto& elem : cst.elements)
            value += elem.factor * x[elem.variable_index];

        if (value < cst.min || value > cst.max)
            return false;
    }

    return true;
}

} // namespace

int
main()
{
    using namespace boost::ut;

    "check_numeric_cast"_test = [] {
        int small_positive = 1;
        int small_negative = -1;
        int large_positive = std::numeric_limits<int>::max();
        int large_negative = std::numeric_limits<int>::min();

        using baryonyx::numeric_cast;

        expect(nothrow([&] { numeric_cast<signed char>(small_positive); }));
        expect(nothrow([&] { numeric_cast<signed char>(small_negative); }));
        expect(throws([&] { numeric_cast<signed char>(large_positive); }));
        expect(throws([&] { numeric_cast<signed char>(large_negative); }));

        expect(throws([&] { numeric_cast<unsigned char>(small_negative); }));
        expect(nothrow([&] { numeric_cast<unsigned int>(large_positive); }));
        expect(throws([&] { numeric_cast<unsigned int>(large_negative); }));
        expect(nothrow([&] { numeric_cast<long long>(large_negative); }));
        expect(throws([&] { numeric_cast<std::size_t>(small_negative); }));

        static_assert(numeric_cast<unsigned char>(255) == 255u);

        std::vector<int> v;

        expect(nothrow([&v] { baryonyx::numeric_cast<short int>(v.size()); }));
        expect(
          nothrow([&v] { baryonyx::numeric_cast<short int>(v.capacity()); }));
        expect(
          throws([&v] { baryonyx::numeric_cast<short int>(v.max_size()); }));

        auto checked_size = baryonyx::numeric_cast<unsigned int>(v.size());
        expect(0 == checked_size);
    };

    "test_subvector_api()"_test = [] {
        baryonyx::itm::shared_subvector s;
        s.reserve(128);

        expect(s.size() == 0);

        bool is_initialized = s.init(4);
        expect(is_initialized);

        auto id = s.emplace();
        expect(s.size() == 4);
        expect(id == 0);

        for (baryonyx::itm::shared_subvector::size_type
               i = id,
               e = id + s.element_size();
             i != e;
             ++i)
            s[i] = static_cast<std::int8_t>(i);

        auto id2 = s.emplace();
        expect(s.size() == 8);
        expect(id2 == 4);

        for (baryonyx::itm::shared_subvector::size_type
               i = id2,
               e = id2 + s.element_size();
             i != e;
             ++i)
            s[i] = static_cast<std::int8_t>(i);

        for (baryonyx::itm::shared_subvector::size_type i = 0, e = s.size();
             i != e;
             ++i)
            expect(s[i] == static_cast<int8_t>(i));

        s.remove(id);

        expect(s.size() == 8);

        auto id3 = s.emplace();
        expect(s.size() == 8);
        expect(id3 == 0);

        is_initialized = s.init(5);
        expect(is_initialized);
        expect(s.size() == 0);

        id = s.emplace();
        id2 = s.emplace();
        id3 = s.emplace();

        expect(s.size() == 3u * 5u);

        for (unsigned int i = 0, e = s.size(); i != e; ++i)
            s[i] = static_cast<std::int8_t>(i);

        std::int8_t cmp = 0;
        for (auto [i, e] = s.element(id); i != e; ++i, ++cmp)
            expect(cmp == s[i]);
        for (auto [i, e] = s.element(id2); i != e; ++i, ++cmp)
            expect(cmp == s[i]);
        for (auto [i, e] = s.element(id3); i != e; ++i, ++cmp)
            expect(cmp == s[i]);

        s.remove(id3);
        s.remove(id2);
        s.remove(id);

        // Elements of 5 bytes store the free list links at unaligned offsets.
        expect(s.emplace() == id);
        expect(s.emplace() == id2);
        expect(s.emplace() == id3);
        expect(s.size() == 3u * 5u);
    };

    "subvector_4096_16"_test = [] {
        const auto size = 4096u;
        const auto element_size = 16u;
        baryonyx::itm::shared_subvector s;
        s.reserve(size);

        expect(s.size() == 0);

        bool is_initialized = s.init(element_size);
        expect(is_initialized);

        unsigned element_total = size / element_size;

        for (unsigned int i = 0; i != element_total; ++i) {
            auto id = s.emplace();

            expect(id < size);
            for (auto [si, se] = s.element(id); si != se; ++si)
                s[si] = 99;
        }

        unsigned element_removed = 0;
        for (unsigned int i = 0; i != element_total / 2; ++i) {
            for (auto [si, se] = s.element(i * element_size); si != se; ++si) {
                expect(s[si] == 99);
                s[si] = 77;
            }

            s.remove(i * element_size);
            ++element_removed;
        }

        expect(element_removed == element_total / 2);

        for (unsigned int i = 0; i != element_removed; ++i) {
            auto id = s.emplace();

            expect(id < size);
            for (auto [si, se] = s.element(id); si != se; ++si)
                s[si] = 0;
        }

        for (unsigned int i = 0, e = s.size(); i != e; ++i)
            expect(s[i] == 99 || s[i] == 0);
    };

    "subvector_8192_512"_test = [] {
        const auto size = 8192u;
        const auto element_size = 512u;
        baryonyx::itm::shared_subvector s;
        s.reserve(size);

        expect(s.size() == 0);

        bool is_initialized = s.init(element_size);
        expect(is_initialized);

        unsigned element_total = size / element_size;

        for (unsigned int i = 0; i != element_total; ++i) {
            auto id = s.emplace();

            expect(id < size);
            for (auto [si, se] = s.element(id); si != se; ++si)
                s[si] = 99;
        }

        unsigned element_removed = 0;
        for (unsigned int i = 0; i != element_total / 2; ++i) {
            for (auto [si, se] = s.element(i * element_size); si != se; ++si) {
                expect(s[si] == 99);
                s[si] = 77;
            }

            s.remove(i * element_size);
            ++element_removed;
        }

        expect(element_removed == element_total / 2);

        for (unsigned int i = 0; i != element_removed; ++i) {
            auto id = s.emplace();

            expect(id < size);
            for (auto [si, se] = s.element(id); si != se; ++si)
                s[si] = 0;
        }

        for (unsigned int i = 0, e = s.size(); i != e; ++i)
            expect(s[i] == 99 || s[i] == 0);
    };

    "check_branch_and_bound_solver"_test = [] {
        struct rc
        {
            baryonyx::real value;
            int id;
            int f;
        };

        baryonyx::itm::branch_and_bound_solver<baryonyx::itm::maximize_tag> bb;
        bb.reserve(10u);

        {
            std::vector<rc> R = {
                { 7.0, 0, 13 }, { 4.0, 1, 12 }, { 3.0, 2, 8 }, { 3.0, 3, 10 }
            };
            const int nb = 4;
            const int b_min = 0;
            const int b_max = 30;

            auto result = bb.solve(R, nb, b_min, b_max);
            expect(result >= -1 && result <= static_cast<int>(R.size()));

            auto [r, f] = sum(R, nb, result);
            expect(f == 25);
            expect(r == 11.0);
        }

        {
            std::vector<rc> R = {
                { 16.0, 0, 2 }, { 19, 1, 3 }, { 23, 2, 4 }, { 26, 3, 5 }
            };

            const int nb = 4;
            const int b_min = 0;
            const int b_max = 7;

            auto result = bb.solve(R, nb, b_min, b_max);
            expect(result >= -1 && result <= static_cast<int>(R.size()));

            auto [r, f] = sum(R, nb, result);
            expect(f == 7);
            expect(r == 42.0);
        }

        {
            std::vector<rc> R = {
                { 1.0, 0, -1 }, { 1.0, 1, -1 }, { 1.0, 2, -1 }, { 1.0, 3, 1 }
            };

            const int nb = 4;
            const int b_min = -2;
            const int b_max = -2;

            auto result = bb.solve(R, nb, b_min, b_max);
            expect(result >= -1 && result <= static_cast<int>(R.size()));

            auto [r, f] = sum(R, nb, result);
            expect(f == -2);
            expect(r == 2.0); // Strange ?!?
        }

        baryonyx::itm::branch_and_bound_solver<baryonyx::itm::minimize_tag>
          bb_min;
        bb_min.reserve(16u);

        {
            std::vector<rc> R = {
                { 5.0, 0, 1 }, { 6.0, 1, 1 }, { 7.0, 2, 1 }, { 8.0, 3, 1 }
            };

            auto result = bb_min.solve(R, 4, 0, 4);
            expect(result == -1);
        }

        {
            constexpr auto inf =
              std::numeric_limits<baryonyx::real>::infinity();

            std::vector<rc> R = { { -1.4e36f, 11, -2 }, { 2.1e32f, 8, -3 },
                                  { 2.8e35f, 5, -3 },   { 1.7e36f, 4, -4 },
                                  { inf, 0, -3 },       { -2.6e36f, 3, 5 },
                                  { -5.6e35f, 1, 2 },   { -9.3e34f, 6, 1 },
                                  { -1.5e32f, 2, 1 },   { 3.7e32f, 10, 1 },
                                  { 2.3e37f, 7, 4 },    { inf, 9, 2 } };
            const int nb = static_cast<int>(R.size());

            auto result = bb_min.solve(R, nb, 12, 16);
            expect(result >= -1 && result <= nb);

            auto [r, f] = sum(R, nb, result);
            expect(12 <= f && f <= 16);
        }
    };

    "check_exhaustive_solver"_test = [] {
        struct rc
        {
            baryonyx::real value;
            int id;
            int f;
        };

        struct cst
        {
            int factor;
        };

        {
            // Only the empty selection satisfies 0 <= 2x + 3y <= 1.
            baryonyx::itm::exhaustive_solver<baryonyx::itm::minimize_tag> ex;
            ex.reserve(10u, 1u);

            std::vector<rc> R = { { -1.0, 0, 2 }, { -2.0, 1, 3 } };
            std::vector<cst> C = { { 2 }, { 3 } };

            ex.build_constraints(0, C, 0, 1);
            expect(ex.solve(0, R, 2) == -1);
        }

        {
            baryonyx::itm::exhaustive_solver<baryonyx::itm::maximize_tag> ex;
            ex.reserve(10u, 1u);

            std::vector<rc> R = {
                { 7.0, 0, 13 }, { 4.0, 1, 12 }, { 3.0, 2, 8 }, { 3.0, 3, 10 }
            };

            std::vector<cst> C = { { 13 }, { 12 }, { 8 }, { 10 } };

            const int nb = 4;
            const int b_min = 0;
            const int b_max = 30;

            ex.build_constraints(0, C, b_min, b_max);
            auto result = ex.solve(0, R, nb);
            expect(result >= -1 && result <= static_cast<int>(R.size()));

            auto [r, f] = sum(R, nb, result);
            expect(f == 25);
            expect(r == 11.0);
        }

        {
            baryonyx::itm::exhaustive_solver<baryonyx::itm::maximize_tag>
              ex;
            ex.reserve(10u, 1u);

            std::vector<rc> R = {
                { 16.0, 0, 2 }, { 19, 1, 3 }, { 23, 2, 4 }, { 26, 3, 5 }
            };

            std::vector<cst> C = { { 2 }, { 3 }, { 4 }, { 5 } };

            const int nb = 4;
            const int b_min = 0;
            const int b_max = 7;

            ex.build_constraints(0, C, b_min, b_max);
            auto result = ex.solve(0, R, nb);
            expect(result >= -1 && result <= static_cast<int>(R.size()));

            auto [r, f] = sum(R, nb, result);
            expect(f == 7);
            expect(r == 42.0);
        }

        {
            baryonyx::itm::exhaustive_solver<baryonyx::itm::maximize_tag>
              ex;
            ex.reserve(10u, 1u);

            std::vector<rc> R = {
                { 1.0, 0, -1 }, { 1.0, 1, -1 }, { 1.0, 2, -1 }, { 1.0, 3, 1 }
            };

            std::vector<cst> C = { { -1 }, { -1 }, { -1 }, { 1 } };

            const int nb = 4;
            const int b_min = -2;
            const int b_max = -2;

            ex.build_constraints(0, C, b_min, b_max);
            auto result = ex.solve(0, R, nb);
            expect(result >= -1 && result <= static_cast<int>(R.size()));

            auto [r, f] = sum(R, nb, result);
            expect(f == -2);
            expect(r == 4.0);
        }

        {
            baryonyx::itm::exhaustive_solver<baryonyx::itm::maximize_tag>
              ex;
            ex.reserve(10u, 2u);

            std::vector<rc> R1 = {
                { 1.0, 0, -1 }, { 1.0, 1, -1 }, { 1.0, 2, -1 }, { 1.0, 3, 1 }
            };
            std::vector<rc> R2 = { { 1.0, 0, -1 }, { -0.5, 1, -1 } };

            std::vector<cst> C1 = { { -1 }, { -1 }, { -1 }, { 1 } };
            std::vector<cst> C2 = { { -1 }, { -1 } };

            ex.build_constraints(0, C1, -2, -2);
            ex.build_constraints(1, C2, -2, -1);

            {
                auto result = ex.solve(0, R1, 4);
                expect(result >= -1 && result <= static_cast<int>(R1.size()));

                auto [r, f] = sum(R1, 4, result);
                expect(f == -2);
                expect(r == 4.0);
            }

            {
                auto result = ex.solve(1, R2, 2);
                expect(result >= -1 && result <= static_cast<int>(R2.size()));

                auto [r, f] = sum(R2, 2, result);
                expect(f == -1);
                expect(r == 1.0);
            }
        }
    };

    "check_generate_seed"_test = [] {
        baryonyx::random_engine rng(42);

        const auto seeds = baryonyx::itm::generate_seed(rng, 256u);
        expect(seeds.size() == 256u);

        auto sorted = seeds;
        std::ranges::sort(sorted);
        expect(std::ranges::adjacent_find(sorted) == sorted.end());

        expect(baryonyx::itm::generate_seed(rng, 0u).empty());
    };

    "check_sparse_matrix"_test = [] {
        struct cst
        {
            std::vector<baryonyx::function_element> elements;
        };

        // Row 0: 2x0 + 3x2, row 1: -x1 + 4x2. The matrix stores indices into
        // the row-ordered elements so both views must reference the same
        // cells.
        const std::vector<cst> csts{ { { { 2, 0 }, { 3, 2 } } },
                                     { { { -1, 1 }, { 4, 2 } } } };

        baryonyx::sparse_matrix<int> m(csts, 2, 3);
        expect(m.length() == 4);

        auto [rb, re] = m.row(1);
        expect(std::distance(rb, re) == 2);
        expect(rb->column == 1 && (rb + 1)->column == 2);

        auto [cb, ce] = m.column(2);
        expect(std::distance(cb, ce) == 2);
        expect(cb->row == 0 && (cb + 1)->row == 1);
        expect(cb->value == 1 && (cb + 1)->value == 3);
        expect(cb->factor == 3 && (cb + 1)->factor == 4);

        auto [c0b, c0e] = m.column(0);
        expect(std::distance(c0b, c0e) == 1);

        static_assert(sizeof(baryonyx::sparse_matrix<>::col_value) ==
                      2 * sizeof(int));

        const baryonyx::sparse_matrix<> unweighted(csts, 2, 3);
        auto [ub, ue] = unweighted.column(2);
        expect(std::distance(ub, ue) == 2);
        expect(ub->row == 0 && (ub + 1)->value == 3);

        const std::vector<cst> empty_rows(2);
        baryonyx::sparse_matrix<int> empty(empty_rows, 2, 3);
        expect(empty.length() == 0);
        auto [eb, ee] = empty.column(1);
        expect(eb == ee);
    };

    "check_xoshiro256pp"_test = [] {
        // Reference values from the authors' C code with splitmix64 seeding.
        baryonyx::xoshiro256pp rng(1);
        expect(rng() == 0xcfc5d07f6f03c29bULL);
        expect(rng() == 0xbf424132963fe08dULL);
        expect(rng() == 0x19a37d5757aaf520ULL);

        baryonyx::xoshiro256pp same(1), other(2);
        rng.seed(1);
        for (int i = 0; i != 100; ++i)
            expect(rng() == same());
        rng.seed(1);
        expect(rng() != other());
    };

    "check_random_distributions"_test = [] {
        // Reference values from an independent transcription of the
        // algorithms: these draws must not depend on the platform.
        {
            baryonyx::xoshiro256pp rng(1);
            expect(baryonyx::uniform_below(rng, 1000) == 387u);
            expect(baryonyx::uniform_below(rng, 1000) == 965u);
            expect(baryonyx::uniform_below(rng, 1000) == 744u);
        }

        {
            baryonyx::xoshiro256pp rng(1);
            std::vector<int> v(10);
            std::iota(v.begin(), v.end(), 0);
            baryonyx::shuffle(v.begin(), v.end(), rng);
            expect(v == std::vector<int>{ 1, 4, 5, 3, 6, 2, 9, 0, 8, 7 });
        }

        constexpr int draws = 200000;
        baryonyx::xoshiro256pp rng(42);

        {
            std::vector<int> buckets(10, 0);
            for (int i = 0; i != draws; ++i)
                ++buckets[baryonyx::uniform_below(rng, 10)];
            for (const int count : buckets)
                expect(std::abs(count - draws / 10) < draws / 100);
        }

        {
            baryonyx::bernoulli_distribution coin(0.3);
            int heads = 0;
            for (int i = 0; i != draws; ++i)
                heads += coin(rng) ? 1 : 0;
            expect(std::abs(heads - draws * 3 / 10) < draws / 100);

            baryonyx::bernoulli_distribution never(0.0), always(1.0);
            expect(!never(rng) && always(rng));
        }

        {
            baryonyx::uniform_int_distribution<int> dice(1, 6);
            bool seen_low = false, seen_high = false;
            for (int i = 0; i != 1000; ++i) {
                const int v = dice(rng);
                expect(1 <= v && v <= 6);
                seen_low = seen_low || v == 1;
                seen_high = seen_high || v == 6;
            }
            expect(seen_low && seen_high);
        }

        {
            baryonyx::uniform_real_distribution<double> real(-2.0, 3.0);
            for (int i = 0; i != 1000; ++i) {
                const double v = real(rng);
                expect(-2.0 <= v && v < 3.0);
            }
        }

        {
            baryonyx::normal_distribution<double> normal(1.0, 2.0);
            expect(normal.mean() == 1.0 && normal.stddev() == 2.0);

            double sum = 0, sum_sq = 0;
            for (int i = 0; i != draws; ++i) {
                const double v = normal(rng);
                sum += v;
                sum_sq += v * v;
            }
            const double mean = sum / draws;
            const double variance = sum_sq / draws - mean * mean;
            expect(std::abs(mean - 1.0) < 0.02);
            expect(std::abs(std::sqrt(variance) - 2.0) < 0.02);
        }
    };

    "check_parallel_for"_test = [] {
        for (const int threads : { 1, 2, 4 }) {
            baryonyx::parallel_for pool(threads);
            expect(pool.thread_number() == threads);

            std::vector<int> visits(1000, 0);
            auto body = [&visits](int i) { ++visits[i]; };

            for (int run = 0; run != 200; ++run)
                pool.run(static_cast<int>(visits.size()), body);

            expect(std::ranges::all_of(visits,
                                       [](int v) { return v == 200; }));

            int calls = 0;
            auto count = [&calls](int) { ++calls; };
            pool.run(0, count);
            pool.run(1, count);
            expect(calls == 1);
        }
    };

    "check_pnm_iterator"_test = [] {
        baryonyx::pnm_array pnm(1, 3);
        expect((!!pnm) >> fatal);

        auto first = pnm.begin();
        auto it = first;
        const auto old = it++;

        expect(old == first);
        expect(it != first);
        expect(first < it);

        *old = baryonyx::colormap::rgb{ 1, 2, 3 };
        *it++ = baryonyx::colormap::rgb{ 4, 5, 6 };
        *it++ = baryonyx::colormap::rgb{ 7, 8, 9 };
        expect(it == pnm.end());

        expect(pnm(0, 0)[0] == 1 && pnm(0, 0)[2] == 3);
        expect(pnm(0, 1)[0] == 4 && pnm(0, 1)[2] == 6);
        expect(pnm(0, 2)[0] == 7 && pnm(0, 2)[2] == 9);
    };

    "check_observer_pnm"_test = [] {
        {
            std::vector<double> v{ 0,   -1, -3,  -2, 0,    3, 3,
                                   1.5, 2,  0.5, 2,  -1.5, 1, 4 };

            baryonyx::pnm_array obs(7, 2);
            expect(!!obs);

            std::transform(
              v.begin(), v.end(), obs.begin(), baryonyx::colormap(-5.f, 5.f));

            obs("test.pnm");
        }

        {
            baryonyx::pnm_vector obs("test2.pnm", 4, 4);
            expect(!!obs);

            std::vector<double> v{ 0, 0.1, 0.2, 0.3 };

            for (int i = 0; i != 4; ++i) {
                std::transform(v.begin(),
                               v.end(),
                               obs.begin(),
                               baryonyx::colormap(-1.0f, 1.0f));
                obs.flush();

                std::ranges::transform(
                  v, v.begin(), [](double value) { return value + 0.1; });
            }
        }
    };

    "check_show_size"_test = [] {
        double size;
        baryonyx::show_size_type type;

        std::tie(size, type) = baryonyx::memory_consumed_size(0);
        expect(size == 0 && type == baryonyx::show_size_type::B);

        std::tie(size, type) = baryonyx::memory_consumed_size(1024);
        expect(size == 1 && type == baryonyx::show_size_type::KB);

        std::tie(size, type) = baryonyx::memory_consumed_size(2048);
        expect(size == 2 && type == baryonyx::show_size_type::KB);

        std::tie(size, type) = baryonyx::memory_consumed_size(2048);
        expect(size == 2 && type == baryonyx::show_size_type::KB);

        std::tie(size, type) = baryonyx::memory_consumed_size(1049088);
        expect(size > 1 && size < 2 && type == baryonyx::show_size_type::MB);

        std::tie(size, type) = baryonyx::memory_consumed_size(1049088);
        expect(size > 1 && size < 2 && type == baryonyx::show_size_type::MB);

        std::tie(size, type) = baryonyx::memory_consumed_size(17179869696);
        expect(size > 16 && size < 17 && type == baryonyx::show_size_type::GB);
    };

    "check_trim_functions"_test = [] {
        {
            std::string str = "    xxxx    ";
            expect(baryonyx::left_trim(str) == "xxxx    ");
            expect(baryonyx::right_trim(str) == "    xxxx");
            expect(baryonyx::trim(str) == "xxxx");
            expect(str == "    xxxx    ");
        }

        {
            std::string str = "xxxx    ";
            expect(baryonyx::left_trim(str) == "xxxx    ");
            expect(baryonyx::right_trim(str) == "xxxx");
            expect(baryonyx::trim(str) == "xxxx");
            expect(str == "xxxx    ");
        }

        {
            std::string str = "    xxxx";
            expect(baryonyx::left_trim(str) == "xxxx");
            expect(baryonyx::right_trim(str) == "    xxxx");
            expect(baryonyx::trim(str) == "xxxx");
            expect(str == "    xxxx");
        }

        {
            std::string str = "    ";
            expect(baryonyx::left_trim(str).empty());
            expect(baryonyx::right_trim(str).empty());
            expect(baryonyx::trim(str).empty());
            expect(str == "    ");
        }

        {
            std::string str = "";
            expect(baryonyx::left_trim(str).empty());
            expect(baryonyx::right_trim(str).empty());
            expect(baryonyx::trim(str).empty());
            expect(str.empty());
        }

        {
            std::string str = "abc";
            expect(baryonyx::left_trim(str) == "abc");
            expect(baryonyx::right_trim(str) == "abc");
            expect(baryonyx::trim(str) == "abc");
            expect(str == "abc");
        }
    };

    "check_bit_array"_test = [] {
        baryonyx::bit_array a(100);

        expect(a.size() == 100);
        expect(a.block_size() ==
               1 + (100 / baryonyx::bit_array::bit_per_block));

        a.zeros();
        for (auto i = 0; i != a.size(); ++i)
            expect(a.get(i) == 0);

        a.ones();
        for (auto i = 0; i != a.size(); ++i)
            expect(a.get(i) == 1);

        for (auto i = 0; i != a.size(); ++i)
            a.unset(i);
        for (auto i = 0; i != a.size(); ++i)
            expect(a.get(i) == 0);

        for (auto i = 0; i != a.size(); ++i)
            a.set(i);
        for (auto i = 0; i != a.size(); ++i)
            expect(a.get(i) == 1);

        for (auto i = 0; i != a.size(); i += 2)
            a.unset(i);
        for (auto i = 0; i != a.size(); i += 2)
            expect(a.get(i) == 0 && a.get(i + 1) == 1);

        for (auto i = 0; i != a.size(); ++i)
            a.invert(i);

        for (auto i = 0; i != a.size(); i += 2)
            expect(a.get(i) == 1 && a.get(i + 1) == 0);

        for (auto i = 0; i != a.size(); ++i) {
            a.set(i, i % 3 == 0);
            expect(a.get(i) == (i % 3 == 0 ? 1 : 0));
        }

        constexpr auto bits = baryonyx::bit_array::bit_per_block;
        expect(baryonyx::bit_array(bits * 2).block_size() == 2);
        expect(baryonyx::bit_array(bits * 2 + 1).block_size() == 3);

        // Bits past size() stay zero, so equal contents compare and hash
        // equal whatever the way they were written.
        baryonyx::bit_array by_ones(100), by_set(100);
        by_ones.ones();
        for (auto i = 0; i != by_set.size(); ++i)
            by_set.set(i);
        expect(by_ones == by_set);
        expect(baryonyx::bit_array_hash()(by_ones) ==
               baryonyx::bit_array_hash()(by_set));

        by_set.set_block(1, baryonyx::bit_array::k_ones);
        expect(by_ones == by_set);

        baryonyx::bit_array moved(std::move(by_ones));
        expect(moved == by_set);
        expect(by_ones.size() == 0 && by_ones.block_size() == 0);
    };

    "check_repair"_test = [] {
        using baryonyx::itm::merged_constraint;
        using fe = baryonyx::function_element;

        // Set partitioning: each variable covers two of the four rows.
        const std::vector<merged_constraint> spp{
            { { fe(1, 0), fe(1, 1), fe(1, 2) }, 1, 1, 0 },
            { { fe(1, 0), fe(1, 3), fe(1, 4) }, 1, 1, 1 },
            { { fe(1, 1), fe(1, 3), fe(1, 5) }, 1, 1, 2 },
            { { fe(1, 2), fe(1, 4), fe(1, 5) }, 1, 1, 3 },
        };

        for (std::uint64_t seed = 1; seed != 20; ++seed) {
            baryonyx::random_engine rng(seed);
            baryonyx::bit_array x(6);
            x.set(0);
            x.set(1);
            x.set(5);

            expect(!satisfies(x, spp));
            expect(baryonyx::itm::repair(x, spp, unit_costs{}, 6, rng, 10000) ==
                   0);
            expect(satisfies(x, spp));
        }

        // Negative and integer coefficients, the only solution is
        // x = (1, 1, 1, 0).
        const std::vector<merged_constraint> z{
            { { fe(2, 0), fe(3, 1), fe(-1, 2) }, 4, 4, 0 },
            { { fe(1, 0), fe(1, 3) }, 0, 1, 1 },
            { { fe(-1, 1), fe(1, 3) }, -1, 0, 2 },
        };

        for (std::uint64_t seed = 1; seed != 20; ++seed) {
            baryonyx::random_engine rng(seed);
            baryonyx::bit_array x(4);
            x.set(2);
            x.set(3);

            expect(baryonyx::itm::repair(x, z, unit_costs{}, 4, rng, 10000) ==
                   0);
            expect(satisfies(x, z));
        }

        // An unreachable bound stops at the flip limit.
        const std::vector<merged_constraint> impossible{
            { { fe(1, 0), fe(1, 1) }, 3, 3, 0 },
        };

        baryonyx::random_engine rng(7);
        baryonyx::bit_array x(2);
        expect(baryonyx::itm::repair(x, impossible, unit_costs{}, 2, rng, 500) ==
               1);

        expect(baryonyx::itm::is_repairable(1, 50));
        expect(baryonyx::itm::is_repairable(10, 50));
        expect(!baryonyx::itm::is_repairable(11, 50));
        expect(!baryonyx::itm::is_repairable(0, 50));
        expect(baryonyx::itm::is_repairable(100, 10000));
        expect(!baryonyx::itm::is_repairable(101, 10000));
    };

    "check_portfolio"_test = [] {
        using baryonyx::itm::portfolio_variant_index;
        using baryonyx::itm::portfolio_variants;

        const baryonyx::solver_parameters user;
        const int variants = baryonyx::length(portfolio_variants);

        expect(portfolio_variant_index(user, 0) == -1);

        std::vector<int> seen;
        for (unsigned thread = 1; thread <= unsigned(variants); ++thread)
            seen.emplace_back(portfolio_variant_index(user, thread));
        std::sort(seen.begin(), seen.end());
        expect(std::adjacent_find(seen.begin(), seen.end()) == seen.end());
        expect(seen.front() == 0 && seen.back() == variants - 1);

        expect(portfolio_variant_index(user, unsigned(variants) + 1) ==
               portfolio_variant_index(user, 1));

        // A variant equal to the user's parameters is skipped.
        auto loo = user;
        loo.cost_norm = baryonyx::solver_parameters::cost_norm_type::loo;
        for (unsigned thread = 1; thread != 20; ++thread) {
            const auto variant = portfolio_variant_index(loo, thread);
            expect(variant >= 0);
            expect(portfolio_variants[variant].name != "norm=loo");
        }
    };

    "check_polish"_test = [] {
        using baryonyx::itm::merged_constraint;
        using baryonyx::itm::polisher;
        using fe = baryonyx::function_element;

        // A single flip: x1 is useless to cover x0 + x1 >= 1.
        {
            const std::vector<merged_constraint> cover{
                { { fe(1, 0), fe(1, 1) }, 1, 2, 0 },
            };
            baryonyx::bit_array x(2);
            x.set(0);
            x.set(1);

            polisher polish(cover, 2);
            expect(polish(x, vector_costs{ { 1, 2 } }) == 1);
            expect(x[0] == 1 && x[1] == 0);
            expect(satisfies(x, cover));
        }

        // Only an exchange improves a partition: x1 replaces x0.
        {
            const std::vector<merged_constraint> partition{
                { { fe(1, 0), fe(1, 1), fe(1, 2) }, 1, 1, 0 },
                { { fe(1, 3), fe(1, 4) }, 1, 1, 1 },
            };
            baryonyx::bit_array x(5);
            x.set(0);
            x.set(3);

            polisher polish(partition, 5);
            expect(polish(x, vector_costs{ { 3, 1, 2, 1, 4 } }) == 1);
            expect(x[0] == 0 && x[1] == 1 && x[2] == 0);
            expect(x[3] == 1 && x[4] == 0);
            expect(satisfies(x, partition));

            // A local optimum stays unchanged.
            expect(polish(x, vector_costs{ { 3, 1, 2, 1, 4 } }) == 0);
            expect(x[1] == 1 && x[3] == 1);
        }

        // Quadratic costs: x0 and x1 are cheap alone, expensive together,
        // and the x2 * x2 term makes x2 worth taking.
        {
            baryonyx::objective_function obj;
            obj.elements.emplace_back(-1.0, 0);
            obj.elements.emplace_back(-1.0, 1);
            obj.elements.emplace_back(1.0, 2);
            obj.qelements.emplace_back(3.0, 0, 1);
            obj.qelements.emplace_back(-3.0, 2, 2);

            const baryonyx::itm::quadratic_cost_type costs(obj, 3);
            const std::vector<merged_constraint> none;
            baryonyx::bit_array x(3);
            x.set(0);
            x.set(1);

            polisher polish(none, 3);
            expect(polish(x, costs) == 2);
            expect(x[0] + x[1] == 1);
            expect(x[2] == 1);
            expect(costs.results(x, 0) == -3);
        }

        // Only the changed variables and their neighbours are checked.
        {
            const std::vector<merged_constraint> rows{
                { { fe(1, 0), fe(1, 1) }, 1, 2, 0 },
                { { fe(1, 2), fe(1, 3) }, 1, 2, 1 },
            };
            baryonyx::bit_array x(4);
            x.set(0);
            x.set(1);
            x.set(2);
            x.set(3);

            polisher polish(rows, 4);
            const std::vector<int> changed{ 0 };
            expect(polish(x, vector_costs{ { 1, 2, 1, 2 } }, changed) == 1);
            expect(x[1] == 0 && x[3] == 1);
        }
    };
}
