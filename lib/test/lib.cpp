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

#include "bit-array.hpp"
#include "branch-and-bound-solver.hpp"
#include "exhaustive-solver.hpp"
#include "fixed-2darray.hpp"
#include "fixed-array.hpp"
#include "memory.hpp"
#include "pnm.hpp"

#include <boost/ut.hpp>

#include <baryonyx/core>

#include <functional>
#include <numeric>

template<typename R>
std::tuple<double, int>
sum(const R& r, int r_size, int result)
{
    double sumr{ 0 };
    int sumf{ 0 };

    for (int i = 0, e = std::min(r_size - 1, result); i <= e; ++i) {
        sumr += r[i].value;
        sumf += r[i].f;
    }

    return std::make_tuple(sumr, sumf);
}

int
main()
{
    using namespace boost::ut;

    "check_numeric_cast"_test = [] {
        int small_positive = 1;
        int small_negative = -1;
        int large_positive = std::numeric_limits<int>::max();
        int large_negative = std::numeric_limits<int>::min();

        expect(baryonyx::is_numeric_castable<signed char>(small_positive));
        expect(baryonyx::is_numeric_castable<signed char>(small_negative));
        expect(!baryonyx::is_numeric_castable<signed char>(large_positive));
        expect(!baryonyx::is_numeric_castable<signed char>(large_negative));

        expect(baryonyx::is_numeric_castable<unsigned char>(small_positive));
        expect(!baryonyx::is_numeric_castable<unsigned char>(small_negative));
        expect(!baryonyx::is_numeric_castable<unsigned char>(large_positive));
        expect(!baryonyx::is_numeric_castable<unsigned char>(large_negative));

        expect(baryonyx::is_numeric_castable<signed int>(small_positive));
        expect(baryonyx::is_numeric_castable<signed int>(small_negative));
        expect(baryonyx::is_numeric_castable<signed int>(large_positive));
        expect(baryonyx::is_numeric_castable<signed int>(large_negative));

        expect(baryonyx::is_numeric_castable<unsigned int>(small_positive));
        expect(!baryonyx::is_numeric_castable<unsigned int>(small_negative));
        expect(baryonyx::is_numeric_castable<unsigned int>(large_positive));
        expect(!baryonyx::is_numeric_castable<unsigned int>(large_negative));

        expect(baryonyx::is_numeric_castable<long long>(small_positive));
        expect(baryonyx::is_numeric_castable<long long>(large_negative));
        expect(baryonyx::is_numeric_castable<long long>(small_positive));
        expect(baryonyx::is_numeric_castable<long long>(large_negative));

        expect(
          baryonyx::is_numeric_castable<unsigned long long>(small_positive));
        expect(
          !baryonyx::is_numeric_castable<unsigned long long>(small_negative));
        expect(
          baryonyx::is_numeric_castable<unsigned long long>(large_positive));
        expect(
          !baryonyx::is_numeric_castable<unsigned long long>(large_negative));

        expect(!baryonyx::is_numeric_castable<size_t>(small_negative));
        expect(!baryonyx::is_numeric_castable<size_t>(large_negative));

        std::vector<int> v;

        expect(nothrow([&v] { baryonyx::numeric_cast<short int>(v.size()); }));
        expect(
          nothrow([&v] { baryonyx::numeric_cast<short int>(v.capacity()); }));
        expect(
          throws([&v] { baryonyx::numeric_cast<short int>(v.max_size()); }));

        auto checked_size = baryonyx::numeric_cast<unsigned int>(v.size());
        expect(0 == checked_size);
    };

    "check_fixed_array"_test = [] {
        baryonyx::fixed_array<int> a(10);

        expect(a.size() == 10);

        std::iota(a.begin(), a.end(), 1);

        expect(a[0] == 1);
        expect(a[1] == 2);
        expect(a[2] == 3);
        expect(a[3] == 4);
        expect(a[4] == 5);
        expect(a[5] == 6);
        expect(a[6] == 7);
        expect(a[7] == 8);
        expect(a[8] == 9);
        expect(a[9] == 10);

        {
            baryonyx::fixed_array<int> copy(a);

            std::iota(copy.rbegin(), copy.rend(), 1);
            expect(copy[9] == 1);
            expect(copy[8] == 2);
            expect(copy[7] == 3);
            expect(copy[6] == 4);
            expect(copy[5] == 5);
            expect(copy[4] == 6);
            expect(copy[3] == 7);
            expect(copy[2] == 8);
            expect(copy[1] == 9);
            expect(copy[0] == 10);
        }

        baryonyx::fixed_array<int> b(a);

        expect(a.data() != b.data());

        baryonyx::fixed_array<int> c(std::move(a));

        expect(a.data() == nullptr);
        expect(b.data() != c.data());

        baryonyx::fixed_array<double> d(15, 3.0);

        expect(d[0] == 3.0);
        expect(d[7] == 3.0);
        expect(d[14] == 3.0);

        baryonyx::fixed_array<double> e;

        std::swap(d, e);

        expect(!d);
        expect(e[0] == 3.0);
        expect(e[7] == 3.0);
        expect(e[14] == 3.0);

        baryonyx::fixed_array<double> x(1000, 123.0);
        if (x) {
            for (int i{ 0 }; i < 1000; ++i)
                expect(x[i] == 123.0);

            auto it = std::find_if_not(x.begin(),
                                       x.end(),
                                       std::bind(std::equal_to<double>(),
                                                 123.0,
                                                 std::placeholders::_1));

            expect(it == x.end());
        }
    };

    "check_fixed_2darray"_test = [] {
        baryonyx::fixed_2darray<int> a(static_cast<size_t>(2),
                                       static_cast<size_t>(10));

        expect(a.size() == 20);
        expect(a.rows() == 2);
        expect(a.columns() == 10);

        std::iota(a.begin(), a.end(), 0);

        expect(a.data()[0] == 0);
        expect(a.data()[19] == 19);

        std::iota(a.rbegin(), a.rend(), 0);

        expect(a.data()[0] == 19);
        expect(a.data()[19] == 0);

        expect(a(0, 0) == 19);
        expect(a(1, 9) == 0);

        a(0, 1) = 100;
        a(1, 0) = 200;

        expect(a.data()[1] == 100);
        expect(a.data()[10] == 200);
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
            double value;
            int id;
            int f;
        };

        baryonyx::itm::branch_and_bound_solver<baryonyx::itm::maximize_tag,
                                               double>
          bb;
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
    };

    "check_exhaustive_solver"_test = [] {
        struct rc
        {
            double value;
            int id;
            int f;
        };

        struct cst
        {
            int factor;
        };

        {
            baryonyx::itm::exhaustive_solver<baryonyx::itm::maximize_tag,
                                             double>
              ex;
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
            baryonyx::itm::exhaustive_solver<baryonyx::itm::maximize_tag,
                                             double>
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
            baryonyx::itm::exhaustive_solver<baryonyx::itm::maximize_tag,
                                             double>
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
            baryonyx::itm::exhaustive_solver<baryonyx::itm::maximize_tag,
                                             double>
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

                std::transform(
                  v.begin(),
                  v.end(),
                  v.begin(),
                  std::bind(std::plus<double>(), 0.1, std::placeholders::_1));
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

        baryonyx::bit_array b(baryonyx::bit_array::bit_per_block * 2);
        baryonyx::bit_array c(baryonyx::bit_array::bit_per_block * 2);

        constexpr auto middle = baryonyx::bit_array::bit_per_block / 2;

        {
            b.ones();
            c.zeros();

            expect(b.block(0) == baryonyx::bit_array::k_ones);
            expect(b.block(1) == baryonyx::bit_array::k_ones);
            expect(c.block(0) == baryonyx::bit_array::k_zeros);
            expect(c.block(1) == baryonyx::bit_array::k_zeros);

            b.assign(c,
                     baryonyx::bit_array::bit_per_block / 2,
                     baryonyx::bit_array::bit_per_block * 3 / 2);

            expect(b.block(0) == baryonyx::bit_array::k_ones << middle);
            expect(b.block(1) == baryonyx::bit_array::k_ones >> middle);
        }

        {
            b.ones();
            c.zeros();

            expect(b.block(0) == baryonyx::bit_array::k_ones);
            expect(b.block(1) == baryonyx::bit_array::k_ones);
            expect(c.block(0) == baryonyx::bit_array::k_zeros);
            expect(c.block(1) == baryonyx::bit_array::k_zeros);

            b.assign(c, 0, baryonyx::bit_array::bit_per_block / 2);

            expect(b.block(0) == baryonyx::bit_array::k_ones >> middle);
            expect(b.block(1) == baryonyx::bit_array::k_ones);
        }

        {
            b.ones();
            c.zeros();

            expect(b.block(0) == baryonyx::bit_array::k_ones);
            expect(b.block(1) == baryonyx::bit_array::k_ones);
            expect(c.block(0) == baryonyx::bit_array::k_zeros);
            expect(c.block(1) == baryonyx::bit_array::k_zeros);

            b.assign(c,
                     baryonyx::bit_array::bit_per_block / 2,
                     baryonyx::bit_array::bit_per_block * 2);

            expect(b.block(0) == baryonyx::bit_array::k_ones << middle);
            expect(b.block(1) == baryonyx::bit_array::k_zeros);
        }

        {
            b.ones();
            c.zeros();

            expect(b.block(0) == baryonyx::bit_array::k_ones);
            expect(b.block(1) == baryonyx::bit_array::k_ones);
            expect(c.block(0) == baryonyx::bit_array::k_zeros);
            expect(c.block(1) == baryonyx::bit_array::k_zeros);

            b.assign(c,
                     baryonyx::bit_array::bit_per_block * 2 - 1,
                     baryonyx::bit_array::bit_per_block * 2);

            expect(b.block(0) == baryonyx::bit_array::k_ones);
            expect(b.block(1) == baryonyx::bit_array::k_ones - 1);
        }
    };
}
