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
#include "knapsack-dp-solver.hpp"
#include "memory.hpp"
#include "pnm.hpp"
#include "unit-test.hpp"

#include <baryonyx/core>

#include <fmt/printf.h>

#include <functional>
#include <numeric>

static void
check_print_api()
{
    fmt::print("\n\nOnly 16 colors:\n");
    fmt::print(fg(fmt::terminal_color::black), "A test of the color black\n");
    fmt::print(fg(fmt::terminal_color::red), "A test of the color red\n");
    fmt::print(fg(fmt::terminal_color::green), "A test of the color green\n");
    fmt::print(fg(fmt::terminal_color::yellow),
               "A test of the color yellow\n");
    fmt::print(fg(fmt::terminal_color::blue), "A test of the color blue\n");
    fmt::print(fg(fmt::terminal_color::magenta),
               "A test of the color magenta\n");
    fmt::print(fg(fmt::terminal_color::cyan), "A test of the color cyan\n");
    fmt::print(fg(fmt::terminal_color::white), "A test of the color white\n");
    fmt::print(fg(fmt::terminal_color::bright_black),
               "A test of the color bright_black\n");
    fmt::print(fg(fmt::terminal_color::bright_red),
               "A test of the color bright_red\n");
    fmt::print(fg(fmt::terminal_color::bright_green),
               "A test of the color bright_green\n");
    fmt::print(fg(fmt::terminal_color::bright_yellow),
               "A test of the color bright_yellow\n");
    fmt::print(fg(fmt::terminal_color::bright_blue),
               "A test of the color bright_blue\n");
    fmt::print(fg(fmt::terminal_color::bright_magenta),
               "A test of the color bright_magenta\n");
    fmt::print(fg(fmt::terminal_color::bright_cyan),
               "A test of the color bright_cyan\n");
    fmt::print(fg(fmt::terminal_color::bright_white),
               "A test of the color bright_white\n");

    fmt::print("\n\nBold + 16 colors:\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::black),
               "A test of the color black\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::red),
               "A test of the color red\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::green),
               "A test of the color green\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::yellow),
               "A test of the color yellow\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::blue),
               "A test of the color blue\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::magenta),
               "A test of the color magenta\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::cyan),
               "A test of the color cyan\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::white),
               "A test of the color white\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::bright_black),
               "A test of the color bright_black\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::bright_red),
               "A test of the color bright_red\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::bright_green),
               "A test of the color bright_green\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::bright_yellow),
               "A test of the color bright_yellow\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::bright_blue),
               "A test of the color bright_blue\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::bright_magenta),
               "A test of the color bright_magenta\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::bright_cyan),
               "A test of the color bright_cyan\n");
    fmt::print(fmt::emphasis::bold | fg(fmt::terminal_color::bright_white),
               "A test of the color bright_white\n");

    fmt::print("\n\nItalic + 16 colors:\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::black),
               "A test of the color black\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::red),
               "A test of the color red\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::green),
               "A test of the color green\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::yellow),
               "A test of the color yellow\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::blue),
               "A test of the color blue\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::magenta),
               "A test of the color magenta\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::cyan),
               "A test of the color cyan\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::white),
               "A test of the color white\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::bright_black),
               "A test of the color bright_black\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::bright_red),
               "A test of the color bright_red\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::bright_green),
               "A test of the color bright_green\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::bright_yellow),
               "A test of the color bright_yellow\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::bright_blue),
               "A test of the color bright_blue\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::bright_magenta),
               "A test of the color bright_magenta\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::bright_cyan),
               "A test of the color bright_cyan\n");
    fmt::print(fmt::emphasis::italic | fg(fmt::terminal_color::bright_white),
               "A test of the color bright_white\n");

    fmt::print("\n\nUnderline + 16 colors:\n");
    fmt::print(fmt::emphasis::underline | fg(fmt::terminal_color::black),
               "A test of the color black\n");
    fmt::print(fmt::emphasis::underline | fg(fmt::terminal_color::red),
               "A test of the color red\n");
    fmt::print(fmt::emphasis::underline | fg(fmt::terminal_color::green),
               "A test of the color green\n");
    fmt::print(fmt::emphasis::underline | fg(fmt::terminal_color::yellow),
               "A test of the color yellow\n");
    fmt::print(fmt::emphasis::underline | fg(fmt::terminal_color::blue),
               "A test of the color blue\n");
    fmt::print(fmt::emphasis::underline | fg(fmt::terminal_color::magenta),
               "A test of the color magenta\n");
    fmt::print(fmt::emphasis::underline | fg(fmt::terminal_color::cyan),
               "A test of the color cyan\n");
    fmt::print(fmt::emphasis::underline | fg(fmt::terminal_color::white),
               "A test of the color white\n");
    fmt::print(fmt::emphasis::underline |
                 fg(fmt::terminal_color::bright_black),
               "A test of the color bright_black\n");
    fmt::print(fmt::emphasis::underline | fg(fmt::terminal_color::bright_red),
               "A test of the color bright_red\n");
    fmt::print(fmt::emphasis::underline |
                 fg(fmt::terminal_color::bright_green),
               "A test of the color bright_green\n");
    fmt::print(fmt::emphasis::underline |
                 fg(fmt::terminal_color::bright_yellow),
               "A test of the color bright_yellow\n");
    fmt::print(fmt::emphasis::underline | fg(fmt::terminal_color::bright_blue),
               "A test of the color bright_blue\n");
    fmt::print(fmt::emphasis::underline |
                 fg(fmt::terminal_color::bright_magenta),
               "A test of the color bright_magenta\n");
    fmt::print(fmt::emphasis::underline | fg(fmt::terminal_color::bright_cyan),
               "A test of the color bright_cyan\n");
    fmt::print(fmt::emphasis::underline |
                 fg(fmt::terminal_color::bright_white),
               "A test of the color bright_white\n");
}

static void
check_numeric_cast()
{
    int small_positive = 1;
    int small_negative = -1;
    int large_positive = std::numeric_limits<int>::max();
    int large_negative = std::numeric_limits<int>::min();

    Ensures(baryonyx::is_numeric_castable<signed char>(small_positive));
    Ensures(baryonyx::is_numeric_castable<signed char>(small_negative));
    Ensures(!baryonyx::is_numeric_castable<signed char>(large_positive));
    Ensures(!baryonyx::is_numeric_castable<signed char>(large_negative));

    Ensures(baryonyx::is_numeric_castable<unsigned char>(small_positive));
    Ensures(!baryonyx::is_numeric_castable<unsigned char>(small_negative));
    Ensures(!baryonyx::is_numeric_castable<unsigned char>(large_positive));
    Ensures(!baryonyx::is_numeric_castable<unsigned char>(large_negative));

    Ensures(baryonyx::is_numeric_castable<signed int>(small_positive));
    Ensures(baryonyx::is_numeric_castable<signed int>(small_negative));
    Ensures(baryonyx::is_numeric_castable<signed int>(large_positive));
    Ensures(baryonyx::is_numeric_castable<signed int>(large_negative));

    Ensures(baryonyx::is_numeric_castable<unsigned int>(small_positive));
    Ensures(!baryonyx::is_numeric_castable<unsigned int>(small_negative));
    Ensures(baryonyx::is_numeric_castable<unsigned int>(large_positive));
    Ensures(!baryonyx::is_numeric_castable<unsigned int>(large_negative));

    Ensures(baryonyx::is_numeric_castable<long long>(small_positive));
    Ensures(baryonyx::is_numeric_castable<long long>(large_negative));
    Ensures(baryonyx::is_numeric_castable<long long>(small_positive));
    Ensures(baryonyx::is_numeric_castable<long long>(large_negative));

    Ensures(baryonyx::is_numeric_castable<unsigned long long>(small_positive));
    Ensures(
      !baryonyx::is_numeric_castable<unsigned long long>(small_negative));
    Ensures(baryonyx::is_numeric_castable<unsigned long long>(large_positive));
    Ensures(
      !baryonyx::is_numeric_castable<unsigned long long>(large_negative));

    Ensures(!baryonyx::is_numeric_castable<size_t>(small_negative));
    Ensures(!baryonyx::is_numeric_castable<size_t>(large_negative));

    std::vector<int> v;

    EnsuresNotThrow(baryonyx::numeric_cast<short int>(v.size()),
                    std::exception);

    EnsuresNotThrow(baryonyx::numeric_cast<short int>(v.capacity()),
                    std::exception);

    EnsuresThrow(baryonyx::numeric_cast<short int>(v.max_size()),
                 std::exception);

    auto checked_size = baryonyx::numeric_cast<unsigned int>(v.size());
    Ensures(0 == checked_size);
}

static void
check_fixed_array()
{
    baryonyx::fixed_array<int> a(10);

    Ensures(a.size() == 10);

    std::iota(a.begin(), a.end(), 1);

    Ensures(a[0] == 1);
    Ensures(a[1] == 2);
    Ensures(a[2] == 3);
    Ensures(a[3] == 4);
    Ensures(a[4] == 5);
    Ensures(a[5] == 6);
    Ensures(a[6] == 7);
    Ensures(a[7] == 8);
    Ensures(a[8] == 9);
    Ensures(a[9] == 10);

    {
        baryonyx::fixed_array<int> copy(a);

        std::iota(copy.rbegin(), copy.rend(), 1);
        Ensures(copy[9] == 1);
        Ensures(copy[8] == 2);
        Ensures(copy[7] == 3);
        Ensures(copy[6] == 4);
        Ensures(copy[5] == 5);
        Ensures(copy[4] == 6);
        Ensures(copy[3] == 7);
        Ensures(copy[2] == 8);
        Ensures(copy[1] == 9);
        Ensures(copy[0] == 10);
    }

    baryonyx::fixed_array<int> b(a);

    Ensures(a.data() != b.data());

    baryonyx::fixed_array<int> c(std::move(a));

    Ensures(a.data() == nullptr);
    Ensures(b.data() != c.data());

    baryonyx::fixed_array<double> d(15, 3.0);

    Ensures(d[0] == 3.0);
    Ensures(d[7] == 3.0);
    Ensures(d[14] == 3.0);

    baryonyx::fixed_array<double> e;

    std::swap(d, e);

    Ensures(!d);
    Ensures(e[0] == 3.0);
    Ensures(e[7] == 3.0);
    Ensures(e[14] == 3.0);

    baryonyx::fixed_array<double> x(1000, 123.0);
    if (x) {
        for (int i{ 0 }; i < 1000; ++i)
            Ensures(x[i] == 123.0);

        auto it = std::find_if_not(
          x.begin(),
          x.end(),
          std::bind(std::equal_to<double>(), 123.0, std::placeholders::_1));

        Ensures(it == x.end());
    }
}

static void
check_fixed_2darray()
{
    baryonyx::fixed_2darray<int> a(static_cast<size_t>(2),
                                   static_cast<size_t>(10));

    Ensures(a.size() == 20);
    Ensures(a.rows() == 2);
    Ensures(a.columns() == 10);

    std::iota(a.begin(), a.end(), 0);

    Ensures(a.data()[0] == 0);
    Ensures(a.data()[19] == 19);

    std::iota(a.rbegin(), a.rend(), 0);

    Ensures(a.data()[0] == 19);
    Ensures(a.data()[19] == 0);

    Ensures(a(0, 0) == 19);
    Ensures(a(1, 9) == 0);

    a(0, 1) = 100;
    a(1, 0) = 200;

    Ensures(a.data()[1] == 100);
    Ensures(a.data()[10] == 200);
}

void
test_subvector_api()
{
    baryonyx::itm::shared_subvector s;
    s.reserve(128);

    Ensures(s.size() == 0);

    bool is_initialized = s.init(4);
    Ensures(is_initialized);

    auto id = s.emplace();
    Ensures(s.size() == 4);
    Ensures(id == 0);

    for (baryonyx::itm::shared_subvector::size_type i = id,
                                                    e = id + s.element_size();
         i != e;
         ++i)
        s[i] = static_cast<std::int8_t>(i);

    auto id2 = s.emplace();
    Ensures(s.size() == 8);
    Ensures(id2 == 4);

    for (baryonyx::itm::shared_subvector::size_type i = id2,
                                                    e = id2 + s.element_size();
         i != e;
         ++i)
        s[i] = static_cast<std::int8_t>(i);

    for (baryonyx::itm::shared_subvector::size_type i = 0, e = s.size();
         i != e;
         ++i)
        Ensures(s[i] == static_cast<int8_t>(i));

    s.remove(id);

    Ensures(s.size() == 8);

    auto id3 = s.emplace();
    Ensures(s.size() == 8);
    Ensures(id3 == 0);

    is_initialized = s.init(5);
    Ensures(is_initialized);
    Ensures(s.size() == 0);

    id = s.emplace();
    id2 = s.emplace();
    id3 = s.emplace();

    Ensures(s.size() == 3u * 5u);

    for (unsigned int i = 0, e = s.size(); i != e; ++i)
        s[i] = static_cast<std::int8_t>(i);

    std::int8_t cmp = 0;
    for (auto [i, e] = s.element(id); i != e; ++i, ++cmp)
        Ensures(cmp == s[i]);
    for (auto [i, e] = s.element(id2); i != e; ++i, ++cmp)
        Ensures(cmp == s[i]);
    for (auto [i, e] = s.element(id3); i != e; ++i, ++cmp)
        Ensures(cmp == s[i]);

    s.remove(id3);
    s.remove(id2);
    s.remove(id);
}

void
test_size_type_greater_than_int8_subvector_impl(const unsigned size,
                                                const unsigned element_size)
{
    baryonyx::itm::shared_subvector s;
    s.reserve(size);

    Ensures(s.size() == 0);

    bool is_initialized = s.init(element_size);
    Ensures(is_initialized);

    unsigned element_total = size / element_size;

    for (unsigned int i = 0; i != element_total; ++i) {
        auto id = s.emplace();

        Ensures(id < size);
        for (auto [si, se] = s.element(id); si != se; ++si)
            s[si] = 99;
    }

    unsigned element_removed = 0;
    for (unsigned int i = 0; i != element_total / 2; ++i) {
        for (auto [si, se] = s.element(i * element_size); si != se; ++si) {
            Ensures(s[si] == 99);
            s[si] = 77;
        }

        s.remove(i * element_size);
        ++element_removed;
    }

    Ensures(element_removed == element_total / 2);

    for (unsigned int i = 0; i != element_removed; ++i) {
        auto id = s.emplace();

        Ensures(id < size);
        for (auto [si, se] = s.element(id); si != se; ++si)
            s[si] = 0;
    }

    for (unsigned int i = 0, e = s.size(); i != e; ++i)
        Ensures(s[i] == 99 || s[i] == 0);
}

void
test_size_type_greater_than_int8_subvector()
{
    test_size_type_greater_than_int8_subvector_impl(4096u, 16u);
    test_size_type_greater_than_int8_subvector_impl(8192u, 512u);
}

static void
check_knapsack_solver()
{
    struct rc
    {
        double value;
        int id;
    };

    struct ap
    {
        int value;
    };

    {
        // An example:
        // maximize: 16x1 + 19x2 + 23x3 + 26x4
        // st: 2x1 + 3x2 +4x3 + 5x4 <= 7

        auto R = std::make_unique<rc[]>(4);
        R[0] = { 16, 0 };
        R[1] = { 19, 1 };
        R[2] = { 23, 2 };
        R[3] = { 26, 3 };

        std::vector<ap> v{ { 0 }, { 1 }, { 2 }, { 3 } };

        std::vector<int> factors{ 2, 3, 4, 5 };
        int selected =
          baryonyx::itm::knapsack_dp_solver<baryonyx::itm::maximize_tag,
                                            double>(
            factors, R, v.begin(), 4, 7);

        Ensures(selected == 2);
        Ensures(R[0].id == 1 || R[1].id == 1);
        Ensures(R[0].id == 2 || R[1].id == 2);

        Ensures(
          std::accumulate(
            R.get(), R.get() + selected, 0.0, [](double init, const rc& elem) {
                return init + elem.value;
            }) == 42.0);
    }
}

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

static void
check_branch_and_bound_solver()
{
    struct rc
    {
        double value;
        int id;
        int f;
    };

    baryonyx::itm::branch_and_bound_solver<baryonyx::itm::maximize_tag, double>
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
        Ensures(result >= -1 && result <= static_cast<int>(R.size()));

        auto [r, f] = sum(R, nb, result);
        Ensures(f == 25);
        Ensures(r == 11.0);
    }

    {
        std::vector<rc> R = {
            { 16.0, 0, 2 }, { 19, 1, 3 }, { 23, 2, 4 }, { 26, 3, 5 }
        };

        const int nb = 4;
        const int b_min = 0;
        const int b_max = 7;

        auto result = bb.solve(R, nb, b_min, b_max);
        Ensures(result >= -1 && result <= static_cast<int>(R.size()));

        auto [r, f] = sum(R, nb, result);
        Ensures(f == 7);
        Ensures(r == 42.0);
    }

    {
        std::vector<rc> R = {
            { 1.0, 0, -1 }, { 1.0, 1, -1 }, { 1.0, 2, -1 }, { 1.0, 3, 1 }
        };

        const int nb = 4;
        const int b_min = -2;
        const int b_max = -2;

        auto result = bb.solve(R, nb, b_min, b_max);
        Ensures(result >= -1 && result <= static_cast<int>(R.size()));

        auto [r, f] = sum(R, nb, result);
        Ensures(f == -2);
        Ensures(r == 2.0); // Strange ?!?
    }
}

static void
check_exhaustive_solver()
{
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
        baryonyx::itm::exhaustive_solver<baryonyx::itm::maximize_tag, double>
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
        Ensures(result >= -1 && result <= static_cast<int>(R.size()));

        auto [r, f] = sum(R, nb, result);
        Ensures(f == 25);
        Ensures(r == 11.0);
    }

    {
        baryonyx::itm::exhaustive_solver<baryonyx::itm::maximize_tag, double>
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
        Ensures(result >= -1 && result <= static_cast<int>(R.size()));

        auto [r, f] = sum(R, nb, result);
        Ensures(f == 7);
        Ensures(r == 42.0);
    }

    {
        baryonyx::itm::exhaustive_solver<baryonyx::itm::maximize_tag, double>
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
        Ensures(result >= -1 && result <= static_cast<int>(R.size()));

        auto [r, f] = sum(R, nb, result);
        Ensures(f == -2);
        Ensures(r == 4.0);
    }

    {
        baryonyx::itm::exhaustive_solver<baryonyx::itm::maximize_tag, double>
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
            Ensures(result >= -1 && result <= static_cast<int>(R1.size()));

            auto [r, f] = sum(R1, 4, result);
            Ensures(f == -2);
            Ensures(r == 4.0);
        }

        {
            auto result = ex.solve(1, R2, 2);
            Ensures(result >= -1 && result <= static_cast<int>(R2.size()));

            auto [r, f] = sum(R2, 2, result);
            Ensures(f == -1);
            Ensures(r == 1.0);
        }
    }
}

void
check_observer_pnm()
{
    {
        std::vector<double> v{ 0,   -1, -3,  -2, 0,    3, 3,
                               1.5, 2,  0.5, 2,  -1.5, 1, 4 };

        baryonyx::pnm_array obs(7, 2);
        Ensures(obs);

        std::transform(
          v.begin(), v.end(), obs.begin(), baryonyx::colormap(-5.f, 5.f));

        for (auto ob : obs)
            fmt::print("{}/{}/{} ", ob.red(), ob.green(), ob.blue());

        obs("test.pnm");
    }

    {
        baryonyx::pnm_vector obs("test2.pnm", 4, 4);
        Ensures(obs);

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
}

void
check_show_size()
{
    double size;
    baryonyx::show_size_type type;

    std::tie(size, type) = baryonyx::memory_consumed_size(0);
    Ensures(size == 0 && type == baryonyx::show_size_type::B);

    std::tie(size, type) = baryonyx::memory_consumed_size(1024);
    Ensures(size == 1 && type == baryonyx::show_size_type::KB);

    std::tie(size, type) = baryonyx::memory_consumed_size(2048);
    Ensures(size == 2 && type == baryonyx::show_size_type::KB);

    std::tie(size, type) = baryonyx::memory_consumed_size(2048);
    Ensures(size == 2 && type == baryonyx::show_size_type::KB);

    std::tie(size, type) = baryonyx::memory_consumed_size(1049088);
    Ensures(size > 1 && size < 2 && type == baryonyx::show_size_type::MB);

    std::tie(size, type) = baryonyx::memory_consumed_size(1049088);
    Ensures(size > 1 && size < 2 && type == baryonyx::show_size_type::MB);

    std::tie(size, type) = baryonyx::memory_consumed_size(17179869696);
    Ensures(size > 16 && size < 17 && type == baryonyx::show_size_type::GB);
}

void
check_trim_functions()
{
    {
        std::string str = "    xxxx    ";
        Ensures(baryonyx::left_trim(str) == "xxxx    ");
        Ensures(baryonyx::right_trim(str) == "    xxxx");
        Ensures(baryonyx::trim(str) == "xxxx");
        Ensures(str == "    xxxx    ");
    }

    {
        std::string str = "xxxx    ";
        Ensures(baryonyx::left_trim(str) == "xxxx    ");
        Ensures(baryonyx::right_trim(str) == "xxxx");
        Ensures(baryonyx::trim(str) == "xxxx");
        Ensures(str == "xxxx    ");
    }

    {
        std::string str = "    xxxx";
        Ensures(baryonyx::left_trim(str) == "xxxx");
        Ensures(baryonyx::right_trim(str) == "    xxxx");
        Ensures(baryonyx::trim(str) == "xxxx");
        Ensures(str == "    xxxx");
    }

    {
        std::string str = "    ";
        Ensures(baryonyx::left_trim(str) == "");
        Ensures(baryonyx::right_trim(str) == "");
        Ensures(baryonyx::trim(str) == "");
        Ensures(str == "    ");
    }

    {
        std::string str = "";
        Ensures(baryonyx::left_trim(str) == "");
        Ensures(baryonyx::right_trim(str) == "");
        Ensures(baryonyx::trim(str) == "");
        Ensures(str == "");
    }

    {
        std::string str = "abc";
        Ensures(baryonyx::left_trim(str) == "abc");
        Ensures(baryonyx::right_trim(str) == "abc");
        Ensures(baryonyx::trim(str) == "abc");
        Ensures(str == "abc");
    }
}

void
check_bit_array()
{
    baryonyx::bit_array a(100);

    Ensures(baryonyx::bit_array::block_t == 4 * 8);
    Ensures(a.size() == 100);
    Ensures(a.block_size() == 1 + (100 / (4 * 8)));
    Ensures(a.block_size() == 4);

    a.zeros();
    for (auto i = 0; i != a.size(); ++i)
        Ensures(a.get(i) == 0);

    a.ones();
    for (auto i = 0; i != a.size(); ++i)
        Ensures(a.get(i) == 1);

    for (auto i = 0; i != a.size(); ++i)
        a.unset(i);
    for (auto i = 0; i != a.size(); ++i)
        Ensures(a.get(i) == 0);

    for (auto i = 0; i != a.size(); ++i)
        a.set(i);
    for (auto i = 0; i != a.size(); ++i)
        Ensures(a.get(i) == 1);

    for (auto i = 0; i != a.size(); i += 2)
        a.unset(i);
    for (auto i = 0; i != a.size(); i += 2)
        Ensures(a.get(i) == 0 && a.get(i + 1) == 1);

    for (auto i = 0; i != a.size(); ++i)
        a.invert(i);

    for (auto i = 0; i != a.size(); i += 2)
        Ensures(a.get(i) == 1 && a.get(i + 1) == 0);
}

int
main(int /* argc */, char* /* argv */ [])
{
    unit_test::checks("check_print_api", check_print_api);
    unit_test::checks("check_numeric_cast", check_numeric_cast);
    unit_test::checks("check_fixed_array", check_fixed_array);
    unit_test::checks("check_fixed_2darray", check_fixed_2darray);
    unit_test::checks("check_knapsack_solver", check_knapsack_solver);

    unit_test::checks("branch_and_bound::subvector api", test_subvector_api);
    unit_test::checks("branch_and_bound::subvector size",
                      test_size_type_greater_than_int8_subvector);
    unit_test::checks("branch_and_bound::subvector size",
                      test_size_type_greater_than_int8_subvector);
    unit_test::checks("check_branch_and_bound_solver",
                      check_branch_and_bound_solver);
    unit_test::checks("check_exhaustive_solver", check_exhaustive_solver);

    unit_test::checks("check_observer_pnm", check_observer_pnm);
    unit_test::checks("check_show_size", check_show_size);
    unit_test::checks("check_trim_functions", check_trim_functions);
    unit_test::checks("check_bit_array", check_bit_array);

    return unit_test::report_errors();
}
