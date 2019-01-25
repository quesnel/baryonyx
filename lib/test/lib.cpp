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

#include "branch-and-bound-solver.hpp"
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
check_clamp()
{
    Ensures(baryonyx::clamp(0.0, 0.0, 1.0) == 0.0);
    Ensures(baryonyx::clamp(1.0, 0.0, 1.0) == 1.0);
    Ensures(baryonyx::clamp(-0.5, 0.0, 1.0) == 0.0);
    Ensures(baryonyx::clamp(1.5, 0.0, 1.0) == 1.0);
    Ensures(baryonyx::clamp(168, -128, +127) == 127);
    Ensures(baryonyx::clamp(168, 0, +255) == 168);
    Ensures(baryonyx::clamp(128, -128, +127) == 127);
    Ensures(baryonyx::clamp(128, 0, +255) == 128);

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

struct rc
{
    double value;
    int id;
};

struct ap
{
    int value;
};

static void
check_knapsack_solver()
{
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

static void
check_branch_and_bound_solver()
{
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
          baryonyx::itm::branch_and_bound_solver<baryonyx::itm::maximize_tag,
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

    // {
    //     std::vector<rc> R{ { 15, 0 }, { 19, 1 }, { 13, 2 }, { 12, 3 } };
    //     std::vector<int> factors{ 2, 1, 3, 2 };

    //     int selected =
    //       baryonyx::itm::branch_and_bound_solver<baryonyx::itm::minimize_tag,
    //                                              double>(
    //         R, factors.begin(), factors.end(), 3);

    //     Ensures(selected == 1);
    //     Ensures(R[0].id == 2 || R[1].id == 2);
    //     Ensures(std::accumulate(R.begin(),
    //                             R.begin() + selected,
    //                             0.0,
    //                             [](double init, const rc& elem) {
    //                                 return init + elem.value;
    //                             }) == 13.0);
    // }

    // {
    //     std::vector<rc> R{
    //         { 16, 0 }, { 19, 1 }, { 23, 2 }, { 28, 3 }, { 5, 4 }
    //     };
    //     std::vector<int> factors{ 2, 3, 4, 5, 7 };

    //     int selected =
    //       baryonyx::itm::branch_and_bound_solver<baryonyx::itm::maximize_tag,
    //                                              double>(
    //         R, factors.begin(), factors.end(), 7);

    //     Ensures(selected == 2);
    //     Ensures(std::accumulate(R.begin(),
    //                             R.begin() + selected,
    //                             0.0,
    //                             [](double init, const rc& elem) {
    //                                 return init + elem.value;
    //                             }) == 44.0);
    // }

    // {
    //     std::vector<rc> R{ { 1000, 0 }, { 16, 1 }, { 19, 2 },
    //                        { 23, 3 },   { 28, 4 }, { 1, 5 } };
    //     std::vector<int> factors{ 1, 2, 3, 4, 5, 6 };

    //     auto selected =
    //       baryonyx::itm::branch_and_bound_solver<baryonyx::itm::maximize_tag,
    //                                              double>(
    //         R, factors.begin(), factors.end(), 7);

    //     Ensures(selected == 3);
    //     Ensures(R[0].id == 0 || R[1].id == 0 || R[2].id == 0);
    //     Ensures(R[0].id == 1 || R[1].id == 1 || R[2].id == 1);
    //     Ensures(R[0].id == 3 || R[1].id == 3 || R[2].id == 3);
    // }

    // {
    //     std::vector<rc> R{ { 1000, 0 }, { 16, 1 }, { 19, 2 },
    //                        { 23, 3 },   { 28, 4 }, { 1, 5 } };
    //     std::vector<int> factors{ 1, 2, 3, 4, 5, 6 };

    //     auto selected =
    //       baryonyx::itm::branch_and_bound_solver<baryonyx::itm::minimize_tag,
    //                                              double>(
    //         R, factors.begin(), factors.end(), 7);

    //     Ensures(selected == 1);
    //     Ensures(R[0].id == 5);
    // }

    // {
    //     std::vector<rc> R{ { 1e-7, 0 }, { 1e-9, 1 }, { .5e-4, 2 },
    //                        { 1e-7, 3 }, { 1e-3, 4 }, { 1e-9, 5 } };
    //     std::vector<int> factors{ 1, 1, 6, 1, 1, 1 };

    //     auto selected =
    //       baryonyx::itm::branch_and_bound_solver<baryonyx::itm::minimize_tag,
    //                                              double>(
    //         R, factors.begin(), factors.end(), 7);

    //     Ensures(selected == 1);
    // }
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

int
main(int /* argc */, char* /* argv */ [])
{
    unit_test::checks("check_clamp", check_clamp);
    unit_test::checks("check_numeric_cast", check_numeric_cast);
    unit_test::checks("check_fixed_array", check_fixed_array);
    unit_test::checks("check_fixed_2darray", check_fixed_2darray);
    unit_test::checks("check_knapsack_solver", check_knapsack_solver);
    unit_test::checks("check_branch_and_bound_solver",
                      check_branch_and_bound_solver);
    unit_test::checks("check_observer_pnm", check_observer_pnm);
    unit_test::checks("check_show_size", check_show_size);

    return unit_test::report_errors();
}
