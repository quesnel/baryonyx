/* Copyright (C) 2016 INRA
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

#include "matrix.hpp"
#include "scoped_array.hpp"
#include "unit-test.hpp"

#include <lpcore>

void
check_numeric_cast()
{
    int small_positive = 1;
    int small_negative = -1;
    int large_positive = std::numeric_limits<int>::max();
    int large_negative = std::numeric_limits<int>::min();

    Ensures(lp::is_numeric_castable<signed char>(small_positive));
    Ensures(lp::is_numeric_castable<signed char>(small_negative));
    Ensures(not lp::is_numeric_castable<signed char>(large_positive));
    Ensures(not lp::is_numeric_castable<signed char>(large_negative));

    Ensures(lp::is_numeric_castable<unsigned char>(small_positive));
    Ensures(not lp::is_numeric_castable<unsigned char>(small_negative));
    Ensures(not lp::is_numeric_castable<unsigned char>(large_positive));
    Ensures(not lp::is_numeric_castable<unsigned char>(large_negative));

    Ensures(lp::is_numeric_castable<signed int>(small_positive));
    Ensures(lp::is_numeric_castable<signed int>(small_negative));
    Ensures(lp::is_numeric_castable<signed int>(large_positive));
    Ensures(lp::is_numeric_castable<signed int>(large_negative));

    Ensures(lp::is_numeric_castable<unsigned int>(small_positive));
    Ensures(not lp::is_numeric_castable<unsigned int>(small_negative));
    Ensures(lp::is_numeric_castable<unsigned int>(large_positive));
    Ensures(not lp::is_numeric_castable<unsigned int>(large_negative));

    Ensures(lp::is_numeric_castable<long long>(small_positive));
    Ensures(lp::is_numeric_castable<long long>(large_negative));
    Ensures(lp::is_numeric_castable<long long>(small_positive));
    Ensures(lp::is_numeric_castable<long long>(large_negative));

    Ensures(lp::is_numeric_castable<unsigned long long>(small_positive));
    Ensures(not lp::is_numeric_castable<unsigned long long>(small_negative));
    Ensures(lp::is_numeric_castable<unsigned long long>(large_positive));
    Ensures(not lp::is_numeric_castable<unsigned long long>(large_negative));

    Ensures(not lp::is_numeric_castable<size_t>(small_negative));
    Ensures(not lp::is_numeric_castable<size_t>(large_negative));

    std::vector<int> v;

    EnsuresNotThrow(lp::numeric_cast<short int>(v.size()), std::exception);

    EnsuresNotThrow(lp::numeric_cast<short int>(v.capacity()), std::exception);

    EnsuresThrow(lp::numeric_cast<short int>(v.max_size()), std::exception);

    unsigned int checked_size = lp::numeric_cast<unsigned int>(v.size());
    Ensures(0 == checked_size);
}

void
check_parameter()
{
    lp::parameter real{ 3.0 };
    Ensures(real.type == lp::parameter::tag::real);

    lp::parameter integer{ 1000l };
    Ensures(integer.type == lp::parameter::tag::integer);

    lp::parameter str{ "hello world" };
    Ensures(str.type == lp::parameter::tag::string);

    str = real;
    Ensures(str.type == lp::parameter::tag::real);
    Ensures(str.d == 3.0);

    str = integer;
    Ensures(str.type == lp::parameter::tag::integer);
    Ensures(str.l == 1000l);

    std::vector<lp::parameter> x(100);
    for (auto& elem : x) {
        Ensures(elem.type == lp::parameter::tag::integer);
        Ensures(elem.l == 0l);
    }

    auto y = lp::parameter(4.0);
    Ensures(y.type == lp::parameter::tag::real);
    Ensures(y.d == 4.0);

    x[0] = lp::parameter(5.0);
    Ensures(x[0].type == lp::parameter::tag::real);
    Ensures(x[0].d == 5.0);

    x[0].swap(x[1]);
    Ensures(x[0].type == lp::parameter::tag::integer);
    Ensures(x[0].l == 0l);
    Ensures(x[1].type == lp::parameter::tag::real);
    Ensures(x[1].d == 5.0);

    x[2] = std::move(x[1]);
    Ensures(x[0].type == lp::parameter::tag::integer);
    Ensures(x[0].l == 0l);
    Ensures(x[1].type == lp::parameter::tag::integer);
    Ensures(x[1].l == 0l);
    Ensures(x[2].type == lp::parameter::tag::real);
    EnsuresEqual(x[2].d, 5.0);
    EnsuresNotEqual(x[2].d, 6.0);

    x[3] = lp::parameter(std::string("hello world!"));
    Ensures(x[3].type == lp::parameter::tag::string);
    Ensures(x[3].s == "hello world!");
}

void
check_matrix()
{
    lp::SparseArray<int, double> m(4, 2);

    EnsuresThrow(m.P(0, 0), std::out_of_range);
    EnsuresThrow(m.P(0, 1), std::out_of_range);
    EnsuresThrow(m.P(1, 0), std::out_of_range);
    EnsuresThrow(m.P(1, 1), std::out_of_range);
    EnsuresThrow(m.P(2, 0), std::out_of_range);
    EnsuresThrow(m.P(2, 1), std::out_of_range);
    EnsuresThrow(m.P(3, 0), std::out_of_range);
    EnsuresThrow(m.P(3, 1), std::out_of_range);
    Ensures(m.size() == 0);

    m.set(1, 0, 1, 1.0);
    m.set(0, 1, 2, 2.0);
    m.set(3, 1, 3, 3.0);
    m.set(2, 1, 4, 4.0);
    m.sort();

    Ensures(m.size() == 4);

    EnsuresThrow(m.P(0, 0), std::out_of_range);
    Ensures(m.A(0, 1) == 2);
    Ensures(m.P(0, 1) == 2.0);

    Ensures(m.A(1, 0) == 1);
    Ensures(m.P(1, 0) == 1.0);
    EnsuresThrow(m.A(1, 1), std::out_of_range);
    EnsuresThrow(m.A(2, 0), std::out_of_range);
    Ensures(m.A(2, 1) == 4);
    Ensures(m.P(2, 1) == 4.0);
    EnsuresThrow(m.A(3, 0), std::out_of_range);
    Ensures(m.A(3, 1) == 3);
    Ensures(m.P(3, 1) == 3.0);
    Ensures(m.size() == 4);

    Ensures(m.rows() == 4);
    Ensures(m.columns() == 2);

    Ensures(m.row(0).size() == 1);
    Ensures(m.row(1).size() == 1);
    Ensures(m.row(2).size() == 1);
    Ensures(m.row(3).size() == 1);
    Ensures(m.column(0).size() == 1);
    Ensures(m.column(1).size() == 3);

    Ensures(m.A().size() == 4);
    Ensures(m.P()[0] == 1.0);
    Ensures(m.P()[1] == 2.0);
    Ensures(m.P()[2] == 3.0);
    Ensures(m.P()[3] == 4.0);
}

void
check_scoped_array()
{
    auto* array = new int[10];

    lp::scoped_array<int> a(array);

    Ensures(a.get() == array);

    array[0] = 1;
    array[1] = 2;
    array[2] = 3;
    array[3] = 4;
    array[4] = 5;
    array[5] = 6;
    array[6] = 7;
    array[7] = 8;
    array[8] = 9;
    array[9] = 0;

    Ensures(a[0] == 1);
    Ensures(a[1] == 2);
    Ensures(a[2] == 3);
    Ensures(a[3] == 4);
    Ensures(a[4] == 5);
    Ensures(a[5] == 6);
    Ensures(a[6] == 7);
    Ensures(a[7] == 8);
    Ensures(a[8] == 9);
    Ensures(a[9] == 0);

    lp::scoped_array<int> b;

    b = std::move(a);

    Ensures(b.get() == array);
    Ensures(a.get() == nullptr);
}

int
main(int /* argc */, char* /* argv */ [])
{
    check_numeric_cast();
    check_parameter();
    check_matrix();
    check_scoped_array();

    return unit_test::report_errors();
}
