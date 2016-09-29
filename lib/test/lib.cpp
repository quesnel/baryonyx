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

#include <lpcore>

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

void check_numeric_cast()
{
    int small_positive = 1;
    int small_negative = -1;
    int large_positive = std::numeric_limits<int>::max();
    int large_negative = std::numeric_limits<int>::min();

    assert(lp::is_numeric_castable<signed char>(small_positive));
    assert(lp::is_numeric_castable<signed char>(small_negative));
    assert(not lp::is_numeric_castable<signed char>(large_positive));
    assert(not lp::is_numeric_castable<signed char>(large_negative));

    assert(lp::is_numeric_castable<unsigned char>(small_positive));
    assert(not lp::is_numeric_castable<unsigned char>(small_negative));
    assert(not lp::is_numeric_castable<unsigned char>(large_positive));
    assert(not lp::is_numeric_castable<unsigned char>(large_negative));

    assert(lp::is_numeric_castable<signed int>(small_positive));
    assert(lp::is_numeric_castable<signed int>(small_negative));
    assert(lp::is_numeric_castable<signed int>(large_positive));
    assert(lp::is_numeric_castable<signed int>(large_negative));

    assert(lp::is_numeric_castable<unsigned int>(small_positive));
    assert(not lp::is_numeric_castable<unsigned int>(small_negative));
    assert(lp::is_numeric_castable<unsigned int>(large_positive));
    assert(not lp::is_numeric_castable<unsigned int>(large_negative));

    assert(lp::is_numeric_castable<long long>(small_positive));
    assert(lp::is_numeric_castable<long long>(large_negative));
    assert(lp::is_numeric_castable<long long>(small_positive));
    assert(lp::is_numeric_castable<long long>(large_negative));

    assert(lp::is_numeric_castable<unsigned long long>(small_positive));
    assert(not lp::is_numeric_castable<unsigned long long>(small_negative));
    assert(lp::is_numeric_castable<unsigned long long>(large_positive));
    assert(not lp::is_numeric_castable<unsigned long long>(large_negative));

    assert(not lp::is_numeric_castable<size_t>(small_negative));
    assert(not lp::is_numeric_castable<size_t>(large_negative));

    try {
        std::vector<int> v;
        unsigned int checked_size = lp::numeric_cast<unsigned int>(v.size());
        assert(0 == checked_size);
    } catch(const std::exception& e) {
        printf("%s\n", e.what());
        assert(false && "bad cast");
    }
}

void check_parameter()
{
    lp::parameter real {3.0};
    assert(real.type == lp::parameter::tag::real);

    lp::parameter integer {1000l};
    assert(integer.type == lp::parameter::tag::integer);

    lp::parameter str {"hello world"};
    assert(str.type == lp::parameter::tag::string);

    str = real;
    assert(str.type == lp::parameter::tag::real);
    assert(str.d == 3.0);

    str = integer;
    assert(str.type == lp::parameter::tag::integer);
    assert(str.l == 1000l);

    std::vector<lp::parameter> x(100);
    for (auto& elem : x) {
        assert(elem.type == lp::parameter::tag::integer);
        assert(elem.l == 0l);
    }

    auto y = lp::parameter(4.0);
    assert(y.type == lp::parameter::tag::real);
    assert(y.d == 4.0);

    x[0] = lp::parameter(5.0);
    assert(x[0].type == lp::parameter::tag::real);
    assert(x[0].d == 5.0);

    x[0].swap(x[1]);
    assert(x[0].type == lp::parameter::tag::integer);
    assert(x[0].l == 0l);
    assert(x[1].type == lp::parameter::tag::real);
    assert(x[1].d == 5.0);

    x[2] = std::move(x[1]);
    assert(x[0].type == lp::parameter::tag::integer);
    assert(x[0].l == 0l);
    assert(x[1].type == lp::parameter::tag::integer);
    assert(x[1].l == 0l);
    assert(x[2].type == lp::parameter::tag::real);
    assert(x[2].d == 5.0);
}

int main(int /* argc */, char */* argv */[])
{
    check_numeric_cast();
    check_parameter();
}
