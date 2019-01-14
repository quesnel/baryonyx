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

#ifndef ORG_VLEPROJECT_UNIT_TEST_HPP
#define ORG_VLEPROJECT_UNIT_TEST_HPP

#include <fmt/color.h>
#include <fmt/format.h>

namespace unit_test {
namespace detail {

struct unit_test_error
{
    const char* what() noexcept
    {
        return "unit test error";
    }
};

struct tester
{
    int errors = 0;
    bool called_report_function = false;

    tester& operator++() noexcept
    {
        errors++;

        return *this;
    }

    ~tester() noexcept
    {
        if (called_report_function == false) {
            fmt::print(fmt::color::red,
                       "unit_test::report_errors() not called.\n\nUsage:\n"
                       "int main(int argc, char*[] argc)\n"
                       "{\n"
                       "    [...]\n"
                       "    return unit_test::report_errors();\n"
                       "}\n");

            std::abort();
        }
    }
};

inline tester&
test_errors()
{
    static tester t;

    return t;
}

inline void
ensures_impl(const char* expr, const char* file, int line, const char* function)
{
    fmt::print(fmt::color::red,
               "{} ({}): test '{}' failed in function '{}'\n",
               file,
               line,
               expr,
               function);

    ++test_errors();
}

inline void
ensures_equal_impl(const char* expr1,
                   const char* expr2,
                   const char* file,
                   int line,
                   const char* function)
{
    fmt::print(fmt::color::red,
               "{} ({}): test '{} == {}' failed in function '{}'\n",
               file,
               line,
               expr1,
               expr2,
               function);

    ++test_errors();
}

inline void
ensures_not_equal_impl(const char* expr1,
                       const char* expr2,
                       const char* file,
                       int line,
                       const char* function)
{
    fmt::print(fmt::color::red,
               "{} ({}): test '{} != {}' failed in function '{}'\n",
               file,
               line,
               expr1,
               expr2,
               function);

    ++test_errors();
}

inline void
ensures_throw_impl(const char* excep,
                   const char* file,
                   int line,
                   const char* function)
{
    fmt::print(fmt::color::red,
               "{} ({}): exception '{}' throw failed in function '{}'\n",
               file,
               line,
               excep,
               function);

    ++test_errors();
}

inline void
ensures_not_throw_impl(const char* excep,
                       const char* file,
                       int line,
                       const char* function)
{
    fmt::print(fmt::color::red,
               "{} ({}): exception '{}' not throw failed in function '{}'\n",
               file,
               line,
               excep,
               function);

    ++test_errors();
}

} // namespace details

inline int
report_errors()
{
    auto& tester = unit_test::detail::test_errors();
    tester.called_report_function = true;
    int errors = tester.errors;

    if (errors == 0)
        fmt::print(fmt::color::green, "\nNo errors detected.\n");
    else
        fmt::print(fmt::color::orange_red,
                   "\n{} error{} detected.\n",
                   errors,
                   (errors == 1 ? "" : "s"));

    return errors;
}

template<typename F>
void
checks(const char* str, F f) noexcept
{
    try {
        f();
    } catch (const detail::unit_test_error&) {
    } catch (const std::exception& e) {
        fmt::print(fmt::color::red, "{}: Exception throw: {}\n", str, e.what());
        ++detail::test_errors();
    } catch (...) {
        fmt::print(fmt::color::red, "{}: Unknown exception throw\n", str);
        ++detail::test_errors();
    }
}

} // namespace unit_test

#define Ensures(expr)                                                          \
    do {                                                                       \
        if (!(expr)) {                                                         \
            unit_test::detail::ensures_impl(                                   \
              #expr, __FILE__, __LINE__, __func__);                            \
            throw unit_test::detail::unit_test_error();                        \
        }                                                                      \
    } while (0)

#define EnsuresEqual(expr1, expr2)                                             \
    do {                                                                       \
        if (!((expr1) == (expr2))) {                                           \
            unit_test::detail::ensures_equal_impl(                             \
              #expr1, #expr2, __FILE__, __LINE__, __func__);                   \
            throw unit_test::detail::unit_test_error();                        \
        }                                                                      \
    } while (0)

#define EnsuresNotEqual(expr1, expr2)                                          \
    do {                                                                       \
        if (!((expr1) != (expr2))) {                                           \
            unit_test::detail::ensures_not_equal_impl(                         \
              #expr1, #expr2, __FILE__, __LINE__, __func__);                   \
            throw unit_test::detail::unit_test_error();                        \
        }                                                                      \
    } while (0)

#define EnsuresThrow(expr, Excep)                                              \
    do {                                                                       \
        try {                                                                  \
            expr;                                                              \
            unit_test::detail::ensures_throw_impl(                             \
              #Excep, __FILE__, __LINE__, __func__);                           \
            throw unit_test::detail::unit_test_error();                        \
        } catch (const Excep&) {                                               \
        } catch (...) {                                                        \
            unit_test::detail::ensures_throw_impl(                             \
              #Excep, __FILE__, __LINE__, __func__);                           \
            throw unit_test::detail::unit_test_error();                        \
        }                                                                      \
    } while (0)

#define EnsuresNotThrow(expr, Excep)                                           \
    do {                                                                       \
        try {                                                                  \
            expr;                                                              \
        } catch (const Excep&) {                                               \
            unit_test::detail::ensures_not_throw_impl(                         \
              #Excep, __FILE__, __LINE__, __func__);                           \
            throw unit_test::detail::unit_test_error();                        \
        } catch (...) {                                                        \
        }                                                                      \
    } while (0)

#endif // #ifndef ORG_VLEPROJECT_UNIT_TEST_HPP
