/* Copyright (C) 2016-2018 INRA
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

#include <fmt/format.h>

#ifndef _WIN32
#include <unistd.h>
#endif

#include <cstdlib>

namespace unit_test {
namespace detail {

#ifndef _WIN32
inline bool
have_color() noexcept
{
    return ::isatty(::fileno(stderr));
}
#else
inline bool
have_color() noexcept
{
    return false;
}
#endif

#define COLOR_RESET "\033[0m"
#define BOLD "\033[1m"
#define BLACK_TEXT "\033[30;1m"
#define RED_TEXT "\033[31;1m"
#define GREEN_TEXT "\033[32;1m"
#define YELLOW_TEXT "\033[33;1m"
#define BLUE_TEXT "\033[34;1m"
#define MAGENTA_TEXT "\033[35;1m"
#define CYAN_TEXT "\033[36;1m"
#define WHITE_TEXT "\033[37;1m"

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
            if (have_color()) {
                fmt::print(stderr,
                           RED_TEXT
                           "unit_test::report_errors() not called.\n\nUsage:\n"
                           "int main(int argc, char*[] argc)\n"
                           "{\n"
                           "    [...]\n"
                           "    return unit_test::report_errors();\n"
                           "}\n" COLOR_RESET);
            } else {
                fmt::print(stderr,
                           "unit_test::report_errors() not called.\n\nUsage:\n"
                           "int main(int argc, char*[] argc)\n"
                           "{\n"
                           "    [...]\n"
                           "    return unit_test::report_errors();\n"
                           "}\n");
            }

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
ensures_impl(const char* expr,
             const char* file,
             int line,
             const char* function)
{
    if (have_color()) {
        fmt::print(stderr,
                   RED_TEXT
                   "{} ({}): test '{}' failed in function '{}'\n" COLOR_RESET,
                   file,
                   line,
                   expr,
                   function);
    } else {
        fmt::print(stderr,
                   "{} ({}): test '{}' failed in function '{}'\n",
                   file,
                   line,
                   expr,
                   function);
    }

    ++test_errors();
}

inline void
ensures_equal_impl(const char* expr1,
                   const char* expr2,
                   const char* file,
                   int line,
                   const char* function)
{
    if (have_color()) {
        fmt::print(
          stderr,
          RED_TEXT
          "{} ({}): test '{} == {}' failed in function '{}'\n" COLOR_RESET,
          file,
          line,
          expr1,
          expr2,
          function);
    } else {
        fmt::print(stderr,
                   "{} ({}): test '{} == {}' failed in function '{}'\n",
                   file,
                   line,
                   expr1,
                   expr2,
                   function);
    }

    ++test_errors();
}

inline void
ensures_not_equal_impl(const char* expr1,
                       const char* expr2,
                       const char* file,
                       int line,
                       const char* function)
{
    if (have_color()) {
        fmt::print(
          stderr,
          RED_TEXT
          "{} ({}): test '{} != {}' failed in function '{}'\n" COLOR_RESET,
          file,
          line,
          expr1,
          expr2,
          function);
    } else {
        fmt::print(stderr,
                   "{} ({}): test '{} != {}' failed in function '{}'\n",
                   file,
                   line,
                   expr1,
                   expr2,
                   function);
    }

    ++test_errors();
}

inline void
ensures_throw_impl(const char* excep,
                   const char* file,
                   int line,
                   const char* function)
{
    if (have_color()) {
        fmt::print(stderr,
                   RED_TEXT "{} ({}): exception '{}' throw failed in function "
                            "'{}'\n" COLOR_RESET,
                   file,
                   line,
                   excep,
                   function);
    } else {
        fmt::print(stderr,
                   "{} ({}): exception '{}' throw failed in function '{}'\n",
                   file,
                   line,
                   excep,
                   function);
    }

    ++test_errors();
}

inline void
ensures_not_throw_impl(const char* excep,
                       const char* file,
                       int line,
                       const char* function)
{
    if (have_color()) {
        fmt::print(stderr,
                   RED_TEXT "{} ({}): exception '{}' not throw failed in "
                            "function '{}'\n" COLOR_RESET,
                   file,
                   line,
                   excep,
                   function);
    } else {
        fmt::print(
          stderr,
          "{} ({}): exception '{}' not throw failed in function '{}'\n",
          file,
          line,
          excep,
          function);
    }

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
        fmt::print(stderr, "No errors detected.\n");
    else {
        if (unit_test::detail::have_color()) {
            fmt::print(stderr,
                       RED_TEXT "{} error{} detected.\n" COLOR_RESET,
                       errors,
                       (errors == 1 ? "" : "s"));
        } else {
            fmt::print(stderr,
                       "{} error{} detected.\n",
                       errors,
                       (errors == 1 ? "" : "s"));
        }
    }

    return errors;
}

} // namespace unit_test

#define Ensures(expr)                                                         \
    do {                                                                      \
        if (!(expr)) {                                                        \
            unit_test::detail::ensures_impl(                                  \
              #expr, __FILE__, __LINE__, __func__);                           \
            return;                                                           \
        }                                                                     \
    } while (0)

#define EnsuresEqual(expr1, expr2)                                            \
    do {                                                                      \
        if (!((expr1) == (expr2))) {                                          \
            unit_test::detail::ensures_equal_impl(                            \
              #expr1, #expr2, __FILE__, __LINE__, __func__);                  \
            return;                                                           \
        }                                                                     \
    } while (0)

#define EnsuresNotEqual(expr1, expr2)                                         \
    do {                                                                      \
        if (!((expr1) != (expr2))) {                                          \
            unit_test::detail::ensures_not_equal_impl(                        \
              #expr1, #expr2, __FILE__, __LINE__, __func__);                  \
            return;                                                           \
        }                                                                     \
    } while (0)

#define EnsuresThrow(expr, Excep)                                             \
    do {                                                                      \
        try {                                                                 \
            expr;                                                             \
            unit_test::detail::ensures_throw_impl(                            \
              #Excep, __FILE__, __LINE__, __func__);                          \
            return;                                                           \
        } catch (const Excep&) {                                              \
        } catch (...) {                                                       \
            unit_test::detail::ensures_throw_impl(                            \
              #Excep, __FILE__, __LINE__, __func__);                          \
            return;                                                           \
        }                                                                     \
    } while (0)

#define EnsuresNotThrow(expr, Excep)                                          \
    do {                                                                      \
        try {                                                                 \
            expr;                                                             \
        } catch (const Excep&) {                                              \
            unit_test::detail::ensures_not_throw_impl(                        \
              #Excep, __FILE__, __LINE__, __func__);                          \
            return;                                                           \
        } catch (...) {                                                       \
        }                                                                     \
    } while (0)

#endif // #ifndef ORG_VLEPROJECT_UNIT_TEST_HPP
