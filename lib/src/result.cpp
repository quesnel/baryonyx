/* Copyright (C) 2016-2026 INRAE
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

#include "result.hpp"
#include "private.hpp"
#include "utils.hpp"

#include <baryonyx/core-utils>

#include <charconv>
#include <fstream>
#include <string>
#include <string_view>
#include <system_error>
#include <unordered_set>
#include <vector>

namespace baryonyx {

result
make_result(const baryonyx::context_ptr& ctx, const std::string& filename)
{
    if (ctx) {
        info(*ctx, "solution reads from file {}\n", filename);
        std::ifstream ifs(filename);
        if (!ifs.is_open())
            return result{ result_status::uninitialized };

        result ret;
        ifs >> ret;
        return ret;
    }

    return result{ result_status::empty_context };
}

namespace {

[[noreturn]] void
fail(file_format_error_tag tag, int line, std::string_view text)
{
    throw file_format_failure(std::string(text), tag, line, 0);
}

} // namespace

std::istream&
operator>>(std::istream& is, result& ret)
{
    ret = result{};
    ret.strings = std::make_shared<baryonyx::string_buffer>();
    ret.status = result_status::success;

    std::unordered_set<std::string_view> names;
    std::vector<std::int8_t> values;
    bool no_solution = false;
    int line = 0;

    for (std::string buffer; std::getline(is, buffer);) {
        ++line;

        const auto text = trim(buffer);
        if (text.empty())
            continue;

        if (text.front() == '\\') {
            // Written in place of the assignments when no solution was found.
            if (text.find("Solution not found") != std::string_view::npos)
                no_solution = true;
            continue;
        }

        const auto equal = text.find('=');
        if (equal == std::string_view::npos)
            fail(file_format_error_tag::bad_value, line, text);

        const auto name = trim(text.substr(0, equal));
        const auto value = trim(text.substr(equal + 1));

        if (name.empty())
            fail(file_format_error_tag::bad_name, line, text);

        int number = -1;
        const auto [end, ec] =
          std::from_chars(value.data(), value.data() + value.size(), number);
        if (ec != std::errc{} || end != value.data() + value.size() ||
            (number != 0 && number != 1))
            fail(file_format_error_tag::bad_value, line, text);

        const auto stored = ret.strings->append(name);
        if (!names.insert(stored).second)
            fail(file_format_error_tag::bad_name, line, text);

        ret.variable_name.emplace_back(stored);
        values.emplace_back(static_cast<std::int8_t>(number));
    }

    if (is.bad())
        fail(file_format_error_tag::bad_end_of_file, line, "read error");

    if (no_solution) {
        if (!values.empty())
            fail(file_format_error_tag::bad_end_of_file,
                 line,
                 "assignments in a file without solution");

        ret.variable_name.clear();
        ret.status = result_status::time_limit_reached;
        return is;
    }

    // An empty or truncated file, for example when the run was killed.
    if (values.empty())
        fail(file_format_error_tag::bad_end_of_file, line, "no assignment");

    ret.solutions.emplace_back(std::move(values), 0.0);

    return is;
}
}
