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

#include "result.hpp"
#include "private.hpp"

#include <baryonyx/core-utils>

#include <fstream>
#include <optional>
#include <string>

#include <cctype>

namespace baryonyx {

result
make_result(const baryonyx::context_ptr& ctx, const std::string& filename)
{
    if (ctx) {
        info(*ctx, "solution reads from file {}\n", filename);
        std::ifstream ifs(filename);
        result ret;
        ifs >> ret;
        return ret;
    }

    return result{ result_status::empty_context };
}

/**
 * @brief Convert a string_view into a integer.
 *
 * @note waiting for std::fron_chars or boost::qi dependencies
 */
constexpr inline std::optional<int>
to_int(std::string_view s) noexcept
{
    constexpr std::size_t size_limit = 512;

    if (s.size() > size_limit)
        return std::nullopt;

    char buffer[size_limit + 1] = { '\0' };
    std::size_t i = 0;
    std::size_t e = std::min(s.size(), size_limit);

    for (i = 0; i != e; ++i)
        buffer[i] = s[i];

    int result = 0;
    if (auto read = std::sscanf(buffer, "%d", &result); read)
        return result;
    else
        return std::nullopt;
}

std::istream&
operator>>(std::istream& is, result& ret)
{
    ret.strings = std::make_shared<baryonyx::string_buffer>();

    ret.status = result_status::success;
    ret.solutions.emplace_back();
    int line{ 0 };

    std::string buffer;
    while (is.good()) {
        std::getline(is, buffer);

        if (buffer.empty() && is.eof())
            return is;

        std::string::size_type i{ 0 };
        std::string::size_type e{ buffer.size() };

        for (; i != e; ++i)
            if (!std::isspace(buffer[i]))
                break;

        if (i != e && buffer[i] == '\\')
            continue;

        auto found = buffer.find('=', 0);
        if (found == std::string::npos)
            throw file_format_failure(
              file_format_error_tag::bad_name, line, 0);

        std::string_view left(buffer.data(), found);
        std::string_view right(buffer.data() + found + 1,
                               buffer.size() - found);

        if (auto value = to_int(right); value.has_value()) {
            ret.variable_name.emplace_back(ret.strings->append(left));
            ret.solutions.back().variables.emplace_back(!!(*value));
        } else {
            throw file_format_failure(
              file_format_error_tag::bad_name, line, 0);
        }
        ++line;
    }

    return is;
}

inline std::ostream&
operator<<(std::ostream& os, const best_solution_writer& writer)
{
    std::size_t i = 0, e = writer.res.variable_name.size();

    for (; i != e; ++i)
        os << writer.res.variable_name[i] << ": "
           << (writer.res.solutions.back().variables[i] ? 1 : 0) << '\n';

    return os;
}
}
