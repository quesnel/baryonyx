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

#include "result.hpp"
#include "private.hpp"

#include <baryonyx/core-utils>

#include <fstream>
#include <string>

#include <cctype>

namespace baryonyx {

result
make_result(const baryonyx::context_ptr& ctx, const std::string& filename)
{
    info(ctx, "solution reads from file {}\n", filename);

    std::ifstream ifs(filename);

    result ret;
    ifs >> ret;

    return ret;
}

std::istream&
operator>>(std::istream& is, result& ret)
{
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

        auto it = buffer.find('=', 0);
        if (it == std::string::npos)
            throw file_format_failure(
              file_format_error_tag::bad_name, line, 0);

        std::string name = buffer.substr(0, it);
        int value;

        try {
            value = std::stoi(buffer.substr(it + 1));
        } catch (...) {
            throw file_format_failure(
              file_format_error_tag::bad_name, line, 0);
        }

        ret.variable_name.emplace_back(name);
        ret.solutions.back().variables.emplace_back(!!value);
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
