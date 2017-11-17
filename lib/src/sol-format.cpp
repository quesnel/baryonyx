/* Copyright (C) 2017 INRA
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

#include <baryonyx/core>

#include <istream>
#include <string>

#include <cctype>

namespace baryonyx_private {

baryonyx::result
read_result(std::istream& is)
{
    baryonyx::result ret;
    ret.status = baryonyx::result_status::success;
    int line{ 0 };

    std::string buffer;
    while (is.good()) {
        std::getline(is, buffer);

        if (buffer.empty() and is.eof())
            return ret;

        std::string::size_type i{ 0 };
        std::string::size_type e{ buffer.size() };

        for (; i != e; ++i)
            if (not std::isspace(buffer[i]))
                break;

        if (i != e and buffer[i] == '\\')
            continue;

        auto it = buffer.find('=', 0);
        if (it == std::string::npos)
            throw baryonyx::file_format_failure(
              baryonyx::file_format_error_tag::bad_name, line, 0);

        std::string name = buffer.substr(0, it);
        int value;

        try {
            value = std::stoi(buffer.substr(it + 1));
        } catch (...) {
            throw baryonyx::file_format_failure(
              baryonyx::file_format_error_tag::bad_name, line, 0);
        }

        ret.variable_name.emplace_back(name);
        ret.variable_value.emplace_back(value);
        ++line;
    }

    return ret;
}
}
