/* Copyright (C) 2018-2021 INRAE
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

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>

#include "debug.hpp"
#include "private.hpp"
#include "utils.hpp"

#include <climits>

constexpr char
to_lower(const char c) noexcept
{
    return (c >= 'A' && c <= 'Z') ? static_cast<char>(c - 'A' + 'a') : c;
}

constexpr bool
are_equal(const std::string_view lhs, const std::string_view rhs) noexcept
{
    if (rhs.size() != lhs.size())
        return false;

    for (std::string_view::size_type i = 0, e = lhs.size(); i != e; ++i)
        if (to_lower(lhs[i]) != to_lower(rhs[i]))
            return false;

    return true;
}

constexpr bool
is_char_operator(const int c) noexcept
{
    return c == '<' || c == '>' || c == '=';
}

constexpr bool
starts_with_operator(const std::string_view str) noexcept
{
    return (!str.empty() && is_char_operator(str[0]));
}

constexpr bool
is_digit(const int c) noexcept
{
    return c >= static_cast<int>('0') && c <= static_cast<int>('9');
}

constexpr bool
is_char_name(const int c) noexcept
{
    if ((c >= static_cast<int>('a') && c <= static_cast<int>('z')) ||
        (c >= static_cast<int>('A') && c <= static_cast<int>('Z')) ||
        is_digit(c))
        return true;

    switch (c) {
    case '!':
    case '"':
    case '#':
    case '$':
    case '%':
    case '&':
    case '(':
    case ')':
    case ',':
    case '.':
    case ';':
    case '?':
    case '@':
    case '_':
    case '{':
    case '}':
    case '~':
        return true;
    default:
        return false;
    }
}

constexpr bool
starts_with_name(const std::string_view str) noexcept
{
    return (!str.empty() && is_char_name(str[0]));
}

constexpr bool
is_number(const int c) noexcept
{
    return is_digit(c) || c == '.' || c == 'e' || c == 'E' || c == '-' ||
           c == '+';
}

constexpr bool
starts_with_a_number(const int c) noexcept
{
    return is_digit(c) || c == '.' || c == '-' || c == '+';
}

constexpr bool
starts_with_number(const std::string_view str) noexcept
{
    if (str.empty())
        return false;

    if (str[0] == 'i' || str[0] == 'I')
        if (are_equal(str, "infinity") || are_equal(str, "inf"))
            return true;

    return starts_with_a_number(str[0]);
}

constexpr bool
is_separator(const int c) noexcept
{
    switch (c) {
    case '<':
    case '=':
    case '>':
    case ':':
    case '-':
    case '+':
    case '[':
    case ']':
    case '*':
    case '^':
        return true;
    default:
        return false;
    }
}

constexpr std::string_view keywords[] = {
    "binary",  "binaries", "bin",  "bound", "bounds", "general",
    "generals", "gen",     "end",  "st",    "subject", "such",
    "s.t.",    "s.t",      "st."
};

constexpr bool
is_keyword(const std::string_view str) noexcept
{
    for (auto i : keywords)
        if (are_equal(i, str))
            return true;

    return false;
}

constexpr bool
is_variable_name(const std::string_view str) noexcept
{
    if (!starts_with_name(str) || is_keyword(str))
        return false;

    for (const char c : str)
        if (!is_char_name(c))
            return false;

    return true;
}

struct operator_token
{
    constexpr operator_token() noexcept = default;

    constexpr operator_token(baryonyx::operator_type op_, int read_) noexcept
      : op(op_)
      , read(read_)
    {}

    baryonyx::operator_type op = baryonyx::operator_type::equal;
    int read = 0;
};

struct real_token
{
    constexpr real_token() noexcept = default;

    constexpr real_token(double value_, int read_) noexcept
      : value(value_)
      , read(read_)
    {}

    double value = 0.0;
    int read = 0;
};

struct function_element_token
{
    constexpr function_element_token() noexcept = default;

    constexpr function_element_token(double factor_,
                                     std::string_view name_,
                                     int read_) noexcept
      : factor(factor_)
      , name(name_)
      , read(read_)
    {}

    double factor = 0.0;
    std::string_view name = {};
    int read = 0;
};

/// A bound statement only updates the sides it mentions: `x <= 4` keeps the
/// current lower bound of `x`.
struct bound_token
{
    std::optional<double> min;
    std::optional<double> max;
    std::string_view name;
    int read = 0;
};

struct raw_problem_status
{
    baryonyx::file_format_error_tag tag;
    unsigned int line = 0u;
    unsigned int column = 0u;

    constexpr raw_problem_status(baryonyx::file_format_error_tag tag_) noexcept
      : tag(tag_)
    {}

    constexpr operator bool() const noexcept
    {
        return tag == baryonyx::file_format_error_tag::success;
    }
};

/// Splits the input stream into LP tokens and keeps a look-ahead window of
/// `stream_buffer_size` tokens.
///
/// Tokens are views into a ring of `stream_buffer_size` words read from the
/// stream. A new word always overwrites the oldest one, which is no longer
/// referenced by the window. A view removed by `pop_front` must therefore
/// never be used afterwards.
class stream_buffer
{
public:
    constexpr static int stream_buffer_size = 10;

    using string_view_array = std::array<std::string_view, stream_buffer_size>;

    stream_buffer(std::istream& is_)
      : is(is_)
    {
        for (int i = 0; i != stream_buffer_size; ++i)
            buffer_ptr[i] = next_token();
    }

    void pop_front()
    {
        pop_front(1);
    }

    void pop_front(int size)
    {
        bx_expects(size >= 1 && size < stream_buffer_size);

        std::copy(buffer_ptr.begin() + size,
                  buffer_ptr.end(),
                  buffer_ptr.begin());

        do_push_back(size);
    }

    const string_view_array& array() const noexcept
    {
        return buffer_ptr;
    }

    std::string_view first() const noexcept
    {
        return buffer_ptr[0];
    }

    std::string_view second() const noexcept
    {
        return buffer_ptr[1];
    }

    std::string_view third() const noexcept
    {
        return buffer_ptr[2];
    }

private:
    struct stream_token
    {
        std::string buffer;
        std::size_t current = 0;
    };

    string_view_array buffer_ptr;
    std::array<stream_token, stream_buffer_size> token;
    std::istream& is;
    int current_token_buffer = 0;

    void do_push_back(int size)
    {
        for (int i = stream_buffer_size - size; i < stream_buffer_size; ++i)
            buffer_ptr[i] = next_token();
    }

    bool read_word(std::string& word)
    {
        word.clear();

        while (is >> word) {
            const auto comment = word.find('\\');
            if (comment != std::string::npos) {
                is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                word.erase(comment);
            }

            if (!word.empty())
                return true;
        }

        word.clear();
        return false;
    }

    std::string_view next_token()
    {
        if (token[current_token_buffer].current >=
            token[current_token_buffer].buffer.size()) {
            current_token_buffer =
              (current_token_buffer + 1) % stream_buffer_size;

            token[current_token_buffer].current = 0;
            if (!read_word(token[current_token_buffer].buffer))
                return std::string_view();
        }

        auto& tok = token[current_token_buffer];
        const std::string_view word = tok.buffer;
        const std::size_t start = tok.current++;

        if (is_separator(word[start]))
            return word.substr(start, 1);

        const bool is_numeric = starts_with_a_number(word[start]);

        for (; tok.current < word.size(); ++tok.current) {
            const char c = word[tok.current];

            const bool is_exponent_sign =
              is_numeric && (c == '-' || c == '+') &&
              (word[tok.current - 1] == 'e' || word[tok.current - 1] == 'E');

            if (is_exponent_sign)
                continue;

            if (is_separator(c))
                break;

            if (is_numeric && !is_number(c))
                break;
        }

        return word.substr(start, tok.current - start);
    }
};

std::optional<int>
to_int_exact(const double value) noexcept
{
    if (!std::isfinite(value) || value != std::trunc(value) ||
        value < static_cast<double>(std::numeric_limits<int>::min()) ||
        value > static_cast<double>(std::numeric_limits<int>::max()))
        return std::nullopt;

    return static_cast<int>(value);
}

int
to_int_saturated(const double value) noexcept
{
    if (value <= static_cast<double>(std::numeric_limits<int>::min()))
        return std::numeric_limits<int>::min();

    if (value >= static_cast<double>(std::numeric_limits<int>::max()))
        return std::numeric_limits<int>::max();

    return static_cast<int>(value);
}

struct problem_parser
{
    baryonyx::raw_problem m_problem;
    std::unordered_map<std::string_view, int> m_cache_variable;

    problem_parser()
    {
        m_problem.strings = std::make_shared<baryonyx::string_buffer>();
    }

    int get_or_assign_variable(const std::string_view value)
    {
        if (auto it = m_cache_variable.find(value);
            it != m_cache_variable.end())
            return it->second;

        if (m_cache_variable.size() >= baryonyx::to_unsigned(INT_MAX))
            return -1;

        auto size = static_cast<int>(m_cache_variable.size());
        auto string = m_problem.strings->append(value);

        m_cache_variable.try_emplace(string, size);

        m_problem.vars.names.emplace_back(string);
        m_problem.vars.values.emplace_back(
          0, std::numeric_limits<int>::max(), baryonyx::variable_type::real);

        return size;
    }

    [[nodiscard]] bool append_to_objective(const function_element_token& elem)
    {
        if (elem.name.empty()) {
            m_problem.objective.value += elem.factor;
            return true;
        }

        auto id = get_or_assign_variable(elem.name);
        if (id == -1)
            return false;

        auto it = std::find_if(
          std::begin(m_problem.objective.elements),
          std::end(m_problem.objective.elements),
          [id](const auto& elem) { return elem.variable_index == id; });

        if (it == std::end(m_problem.objective.elements))
            m_problem.objective.elements.emplace_back(elem.factor, id);
        else
            it->factor += elem.factor;

        return true;
    }

    void append_to_quadratic_objective(int id1, int id2, double factor)
    {
        auto it = std::find_if(m_problem.objective.qelements.begin(),
                               m_problem.objective.qelements.end(),
                               [id1, id2](const auto& elem) {
                                   return (elem.variable_index_a == id1 &&
                                           elem.variable_index_b == id2) ||
                                          (elem.variable_index_a == id2 &&
                                           elem.variable_index_b == id1);
                               });

        if (it == m_problem.objective.qelements.end()) {
            auto& ret = m_problem.objective.qelements.emplace_back();
            ret.factor = factor;
            ret.variable_index_a = id1;
            ret.variable_index_b = id2;
        } else {
            it->factor += factor;
        }
    }

    [[nodiscard]] bool append_to_constraint(
      std::vector<baryonyx::objective_function_element>& elements,
      const function_element_token& elem)
    {
        auto id = get_or_assign_variable(elem.name);
        if (id == -1)
            return false;

        auto it = std::find_if(
          elements.begin(), elements.end(), [id](const auto& elem) {
              return elem.variable_index == id;
          });

        if (it == elements.end())
            elements.emplace_back(elem.factor, id);
        else
            it->factor += elem.factor;

        return true;
    }

    [[nodiscard]] bool set_boolean_variable(const std::string_view value)
    {
        auto id = get_or_assign_variable(value);
        if (id < 0)
            return false;

        m_problem.vars.values[id].type = baryonyx::variable_type::binary;
        m_problem.vars.values[id].min = 0;
        m_problem.vars.values[id].max = 1;

        return true;
    }

    [[nodiscard]] bool set_integer_variable(const std::string_view value)
    {
        auto id = get_or_assign_variable(value);
        if (id < 0)
            return false;

        m_problem.vars.values[id].type = baryonyx::variable_type::general;

        return true;
    }

    /// Bounds are rounded inward (`ceil` for the lower bound, `floor` for the
    /// upper bound) because the solver only handles integer variables.
    [[nodiscard]] bool set_bound_variable(const bound_token& bound)
    {
        auto id = get_or_assign_variable(bound.name);
        if (id < 0)
            return false;

        if (bound.min)
            m_problem.vars.values[id].min =
              to_int_saturated(std::ceil(*bound.min));

        if (bound.max)
            m_problem.vars.values[id].max =
              to_int_saturated(std::floor(*bound.max));

        return true;
    }
};

std::optional<double>
read_real(const std::string_view buf) noexcept
{
    if (are_equal(buf, "inf") || are_equal(buf, "infinity"))
        return std::numeric_limits<double>::infinity();

    double result = 0.0;
    const auto* const end = buf.data() + buf.size();
    const auto [ptr, ec] = std::from_chars(buf.data(), end, result);

    if (ec != std::errc{} || ptr != end || std::isnan(result))
        return std::nullopt;

    return result;
}

std::optional<real_token>
read_signed_real(const std::string_view buf_1,
                 const std::string_view buf_2) noexcept
{
    if (buf_1 == "-" || buf_1 == "+") {
        const double sign = buf_1 == "-" ? -1.0 : 1.0;

        if (auto value = read_real(buf_2); value)
            return real_token(sign * *value, 2);

        return std::nullopt;
    }

    if (auto value = read_real(buf_1); value)
        return real_token(*value, 1);

    return std::nullopt;
}

std::optional<operator_token>
read_operator(const std::string_view buf_1,
              const std::string_view buf_2) noexcept
{
    if (buf_1 == "<")
        return operator_token{ baryonyx::operator_type::less,
                               buf_2 == "=" ? 2 : 1 };

    if (buf_1 == ">")
        return operator_token{ baryonyx::operator_type::greater,
                               buf_2 == "=" ? 2 : 1 };

    if (buf_1 == "=") {
        if (buf_2 == "<")
            return operator_token{ baryonyx::operator_type::less, 2 };

        if (buf_2 == "=")
            return operator_token{ baryonyx::operator_type::equal, 2 };

        if (buf_2 == ">")
            return operator_token{ baryonyx::operator_type::greater, 2 };

        return operator_token{ baryonyx::operator_type::equal, 1 };
    }

    return std::nullopt;
}

/// Reads `[sign] [number] [name]`. Without a name the element is a constant
/// term, which requires a number.
std::optional<function_element_token>
read_function_element(const std::string_view buf_1,
                      const std::string_view buf_2,
                      const std::string_view buf_3) noexcept
{
    const std::array<std::string_view, 3> tokens{ buf_1, buf_2, buf_3 };

    double factor = 1.0;
    std::size_t read = 0;

    if (tokens[0] == "-" || tokens[0] == "+") {
        factor = tokens[0] == "-" ? -1.0 : 1.0;
        ++read;
    }

    bool has_number = false;
    if (auto value = read_real(tokens[read]); value) {
        factor *= *value;
        has_number = true;
        ++read;
    }

    if (read < tokens.size() && is_variable_name(tokens[read]))
        return function_element_token(
          factor, tokens[read], static_cast<int>(read + 1));

    if (has_number)
        return function_element_token(
          factor, std::string_view(), static_cast<int>(read));

    return std::nullopt;
}

std::optional<bound_token>
read_bound(const stream_buffer::string_view_array& tokens) noexcept
{
    const auto apply = [](bound_token& bound,
                          baryonyx::operator_type op,
                          double value,
                          bool variable_on_left) {
        switch (op) {
        case baryonyx::operator_type::equal:
            bound.min = value;
            bound.max = value;
            break;
        case baryonyx::operator_type::less:
            (variable_on_left ? bound.max : bound.min) = value;
            break;
        case baryonyx::operator_type::greater:
            (variable_on_left ? bound.min : bound.max) = value;
            break;
        }
    };

    bound_token bound;
    std::size_t read = 0;

    if (starts_with_number(tokens[0])) {
        auto value = read_signed_real(tokens[0], tokens[1]);
        if (!value)
            return std::nullopt;
        read += static_cast<std::size_t>(value->read);

        auto op = read_operator(tokens[read], tokens[read + 1]);
        if (!op)
            return std::nullopt;
        read += static_cast<std::size_t>(op->read);

        if (!is_variable_name(tokens[read]))
            return std::nullopt;
        bound.name = tokens[read++];

        apply(bound, op->op, value->value, false);
    } else if (is_variable_name(tokens[0])) {
        bound.name = tokens[read++];

        if (are_equal(tokens[read], "free")) {
            bound.min = -std::numeric_limits<double>::infinity();
            bound.max = std::numeric_limits<double>::infinity();
            bound.read = static_cast<int>(read + 1);
            return bound;
        }

        if (!starts_with_operator(tokens[read])) {
            bound.min = -std::numeric_limits<double>::infinity();
            bound.max = std::numeric_limits<double>::infinity();
            bound.read = static_cast<int>(read);
            return bound;
        }
    } else {
        return std::nullopt;
    }

    if (starts_with_operator(tokens[read])) {
        auto op = read_operator(tokens[read], tokens[read + 1]);
        if (!op)
            return std::nullopt;
        read += static_cast<std::size_t>(op->read);

        auto value = read_signed_real(tokens[read], tokens[read + 1]);
        if (!value)
            return std::nullopt;
        read += static_cast<std::size_t>(value->read);

        apply(bound, op->op, value->value, true);
    }

    if (bound.min && *bound.min == std::numeric_limits<double>::infinity())
        return std::nullopt;

    if (bound.max && *bound.max == -std::numeric_limits<double>::infinity())
        return std::nullopt;

    if (bound.min && bound.max && *bound.min > *bound.max)
        return std::nullopt;

    bound.read = static_cast<int>(read);
    return bound;
}

std::optional<baryonyx::objective_function_type>
read_objective_type(const std::string_view token) noexcept
{
    if (are_equal(token, "maximize") || are_equal(token, "maximum") ||
        are_equal(token, "max"))
        return baryonyx::objective_function_type::maximize;

    if (are_equal(token, "minimize") || are_equal(token, "minimum") ||
        are_equal(token, "min"))
        return baryonyx::objective_function_type::minimize;

    return std::nullopt;
}

constexpr int
read_subject_to(const std::string_view buf_1,
                const std::string_view buf_2,
                const std::string_view buf_3) noexcept
{
    if (are_equal(buf_1, "st") || are_equal(buf_1, "st.") ||
        are_equal(buf_1, "s.t") || are_equal(buf_1, "s.t."))
        return are_equal(buf_2, ":") ? 2 : 1;

    if (are_equal(buf_1, "subject"))
        if (are_equal(buf_2, "to"))
            return are_equal(buf_3, ":") ? 3 : 2;

    if (are_equal(buf_1, "such"))
        if (are_equal(buf_2, "that"))
            return are_equal(buf_3, ":") ? 3 : 2;

    return 0;
}

constexpr int
read_bounds(const std::string_view buf_1,
            const std::string_view buf_2) noexcept
{
    if (are_equal(buf_1, "bounds") || are_equal(buf_1, "bound"))
        return are_equal(buf_2, ":") ? 2 : 1;

    return 0;
}

constexpr int
read_binary(const std::string_view buf_1,
            const std::string_view buf_2) noexcept
{
    if (are_equal(buf_1, "binary") || are_equal(buf_1, "binaries") ||
        are_equal(buf_1, "bin"))
        return are_equal(buf_2, ":") ? 2 : 1;

    return 0;
}

constexpr int
read_general(const std::string_view buf_1,
             const std::string_view buf_2) noexcept
{
    if (are_equal(buf_1, "general") || are_equal(buf_1, "generals") ||
        are_equal(buf_1, "gen"))
        return are_equal(buf_2, ":") ? 2 : 1;

    return 0;
}

constexpr int
read_end(const std::string_view buf_1, const std::string_view buf_2) noexcept
{
    if (are_equal(buf_1, "end"))
        return are_equal(buf_2, ":") ? 2 : 1;

    return 0;
}

bool
starts_with_quadratic(const std::string_view buf_1,
                      const std::string_view buf_2) noexcept
{
    return buf_1 == "[" || (buf_1 == "+" && buf_2 == "[") ||
           (buf_1 == "-" && buf_2 == "[");
}

/// Reads `[ q_1 + q_2 ... ] / 2` where each term is `[c] x * y` or `[c] x ^ 2`.
raw_problem_status
read_quadratic_objective(stream_buffer& buf,
                         problem_parser& p,
                         const double sign_factor)
{
    if (buf.first() != "[")
        return baryonyx::file_format_error_tag::bad_objective_quadratic;

    buf.pop_front();

    while (!buf.first().empty() && buf.first() != "]") {
        const auto fct =
          read_function_element(buf.first(), buf.second(), buf.third());
        if (!fct || fct->name.empty())
            return baryonyx::file_format_error_tag::bad_objective_quadratic;

        const auto id1 = p.get_or_assign_variable(fct->name);
        if (id1 == -1)
            return baryonyx::file_format_error_tag::too_many_variables;

        const double factor = fct->factor * sign_factor / 2.0;
        buf.pop_front(fct->read);

        if (buf.first() == "*") {
            if (!is_variable_name(buf.second()))
                return baryonyx::file_format_error_tag::
                  bad_objective_quadratic;

            const auto id2 = p.get_or_assign_variable(buf.second());
            if (id2 == -1)
                return baryonyx::file_format_error_tag::too_many_variables;

            p.append_to_quadratic_objective(id1, id2, factor);
            buf.pop_front(2);
        } else if (buf.first() == "^" && buf.second() == "2") {
            p.append_to_quadratic_objective(id1, id1, factor);
            buf.pop_front(2);
        } else {
            return baryonyx::file_format_error_tag::bad_objective_quadratic;
        }
    }

    if (buf.first() != "]")
        return baryonyx::file_format_error_tag::bad_objective_quadratic;

    buf.pop_front();

    if (buf.first() == "/" && buf.second() == "2")
        buf.pop_front(2);
    else if (buf.first() == "/2")
        buf.pop_front(1);
    else
        return baryonyx::file_format_error_tag::bad_objective_quadratic;

    return baryonyx::file_format_error_tag::success;
}

raw_problem_status
read_objective(stream_buffer& buf, problem_parser& p)
{
    if (auto obj_type = read_objective_type(buf.first()); !obj_type)
        return baryonyx::file_format_error_tag::bad_objective_function_type;
    else
        p.m_problem.type = *obj_type;

    buf.pop_front();

    if (buf.first() == ":")
        buf.pop_front();

    if (is_variable_name(buf.first()) && buf.second() == ":")
        buf.pop_front(2);

    while (!buf.first().empty() && !is_keyword(buf.first())) {
        if (starts_with_quadratic(buf.first(), buf.second())) {
            double factor = 1.0;

            if (buf.first() == "-") {
                factor = -1.0;
                buf.pop_front();
            } else if (buf.first() == "+") {
                buf.pop_front();
            }

            if (auto ret = read_quadratic_objective(buf, p, factor); !ret)
                return ret;

            continue;
        }

        auto elem =
          read_function_element(buf.first(), buf.second(), buf.third());
        if (!elem)
            return baryonyx::file_format_error_tag::bad_objective;

        if (!p.append_to_objective(*elem))
            return baryonyx::file_format_error_tag::too_many_variables;

        buf.pop_front(elem->read);
    }

    return baryonyx::file_format_error_tag::success;
}

/// The coefficients and the right-hand side are stored as integers. Since the
/// variables and the coefficients are integers, the left-hand side is an
/// integer and a fractional right-hand side can be rounded without changing
/// the feasible set: `floor` for `<=`, `ceil` for `>=`.
std::optional<int>
round_constraint_value(const baryonyx::operator_type op, const double value)
{
    switch (op) {
    case baryonyx::operator_type::less:
        return to_int_exact(std::floor(value));
    case baryonyx::operator_type::greater:
        return to_int_exact(std::ceil(value));
    case baryonyx::operator_type::equal:
        return to_int_exact(value);
    }

    return std::nullopt;
}

raw_problem_status
read_constraint(stream_buffer& buf, problem_parser& p, int id)
{
    std::string_view label;
    std::vector<baryonyx::objective_function_element> terms;
    double constant = 0.0;

    if (is_variable_name(buf.first()) && buf.second() == ":") {
        label = p.m_problem.strings->append(buf.first());
        buf.pop_front(2);
    }

    do {
        auto elem =
          read_function_element(buf.first(), buf.second(), buf.third());
        if (!elem)
            return baryonyx::file_format_error_tag::bad_constraint;

        if (elem->name.empty())
            constant += elem->factor;
        else if (!p.append_to_constraint(terms, *elem))
            return baryonyx::file_format_error_tag::too_many_variables;

        buf.pop_front(elem->read);
    } while (!buf.first().empty() && !starts_with_operator(buf.first()));

    if (terms.empty())
        return baryonyx::file_format_error_tag::bad_constraint;

    std::vector<baryonyx::function_element> elements;
    elements.reserve(terms.size());
    for (const auto& term : terms) {
        const auto factor = to_int_exact(term.factor);
        if (!factor)
            return baryonyx::file_format_error_tag::bad_constraint;

        elements.emplace_back(*factor, term.variable_index);
    }

    auto operator_opt = read_operator(buf.first(), buf.second());
    if (!operator_opt)
        return baryonyx::file_format_error_tag::bad_constraint;
    buf.pop_front(operator_opt->read);

    auto value_opt = read_signed_real(buf.first(), buf.second());
    if (!value_opt)
        return baryonyx::file_format_error_tag::bad_constraint;
    buf.pop_front(value_opt->read);

    if (starts_with_operator(buf.first()))
        return baryonyx::file_format_error_tag::bad_constraint;

    auto value =
      round_constraint_value(operator_opt->op, value_opt->value - constant);
    if (!value)
        return baryonyx::file_format_error_tag::bad_constraint;

    switch (operator_opt->op) {
    case baryonyx::operator_type::equal:
        p.m_problem.equal_constraints.emplace_back(
          label, std::move(elements), *value, id);
        break;
    case baryonyx::operator_type::greater:
        p.m_problem.greater_constraints.emplace_back(
          label, std::move(elements), *value, id);
        break;
    case baryonyx::operator_type::less:
        p.m_problem.less_constraints.emplace_back(
          label, std::move(elements), *value, id);
        break;
    }

    return baryonyx::file_format_error_tag::success;
}

raw_problem_status
read_constraints(stream_buffer& buf, problem_parser& p)
{
    int constraint_next_id = 0;

    while (!buf.first().empty() && !is_keyword(buf.first())) {
        if (constraint_next_id == INT_MAX)
            return baryonyx::file_format_error_tag::too_many_constraints;

        if (auto ret = read_constraint(buf, p, constraint_next_id); !ret)
            return ret;

        ++constraint_next_id;
    }

    return baryonyx::file_format_error_tag::success;
}

raw_problem_status
read_bounds_section(stream_buffer& buf, problem_parser& p)
{
    while (!buf.first().empty() && !is_keyword(buf.first())) {
        auto bound_opt = read_bound(buf.array());
        if (!bound_opt)
            return baryonyx::file_format_error_tag::bad_bound;

        if (!p.set_bound_variable(*bound_opt))
            return baryonyx::file_format_error_tag::too_many_variables;

        buf.pop_front(bound_opt->read);
    }

    return baryonyx::file_format_error_tag::success;
}

template<typename Function>
raw_problem_status
read_variables_section(stream_buffer& buf,
                       const baryonyx::file_format_error_tag error,
                       Function set_variable)
{
    while (!buf.first().empty() && !is_keyword(buf.first())) {
        if (!is_variable_name(buf.first()))
            return error;

        if (!set_variable(buf.first()))
            return baryonyx::file_format_error_tag::too_many_variables;

        buf.pop_front();
    }

    return baryonyx::file_format_error_tag::success;
}

raw_problem_status
parse([[maybe_unused]] const baryonyx::context& ctx,
      stream_buffer& buf,
      problem_parser& p)
{
    if (auto ret = read_objective(buf, p); !ret)
        return ret;

    if (auto read = read_subject_to(buf.first(), buf.second(), buf.third());
        read) {
        buf.pop_front(read);

        if (auto ret = read_constraints(buf, p); !ret)
            return ret;
    }

    for (;;) {
        if (auto read = read_bounds(buf.first(), buf.second()); read) {
            buf.pop_front(read);

            if (auto ret = read_bounds_section(buf, p); !ret)
                return ret;

            continue;
        }

        if (auto read = read_binary(buf.first(), buf.second()); read) {
            buf.pop_front(read);

            if (auto ret = read_variables_section(
                  buf,
                  baryonyx::file_format_error_tag::bad_binary,
                  [&p](auto name) { return p.set_boolean_variable(name); });
                !ret)
                return ret;

            continue;
        }

        if (auto read = read_general(buf.first(), buf.second()); read) {
            buf.pop_front(read);

            if (auto ret = read_variables_section(
                  buf,
                  baryonyx::file_format_error_tag::bad_general,
                  [&p](auto name) { return p.set_integer_variable(name); });
                !ret)
                return ret;

            continue;
        }

        break;
    }

    if (auto read = read_end(buf.first(), buf.second()); read) {
        buf.pop_front(read);

        return buf.first().empty()
                 ? baryonyx::file_format_error_tag::success
                 : baryonyx::file_format_error_tag::bad_end_of_file;
    }

    return baryonyx::file_format_error_tag::bad_end;
}

namespace baryonyx {

raw_problem
make_problem(const context& ctx, std::istream& is) noexcept
{
    stream_buffer sb(is);
    problem_parser p;

    if (auto ret = parse(ctx, sb, p); ret)
        return std::move(p.m_problem);
    else
        return raw_problem(ret.tag);
}

} // namespace baryonyx
