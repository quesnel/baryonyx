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

#include "debug.hpp"
#include "dynarray.hpp"
#include "main.hpp"
#include "utils.hpp"

#include <fstream>
#include <limits>
#include <locale>
#include <numeric>
#include <optional>
#include <unordered_map>
#include <variant>

#include <fmt/color.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include <cmath>
#include <utility>

struct csv_whitespace : std::ctype<char>
{
    static mask* make_table()
    {
        auto* v = new mask[table_size];
        std::copy_n(classic_table(), table_size, v);

        // Comma will be classified as whitespace.
        v[static_cast<int>(',')] |= space;

        // Space will be classified as character with graphic representation.
        v[static_cast<int>(' ')] = graph;

        return v;
    }

    csv_whitespace(std::size_t refs = 0)
      : ctype(make_table(), true, refs)
    {}
};

static bool
have_lp_extension(std::string filename)
{
    auto pos = filename.rfind('.');

    if (pos == std::string::npos)
        return false;

    return filename.compare(pos, 3, ".lp") == 0;
}

static inline int
get_digits_number(double x) noexcept
{
    bx_expects(x < static_cast<double>(std::numeric_limits<int>::max()));

    return static_cast<int>(std::floor(std::log10(std::abs(x)))) + 5;
}

constexpr bool
is_valid(double d) noexcept
{
    return d != -std::numeric_limits<double>::infinity() &&
           d != std::numeric_limits<double>::infinity();
}

struct success
{};

struct open_file_error
{
    std::string file_name;
    int errc;
};

struct read_header_error
{
    int line;
};

struct read_data_error
{
    int line;
};

struct model_name_error
{
    std::string token;
    int line;
};

struct solver_name_error
{
    std::string token;
    int line;
    std::string::size_type column;
};

struct open_result_file_error
{
    std::string file_name;
    int errc;
};

using benchmark_result = std::variant<success,
                                      open_file_error,
                                      read_header_error,
                                      read_data_error,
                                      model_name_error,
                                      solver_name_error,
                                      open_result_file_error>;

struct bench
{
    class model
    {
        model(std::string_view name_)
          : m_name(name_)
        {}

        std::string m_name;

    public:
        static std::optional<model> make_model(std::string_view name)
        {
            auto trim_name = baryonyx::trim(name);
            if (trim_name.empty())
                return std::nullopt;

            return model(name);
        }

        std::string_view name() const
        {
            return m_name;
        }

        std::size_t size() const noexcept
        {
            return m_name.size();
        }

        std::string filename() const
        {
            return have_lp_extension(m_name) ? m_name
                                             : fmt::format("{}.lp", m_name);
        }
    };

    class solver
    {
        solver(std::string_view name_)
          : m_name(name_)
        {}

        std::string m_name;

    public:
        static std::optional<solver> make_solver(std::string_view name)
        {
            auto trim_name = baryonyx::trim(name);
            if (trim_name.empty())
                return std::nullopt;

            return solver(name);
        }

        std::string_view name() const
        {
            return m_name;
        }
    };

    struct element
    {
        element() = default;

        element(double solution_)
          : solution(solution_)
          , success(true)
        {}

        double solution = std::numeric_limits<double>::infinity();
        bool success = false;
    };

    std::vector<model> models;
    std::vector<solver> solvers;
    baryonyx::dynamic_array<double> array;

    bool push_back_model(std::string_view name)
    {
        auto mdl = model::make_model(name);

        if (!mdl)
            return false;

        models.emplace_back(*mdl);

        return true;
    }

    bool push_back_solver(std::string_view name)
    {
        auto slv = solver::make_solver(name);
        if (!slv)
            return false;

        solvers.emplace_back(*slv);

        return true;
    }

    void save(std::ostream& os,
              std::string name,
              const std::vector<element>& c)
    {
        fmt::print(os, "file");

        for (const auto& elem : solvers)
            fmt::print(os, ",{}", elem.name());
        fmt::print(os, ",{}\n", name);

        for (std::size_t i{ 0 }, e{ models.size() }; i != e; ++i) {
            fmt::print(os, "{},", models[i].name());

            for (std::size_t j{ 0 }, end_j{ solvers.size() }; j != end_j; ++j)
                fmt::print(os, "{:.10g},", array(i, j));

            fmt::print(os, "{:.10g}\n", c[i].solution);
        }
    }

    benchmark_result load(std::istream& is)
    {
        std::string buffer;
        buffer.reserve(BUFSIZ);

        if (!std::getline(is, buffer))
            return read_header_error{ 1 };

        int line_pos = 1;

        {
            // Tries to read the first column (lp file) and forget the string
            std::string::size_type begin{ 0 };
            std::string::size_type end{ buffer.find(',', begin) };

            // If more than one column exists, we read solvers.
            if (end != std::string::npos) {
                do {
                    begin = end + 1;
                    end = buffer.find(',', begin);

                    if (end == std::string::npos) {
                        if (!push_back_solver(std::string_view(
                              buffer.data() + begin, buffer.size() - begin)))
                            return solver_name_error{ std::string(
                                                        buffer.data() + begin,
                                                        buffer.size() - begin),
                                                      1,
                                                      begin };

                    } else if (begin + 1 != end) {
                        if (!push_back_solver(std::string_view(
                              buffer.data() + begin, end - begin)))
                            return solver_name_error{ std::string(
                                                        buffer.data() + begin,
                                                        end - begin),
                                                      1,
                                                      begin };
                    } else {
                        return solver_name_error{
                            std::string(buffer.data() + begin), 1, begin
                        };
                    }
                } while (end != std::string::npos);
            }
        }

        ++line_pos;

        if (!solvers.empty())
            array.init(solvers.size());

        is.imbue(std::locale(is.getloc(), new csv_whitespace));

        while (is) {
            is >> buffer;
            if (!is.good())
                break;

            if (buffer.empty())
                return model_name_error{ buffer, line_pos };

            if (!push_back_model(std::string_view(buffer)))
                return model_name_error{ buffer, line_pos };

            if (const auto length{ solvers.size() }; length > 0) {
                array.push_line();

                for (std::size_t i{ 0 }; i != length; ++i) {
                    is >> buffer;
                    if (!is || buffer.empty())
                        return read_data_error{ line_pos };

                    array(array.rows() - 1, i) =
                      to_double(baryonyx::trim(std::string_view(buffer)),
                                std::numeric_limits<double>::infinity());
                }
            }

            ++line_pos;
        }

        return success{};
    }

    void show_to_console(const std::vector<element>& current, std::string name)
    {
        const auto model_size{ models.size() };
        const auto solver_size{ solvers.size() };

        std::vector<int> row_length(solvers.size(), 5);
        int current_row_length = { 10 };

        for (size_t i{ 0 }; i != model_size; ++i) {
            for (size_t j{ 0 }; j != solver_size; ++j) {
                if (is_valid(array(i, j))) {
                    auto digits = get_digits_number(array(i, j));

                    if (digits > row_length[j])
                        row_length[j] = digits;
                }
            }
        }

        for (size_t i{ 0 }; i != model_size; ++i) {
            if (current[i].success && is_valid(current[i].solution)) {
                auto digits = get_digits_number(current[i].solution);

                if (digits > current_row_length)
                    current_row_length = digits;
            }
        }

        auto modelname_lenght_max =
          std::max(static_cast<std::string::size_type>(4),
                   std::max_element(models.cbegin(),
                                    models.cend(),
                                    [](const auto& lhs, const auto& rhs) {
                                        return lhs.size() < rhs.size();
                                    })
                     ->size());

        fmt::print("{:^{}} | ", "file", modelname_lenght_max);

        for (std::size_t solver{ 0 }; solver != solver_size; ++solver)
            fmt::print("{:>{}} ",
                       solvers[solver].name().substr(0, 10),
                       row_length[solver]);

        fmt::print("{:>{}}\n", name, current_row_length);

        for (size_t i{ 0 }; i != model_size; ++i) {
            fmt::print("{:^{}} | ", models[i].name(), modelname_lenght_max);

            auto lower = std::numeric_limits<double>::max();
            auto upper = std::numeric_limits<double>::lowest();

            for (std::size_t solver{ 0 }, end_solver(solvers.size());
                 solver != end_solver;
                 ++solver) {
                if (is_valid(array(i, solver)) && lower > array(i, solver))
                    lower = array(i, solver);

                if (is_valid(array(i, solver)) && upper < array(i, solver))
                    upper = array(i, solver);
            }

            if (current[i].success) {
                if (is_valid(current[i].solution) &&
                    lower > current[i].solution)
                    lower = current[i].solution;

                if (is_valid(current[i].solution) &&
                    upper < current[i].solution)
                    upper = current[i].solution;
            }

            for (std::size_t solver{ 0 }, end_solver(solvers.size());
                 solver != end_solver;
                 ++solver) {

                if (array(i, solver) == lower)
                    fmt::print(fg(fmt::color::green),
                               "{:>{}.3f} ",
                               array(i, solver),
                               row_length[solver]);
                else if (array(i, solver) == upper)
                    fmt::print(fg(fmt::color::blue),
                               "{:>{}.3f} ",
                               array(i, solver),
                               row_length[solver]);
                else
                    fmt::print(
                      "{:>{}.3f} ", array(i, solver), row_length[solver]);
            }

            if (current[i].solution == lower)
                fmt::print(fg(fmt::color::green),
                           "{:>{}.3f} ",
                           current[i].solution,
                           current_row_length);
            else if (current[i].solution == upper)
                fmt::print(fg(fmt::color::blue),
                           "{:>{}.3f} ",
                           current[i].solution,
                           current_row_length);
            else
                fmt::print(
                  "{:>{}.3f} ", current[i].solution, current_row_length);

            fmt::print("\n");
        }
    }
};

static benchmark_result
try_benchmark(const baryonyx::context_ptr& ctx,
              std::string filepath,
              std::string name)
{
    std::ifstream ifs(filepath);
    if (!ifs.is_open())
        return open_file_error{ filepath, errno };

    bench b;
    auto ret = b.load(ifs);
    if (!std::holds_alternative<success>(ret))
        return ret;

    std::vector<bench::element> current(b.models.size());
    std::string dirname = filepath.substr(0, filepath.rfind(".csv")) + '/';

    for (std::size_t i{ 0 }, e{ b.models.size() }; i != e; ++i) {
        fmt::print("- optimizing {}\n", b.models[i].name());

        auto filename = b.models[i].filename();

        try {
            fmt::print("- optimizing: {} filename: {} in dirname: {}\n",
                       b.models[i].name(),
                       filename,
                       dirname);

            auto rawpb = baryonyx::make_problem(ctx, dirname + filename);
            auto result = baryonyx::optimize(ctx, rawpb);

            if (result)
                current[i] = result.solutions.back().value;
            else
                current[i] = std::numeric_limits<double>::infinity();

        } catch (const baryonyx::precondition_failure& e) {
            fmt::print(stderr, "internal failure: {}\n", e.what());
        } catch (const baryonyx::postcondition_failure& e) {
            fmt::print(stderr, "internal failure: {}\n", e.what());
        } catch (const baryonyx::numeric_cast_failure& e) {
            fmt::print(
              stderr, "numeric cast internal failure: {}\n", e.what());
        } catch (const baryonyx::file_access_failure& e) {
            fmt::print(stderr,
                       "file `{}' fail {}: {}\n",
                       e.file(),
                       e.error(),
                       std::strerror(e.error()));
        } catch (const baryonyx::file_format_failure& e) {
            fmt::print(stderr,
                       "file format error at line {} column {} "
                       "{}\n",
                       e.line(),
                       e.column(),
                       file_format_error_format(e.failure()));
        } catch (const baryonyx::problem_definition_failure& e) {
            fmt::print(stderr,
                       "definition problem error at {}: {}\n",
                       e.element(),
                       problem_definition_error_format(e.failure()));
        } catch (const baryonyx::solver_failure& e) {
            fmt::print(
              stderr, "solver error: {}\n", solver_error_format(e.failure()));
        } catch (const std::exception& e) {
            fmt::print(stderr, "failure: {}.\n", e.what());
        }
    }

    {
        // Compute the number of feasible solution found by the new solver.

        auto all_found =
          std::accumulate(current.cbegin(),
                          current.cend(),
                          0,
                          [](int cumul, const auto& elem) -> int {
                              return elem.success ? cumul + 1 : cumul;
                          });

        fmt::print("{} found {} feasible solution(s) / {}\n",
                   name,
                   all_found,
                   current.size());
    }

    {
        double mean_distance{ 0 };
        int optimum_number{ 0 };

        for (std::size_t i{ 1 }, e{ current.size() }; i < e; ++i) {
            if (is_valid(current[0].solution) &&
                is_valid(current[i].solution)) {
                ++optimum_number;

                mean_distance += ((current[0].solution - current[i].solution) /
                                  current[0].solution) *
                                 100.0;
            }
        }

        if (optimum_number > 0)
            fmt::print("Average distance from the optimum: {:.10g}%\n",
                       mean_distance / optimum_number);
    }

    {
        for (std::size_t solver{ 0 }, end_solver(b.solvers.size());
             solver != end_solver;
             ++solver) {

            double mean_distance{ 0 };
            int number{ 0 };

            for (std::size_t i{ 0 }, e{ current.size() }; i != e; ++i) {
                if (is_valid(b.array(i, solver)) &&
                    is_valid(current[i].solution)) {
                    ++number;

                    mean_distance +=
                      ((b.array(i, solver) - current[i].solution) /
                       b.array(i, solver)) *
                      100.0;
                }
            }

            if (number > 0)
                fmt::print("Average distance from the solution of {:.10g}: {:.10g}%\n",
                           b.solvers[solver].name(),
                           mean_distance /
                             static_cast<double>(current.size()));
        }
    }

    std::ofstream ofs(filepath);
    if (!ofs.is_open())
        return open_result_file_error{ filepath, errno };

    b.save(ofs, name, current);
    b.show_to_console(current, name);

    return success{};
}

template<class T>
struct always_false : std::false_type
{};

bool
benchmark(const baryonyx::context_ptr& ctx,
          std::string filepath,
          std::string name)
{
    auto ret = try_benchmark(ctx, filepath, name);

    return std::visit(
      [](auto&& err) -> bool {
          using T = std::decay_t<decltype(err)>;

          if constexpr (std::is_same_v<T, success>) {
              return true;
          } else if constexpr (std::is_same_v<T, open_file_error>) {
              fmt::print(stderr,
                         fmt::emphasis::bold |
                           fmt::fg(fmt::terminal_color::red),
                         "Can not open `{}' to start benchmark\n",
                         err.file_name);
              return false;
          } else if constexpr (std::is_same_v<T, read_header_error>) {
              fmt::print(stderr,
                         fmt::emphasis::bold |
                           fmt::fg(fmt::terminal_color::red),
                         "Can not read header at line `{}'\n",
                         err.line);
              return false;
          } else if constexpr (std::is_same_v<T, read_data_error>) {
              fmt::print(stderr,
                         fmt::emphasis::bold |
                           fmt::fg(fmt::terminal_color::red),
                         "Can not read data at line `{}'.\n",
                         err.line);

              return false;
          } else if constexpr (std::is_same_v<T, model_name_error>) {
              fmt::print(stderr,
                         fmt::emphasis::bold |
                           fmt::fg(fmt::terminal_color::red),
                         "Can not read solver name model at line `{}'. Token "
                         "read is `{}'\n",
                         err.line,
                         err.token);
              return false;
          } else if constexpr (std::is_same_v<T, solver_name_error>) {
              fmt::print(stderr,
                         fmt::emphasis::bold |
                           fmt::fg(fmt::terminal_color::red),
                         "Can not read solver name model at line `{}'. Token "
                         "read is `{}'\n",
                         err.line,
                         err.token);
              return false;
          } else if constexpr (std::is_same_v<T, open_result_file_error>) {
              fmt::print(stderr,
                         fmt::emphasis::bold |
                           fmt::fg(fmt::terminal_color::red),
                         "Can not open `{}' to record benchmark results\n",
                         err.file_name);
              return false;
          } else
              static_assert(always_false<T>::value, "non-exhaustive visitor!");
      },
      ret);
}
