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

#include <fstream>
#include <limits>
#include <locale>
#include <numeric>
#include <unordered_map>

#include <fmt/color.h>
#include <fmt/format.h>
#include <fmt/ostream.h>

#include <cmath>
#include <utility>

#if 0
struct csv_whitespace : std::ctype<char>
{
    static const mask* make_table()
    {
        static std::vector<mask> v(classic_table(),
                                   classic_table() + table_size);

        v[','] |= space;  // comma will be classified as whitespace
        v[' '] &= ~space; // space will not be classified as whitespace

        return &v[0];
    }

    csv_whitespace(std::size_t refs = 0)
      : ctype(make_table(), false, refs)
    {}
};
#endif

struct csv_whitespace : std::ctype<char>
{
    static mask* make_table()
    {
        auto* v = new mask[table_size];
        std::copy_n(classic_table(), table_size, v);

        v[static_cast<int>(',')] |=
          space; // comma will be classified as whitespace
        v[static_cast<int>(' ')] &=
          ~space; // space will not be classified as whitespace

        return &v[0];
    }

    csv_whitespace(std::size_t refs = 0)
      : ctype(make_table(), true, refs)
    {}
};

enum class read_status
{
    success,
    eof,
    error
};

static bool
have_lp_extension(std::string filename)
{
    auto pos = filename.rfind('.');

    if (pos == std::string::npos)
        return false;

    return filename.compare(pos, 3, ".lp");
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

struct bench
{
    struct model
    {
        model() = default;

        model(std::string name_)
          : name(name_)
          , best_solution(std::numeric_limits<double>::infinity())
        {}

        model(std::string name_, double best_solution_)
          : name(name_)
          , best_solution(best_solution_)
        {}

        std::string name;     // e.g. scpnrf3, scpd1 etc.
        double best_solution; // e.g. the best or optimum solution.

        std::string filename() const
        {
            return have_lp_extension(name) ? name : fmt::format("{}.lp", name);
        }

        std::size_t name_length() const noexcept
        {
            return name.size();
        }

        bool operator<(const model& other)
        {
            return name < other.name;
        }

        bool operator==(const model& other)
        {
            return name == other.name;
        }
    };

    struct solver
    {
        solver() = default;

        solver(std::string name_)
          : name(std::move(name_))
        {}

        std::string name; // e.g. cplex-10.0.1, baryonyx-0.2

        bool operator<(const solver& other)
        {
            return name < other.name;
        }

        bool operator==(const solver& other)
        {
            return name == other.name;
        }
    };

    struct element
    {
        element() = default;

        element(double duration_, double solution_)
          : duration(duration_)
          , solution(solution_)
          , success(true)
        {}

        element(double duration_)
          : duration(duration_)
        {}

        double duration = std::numeric_limits<double>::infinity();
        double solution = std::numeric_limits<double>::infinity();
        bool success = false;
    };

    std::vector<model> models;
    std::vector<solver> solvers;
    std::unordered_map<std::string, int> filename_id;
    std::unordered_map<std::string, int> solver_id;
    baryonyx::dynamic_array<double> array;

    int push_back_model(std::string name, double best_solution)
    {
        fmt::print("push back model {} - {}\n", name, best_solution);

        auto ret = filename_id.emplace(name, 0);
        if (ret.second) {
            fmt::print(" => does not exists.\n");
            models.emplace_back(name, best_solution);
            ret.first->second = static_cast<int>(models.size());
        }

        return ret.first->second;
    }

    int push_back_solver(std::string name)
    {
        fmt::print("push back solver {}\n", name);
        auto ret = solver_id.emplace(name, 0);
        if (ret.second) {
            fmt::print(" => does not exists.\n");
            solvers.emplace_back(name);
            ret.first->second = static_cast<int>(solvers.size());
        }

        return ret.first->second;
    }

    void save(std::ostream& os,
              std::string name,
              const std::vector<element>& c)
    {
        fmt::print(os, "file,opt");

        for (const auto& elem : solvers)
            fmt::print(os, ",{}", elem.name);
        fmt::print(os, ",{}\n", name);

        for (std::size_t i{ 0 }, e{ models.size() }; i != e; ++i) {
            fmt::print(os, "{},{},", models[i].name, models[i].best_solution);

            for (std::size_t j{ 0 }, end_j{ solvers.size() }; j != end_j; ++j)
                fmt::print(os, "{},", array(i, j));

            fmt::print(os, "{}\n", c[i].solution);
        }
    }

    bool load(std::istream& is)
    {
        std::string line;
        line.reserve(BUFSIZ);

        if (!std::getline(is, line)) {
          fmt::print(stderr, fmt::emphasis::bold | fmt::fg(fmt::terminal_color::red),
                       "benchmark: fail to read header\n");
            return false;
        }

        int line_pos = 1;
        {
            std::string::size_type begin{ 0 };
            std::string::size_type end = line.find(',', begin);
            int i{ 0 };

            do {
                ++i;

                fmt::print("load header: {}\n",
                           line.substr(begin, end - begin));

                if (i > 2)
                    push_back_solver(line.substr(begin, end - begin));

                begin = end == std::string::npos ? end : end + 1;
                end = line.find(',', begin);
            } while (begin < line.size());

            fmt::print("{} solvers\n", solvers.size());
        }

        ++line_pos;
        array.init(solvers.size());

        is.imbue(std::locale(is.getloc(), new csv_whitespace));

        {
            std::string name;
            std::string best_solution;
            std::string current;

            while (is) {
                is >> name >> best_solution;
                if (!is.good())
                    break;

                fmt::print("load array: {} {}\n", name, best_solution);

                push_back_model(
                  name,
                  to_double(best_solution,
                            std::numeric_limits<double>::infinity()));

                array.push_line();

                auto length{ solvers.size() };
                for (std::size_t i{ 0 }; i != length; ++i) {
                    fmt::print("rest solver\n");
                    is >> current;
                    if (!is) {
                        fmt::print(fmt::emphasis::bold |
                                     fmt::fg(fmt::terminal_color::red),
                                   "benchmark: fail to read data at line {}\n",
                                   line_pos);
                        return false;
                    }

                    array(array.rows() - 1, i) = to_double(
                      current, std::numeric_limits<double>::infinity());
                }

                ++line_pos;
            }
        }

        return !models.empty();
    }

    void show_to_console(const std::vector<element>& current, std::string name)
    {
        const auto model_size{ models.size() };
        const auto solver_size{ solvers.size() };

        std::vector<int> row_length(solvers.size(), 5);
        int best_solution_row_length = { 10 };
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

            if (is_valid(models[i].best_solution)) {
                auto digits = get_digits_number(models[i].best_solution);

                if (digits > best_solution_row_length)
                    best_solution_row_length = digits;
            }
        }

        auto modelname_lenght_max = std::max(
          static_cast<std::string::size_type>(4),
          std::max_element(models.cbegin(),
                           models.cend(),
                           [](const auto& lhs, const auto& rhs) {
                               return lhs.name.size() < rhs.name.size();
                           })
            ->name.size());

        fmt::print("{:^{}} {:>{}} | ",
                   "file",
                   modelname_lenght_max,
                   "best",
                   best_solution_row_length);

        for (std::size_t solver{ 0 }; solver != solver_size; ++solver)
            fmt::print("{:>{}} ",
                       solvers[solver].name.substr(0, 10),
                       row_length[solver]);

        fmt::print("{:>{}}\n", name, current_row_length);

        for (size_t i{ 0 }; i != model_size; ++i) {
            fmt::print("{:^{}} {:>{}.3f} | ",
                       models[i].name,
                       modelname_lenght_max,
                       models[i].best_solution,
                       best_solution_row_length);

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

bool
benchmark(const baryonyx::context_ptr& ctx,
          std::string filepath,
          std::string name)
{
    std::ifstream ifs(filepath);
    if (!ifs.is_open()) {
        fmt::print(stderr,
                   fmt::emphasis::bold | fmt::fg(fmt::terminal_color::red),
                   "Can not open {} to start benchmark\n",
                   filepath);
        return false;
    }

    bench b;
    if (!b.load(ifs)) {
        fmt::print(stderr,
                   fmt::emphasis::bold | fmt::fg(fmt::terminal_color::red),
                   "Error reading the benchmark file: {}\n",
                   filepath);
        return false;
    }

    std::vector<bench::element> current(b.models.size());

    std::string dirname = filepath.substr(0, filepath.rfind(".csv")) + '/';

    for (std::size_t i{ 0 }, e{ b.models.size() }; i != e; ++i) {
        fmt::print("- optimizing {}\n", b.models[i].name);

        auto filename = b.models[i].filename();

        try {
            fmt::print("- optimizing: {} filename: {} in dirname: {}\n",
                       b.models[i].name,
                       filename,
                       dirname);

            auto rawpb = baryonyx::make_problem(ctx, dirname + filename);
            auto result = baryonyx::optimize(ctx, rawpb);

            if (result)
                current[i] = { result.duration,
                               result.solutions.back().value };
            else
                current[i] = { result.duration };

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

        for (std::size_t i{ 0 }, e{ current.size() }; i != e; ++i) {
            if (is_valid(b.models[i].best_solution) &&
                is_valid(current[i].solution)) {
                ++optimum_number;

                mean_distance +=
                  ((b.models[i].best_solution - current[i].solution) /
                   b.models[i].best_solution) *
                  100.0;
            }
        }

        if (optimum_number > 0)
            fmt::print("Average distance from the optimum: {}%\n",
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
                fmt::print("Average distance from the solution of {}: {}%\n",
                           b.solvers[solver].name,
                           mean_distance / current.size());
        }
    }

    std::ofstream ofs(filepath);
    if (ofs.is_open())
        b.save(ofs, name, current);

    b.show_to_console(current, name);

    return true;
}
