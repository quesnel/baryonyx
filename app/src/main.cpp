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

#include <baryonyx/core-out>
#include <baryonyx/core>

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

#include <cmath>
#include <cstring>
#include <utility>

#ifndef __WIN32
#include <unistd.h>
#endif

namespace {

struct write_parameters
{
    std::shared_ptr<baryonyx::context> ctx;

    write_parameters(std::shared_ptr<baryonyx::context> ctx_)
      : ctx(std::move(ctx_))
    {
    }
};

} // anonymous namespace

static const char* file_format_error_format(
  baryonyx::file_format_error_tag) noexcept;
static const char* problem_definition_error_format(
  baryonyx::problem_definition_error_tag) noexcept;
static const char* solver_error_format(baryonyx::solver_error_tag) noexcept;
static std::ostream&
operator<<(std::ostream& os, const write_parameters& wp);
static baryonyx::result
solve_or_optimize(std::shared_ptr<baryonyx::context> ctx,
                  baryonyx::problem& pb);

int
main(int argc, char* argv[])
{
    auto ctx = std::make_shared<baryonyx::context>();
    ctx->set_standard_stream_logger();
    int i = ctx->parse(argc, argv);

    if (i < 0)
        return EXIT_FAILURE;

    for (; i < argc; ++i) {
        try {
            auto pb = baryonyx::make_problem(ctx, argv[i]);

            std::string filename(argv[i]);
            filename += '-';
            filename += std::to_string(::getpid());
            filename += '_';
            filename += std::to_string(i);
            filename += ".sol";

            ctx->info("Output file: %s\n", filename.c_str());

            if (ctx->check()) {
                auto filename =
                  ctx->get_string_parameter("check-filename", std::string());

                auto result = baryonyx::make_result(ctx, filename);
                if (result) {
                    auto valid = baryonyx::is_valid_solution(pb, result);
                    ctx->info("Check %s with %s: %s",
                              argv[i],
                              filename.c_str(),
                              (valid ? "success" : "failure"));
                }
            }

            std::ofstream ofs(filename);
            ofs << std::boolalpha
                << std::setprecision(std::floor(
                     std::numeric_limits<double>::digits * std::log10(2) + 2));

            auto now = std::chrono::system_clock::now();
            auto in_time_t = std::chrono::system_clock::to_time_t(now);

            ofs << baryonyx::resume(pb) << R"(\ solver starts: )"
                << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X")
                << "\n\\ parameters:\n"
                << ::write_parameters(ctx);

            auto ret = solve_or_optimize(ctx, pb);

            in_time_t = std::chrono::system_clock::to_time_t(now);
            ofs << R"(\ solver finishes: )"
                << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X")
                << '\n';

            if (ret.status == baryonyx::result_status::success) {
                ofs << R"(\ Solution found: )" << ret.value << '\n' << ret;
            } else {
                ofs << R"(\ Solution not found. Missing constraints: )"
                    << ret.remaining_constraints << '\n';
            }
        } catch (const baryonyx::precondition_failure& e) {
            ctx->error("internal failure\n");
        } catch (const baryonyx::postcondition_failure& e) {
            ctx->error("internal failure\n");
        } catch (const baryonyx::numeric_cast_failure& e) {
            ctx->error("numeric cast interal failure\n");
        } catch (const baryonyx::file_access_failure& e) {
            ctx->error("file `%s' fail %d: %s\n",
                       e.file().c_str(),
                       e.error(),
                       std::strerror(e.error()));
        } catch (const baryonyx::file_format_failure& e) {
            ctx->error("file format error at line %d column %d "
                       "%s\n",
                       e.line(),
                       e.column(),
                       file_format_error_format(e.failure()));
        } catch (const baryonyx::problem_definition_failure& e) {
            ctx->error("definition problem error at %s: %s\n",
                       e.element().c_str(),
                       problem_definition_error_format(e.failure()));
        } catch (const baryonyx::solver_failure& e) {
            ctx->error("solver error: %s\n", solver_error_format(e.failure()));
        } catch (const std::exception& e) {
            ctx->error("failure: %s.\n", e.what());
        }
    }

    return EXIT_SUCCESS;
}

static const char*
file_format_error_format(baryonyx::file_format_error_tag failure) noexcept
{
    static const char* const tag[] = {
        "end of file",     "unknown",
        "already defined", "incomplete",
        "bad name",        "bad operator",
        "bad integer",     "bad objective function type",
        "bad bound",       "bad function element",
        "bad constraint"
    };

    switch (failure) {
    case baryonyx::file_format_error_tag::end_of_file:
        return tag[0];
    case baryonyx::file_format_error_tag::unknown:
        return tag[1];
    case baryonyx::file_format_error_tag::already_defined:
        return tag[2];
    case baryonyx::file_format_error_tag::incomplete:
        return tag[3];
    case baryonyx::file_format_error_tag::bad_name:
        return tag[4];
    case baryonyx::file_format_error_tag::bad_operator:
        return tag[5];
    case baryonyx::file_format_error_tag::bad_integer:
        return tag[6];
    case baryonyx::file_format_error_tag::bad_objective_function_type:
        return tag[7];
    case baryonyx::file_format_error_tag::bad_bound:
        return tag[8];
    case baryonyx::file_format_error_tag::bad_function_element:
        return tag[9];
    case baryonyx::file_format_error_tag::bad_constraint:
        return tag[10];
    }

    return nullptr;
}

static const char*
problem_definition_error_format(
  baryonyx::problem_definition_error_tag failure) noexcept
{
    static const char* const tag[] = {
        "empty variables",
        "empty objective function",
        "variable not used",
        "bad bound",
        "multiple constraints with different value"
    };

    switch (failure) {
    case baryonyx::problem_definition_error_tag::empty_variables:
        return tag[0];
    case baryonyx::problem_definition_error_tag::empty_objective_function:
        return tag[1];
    case baryonyx::problem_definition_error_tag::variable_not_used:
        return tag[2];
    case baryonyx::problem_definition_error_tag::bad_bound:
        return tag[3];
    case baryonyx::problem_definition_error_tag::multiple_constraint:
        return tag[4];
    }

    return nullptr;
}

static const char*
solver_error_format(baryonyx::solver_error_tag failure) noexcept
{
    static const char* const tag[] = { "no solver available",
                                       "unrealisable constraint",
                                       "not enough memory" };

    switch (failure) {
    case baryonyx::solver_error_tag::no_solver_available:
        return tag[0];
    case baryonyx::solver_error_tag::unrealisable_constraint:
        return tag[1];
    case baryonyx::solver_error_tag::not_enough_memory:
        return tag[2];
    }

    return nullptr;
}

static std::ostream&
operator<<(std::ostream& os, const write_parameters& wp)
{
    if (not wp.ctx)
        return os;

    const auto& params = wp.ctx->get_parameters();
    for (const auto & param : params) {
        os << R"(\ )" << param.first << " = ";

        switch (param.second.type) {
        case baryonyx::parameter::tag::string:
            os << param.second.s;
            break;
        case baryonyx::parameter::tag::integer:
            os << param.second.l;
            break;
        case baryonyx::parameter::tag::real:
            os << param.second.d;
            break;
        }

        os << '\n';
    }

    return os;
}

static baryonyx::result
solve_or_optimize(std::shared_ptr<baryonyx::context> ctx,
                  baryonyx::problem& pb)
{
    if (ctx->optimize())
        return baryonyx::optimize(ctx, pb);

    return baryonyx::solve(ctx, pb);
}
