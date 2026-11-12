import os
import re
import threading
import time

import pytest

import baryonyx as bx

DATA = os.path.join(os.path.dirname(__file__), "..", "..", "lib", "test")


def path(name):
    return os.path.join(DATA, name)


TERM = re.compile(r"([+-])?\s*(\d+(?:\.\d*)?)?\s*([A-Za-z_]\w*)")


def terms(text):
    return [((-1.0 if sign == "-" else 1.0) * float(factor or 1), name)
            for sign, factor, name in TERM.findall(text)]


def read_lp(name):
    """Objective and constraints of a linear LP file of the tests."""
    text = open(path(name)).read()
    sense, rest = text.split(None, 1)
    objective, rest = re.split(r"subject to:?", rest, maxsplit=1,
                               flags=re.IGNORECASE)
    body = re.split(r"\n(?:bounds?|binary|binaries|general|end)\b", rest,
                    flags=re.IGNORECASE)[0]
    constraints = []
    # A constraint may span several lines: split on the right-hand sides.
    for m in re.finditer(r"(.*?)(<=|>=|=)\s*(-?\d+)", body, re.DOTALL):
        lhs = re.sub(r"^\s*\w+:", "", m.group(1).strip())
        constraints.append((terms(lhs), m.group(2), int(m.group(3))))
    return sense.lower(), terms(objective), constraints


def check(name, r):
    sense, objective, constraints = read_lp(name)
    x = r.solution
    assert x is not None
    assert r.minimize == (sense == "minimize")
    for lhs, op, rhs in constraints:
        v = sum(f * x[n] for f, n in lhs)
        assert {"=": v == rhs, "<=": v <= rhs, ">=": v >= rhs}[op], (lhs, op, rhs)
    assert sum(f * x[n] for f, n in objective) == pytest.approx(r.value)


def test_optimize_assignment():
    r = bx.optimize(path("assignment_problem_1.lp"), time_limit=1, thread=2,
                    seed=1)
    assert r and r.success and r.status == "success"
    assert len(r.solution) == 16
    assert r.remaining_constraints == 0
    assert r.values[-1] == r.value
    check("assignment_problem_1.lp", r)


@pytest.mark.parametrize("optimizer", ["genetic", "restart", "tuned", "lns"])
def test_optimizers(optimizer):
    r = bx.optimize(path("assignment_problem_1.lp"), time_limit=0.5,
                    thread=2, seed=1, optimizer=optimizer, norm="loo",
                    constraint_order="random-sorting")
    assert r
    check("assignment_problem_1.lp", r)


def test_solve():
    r = bx.solve(path("8_queens_puzzle.lp"), limit=10000, seed=1)
    assert r
    check("8_queens_puzzle.lp", r)


def test_preprocessor_solution():
    r = bx.solve(path("prepro.lp"))
    assert r
    check("prepro.lp", r)
    assert r.value == pytest.approx(6.5)


def test_without_preprocessing():
    r = bx.optimize(path("assignment_problem_1.lp"), time_limit=0.5, seed=1,
                    preprocessing=False)
    assert r
    check("assignment_problem_1.lp", r)


def test_errors():
    with pytest.raises(ValueError, match="norm.*loo"):
        bx.solve(path("assignment_problem_1.lp"), norm="l3")
    with pytest.raises(FileNotFoundError):
        bx.solve(path("missing.lp"))
    with pytest.raises(RuntimeError, match="not binary variable"):
        bx.solve(path("sudoku.lp"))
    with pytest.raises(TypeError):
        bx.solve(path("assignment_problem_1.lp"), 1000)


def test_progress():
    calls = []
    r = bx.optimize(path("capmo1_direct.lp"), time_limit=2, thread=2, seed=1,
                    progress=lambda *args: calls.append(args))
    assert r and calls
    remaining, value, loop, duration = calls[-1]
    assert remaining == 0 and duration >= 0


def test_progress_exception():
    def fail(*args):
        raise KeyError("stop")

    start = time.monotonic()
    with pytest.raises(KeyError):
        bx.optimize(path("capmo1_direct.lp"), time_limit=30, thread=2,
                    seed=1, progress=fail)
    assert time.monotonic() - start < 10


def test_request_stop():
    # Without time limit, optimize only returns on request; the GIL is
    # released, so another Python thread can ask it.
    timer = threading.Timer(1.0, bx.request_stop)
    timer.start()
    start = time.monotonic()
    r = bx.optimize(path("capmo1_direct.lp"), time_limit=0, thread=2, seed=1)
    timer.join()
    assert time.monotonic() - start < 10
    assert r
    check("capmo1_direct.lp", r)
