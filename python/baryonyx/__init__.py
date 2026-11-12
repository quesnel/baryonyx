"""Baryonyx: 0-1 integer linear programming with the Wedelin heuristic.

    >>> import baryonyx as bx
    >>> r = bx.optimize("problem.lp", time_limit=10, thread=4, seed=1)
    >>> if r:
    ...     print(r.value, r.solution)

``solve`` searches one solution, ``optimize`` the best one within the time
limit. The parameters have the names of the command line ``-p`` parameters
with ``_`` instead of ``-``; the enumerations take the same strings:

- ``constraint_order``: none, reversing, random-sorting, infeasibility-decr,
  infeasibility-incr, lagrangian-decr, lagrangian-incr, pi-sign-change,
  cycle;
- ``norm``: none, random, l1, l2, loo;
- ``optimizer``: genetic, restart, tuned, lns.

``time_limit`` <= 0 means no limit: ``optimize`` then runs until Ctrl-C or
``request_stop()``, and returns its best result. ``progress`` is called with
``(remaining_constraints, value, loop, duration)`` at each improvement.
"""

from ._baryonyx import Result, optimize, request_stop, solve

__all__ = ["Result", "optimize", "request_stop", "solve"]
__version__ = "0.6.0"
