from typing import Callable, Literal

ConstraintOrder = Literal[
    "none", "reversing", "random-sorting", "infeasibility-decr",
    "infeasibility-incr", "lagrangian-decr", "lagrangian-incr",
    "pi-sign-change", "cycle"]
Norm = Literal["none", "random", "l1", "l2", "loo"]
Optimizer = Literal["genetic", "restart", "tuned", "lns"]
Progress = Callable[[int, float, int, float], None]

__version__: str

class Result:
    status: str
    success: bool
    minimize: bool
    value: float | None
    duration: float
    loop: int
    variables: int
    constraints: int
    remaining_constraints: int
    values: list[float]
    solution: dict[str, int] | None
    def __bool__(self) -> bool: ...

def solve(file_path: str, *, limit: int = 1000, time_limit: float = -1.0,
          thread: int = 0, seed: int = 0, theta: float = 0.5,
          delta: float = -1.0, kappa_min: float = 0.0,
          kappa_step: float = 1e-3, alpha: float = 1.0,
          constraint_order: ConstraintOrder = "none", norm: Norm = "l2",
          optimizer: Optimizer = "lns", preprocessing: bool = True,
          verbose: bool | int = False,
          progress: Progress | None = None) -> Result: ...

def optimize(file_path: str, *, limit: int = 1000, time_limit: float = -1.0,
             thread: int = 0, seed: int = 0, theta: float = 0.5,
             delta: float = -1.0, kappa_min: float = 0.0,
             kappa_step: float = 1e-3, alpha: float = 1.0,
             constraint_order: ConstraintOrder = "none", norm: Norm = "l2",
             optimizer: Optimizer = "lns", preprocessing: bool = True,
             verbose: bool | int = False,
             progress: Progress | None = None) -> Result: ...

def request_stop() -> None: ...
