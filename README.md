# Baryonyx-0.2.1

[![Build Status](https://travis-ci.org/quesnel/baryonyx.png?branch=master)](https://travis-ci.org/quesnel/baryonyx)

[Baryonyx](https://en.wikipedia.org/wiki/Baryonyx) is an integer and
binary linear programming solver based on
the [Dag Wedelin](http://www.cse.chalmers.se/~dag/) heuristic.

Copyright © 2017 [INRA](http://www.inra.fr/en)

Author [Gauthier Quesnel](https://mia.toulouse.inra.fr/Gauthier_QUESNEL) <gauthier.quesnel@inra.fr>

The software is released under the MIT license. See the LICENSE file.

## Requirements

* `cmake` (≥ 3)
* c++ compiler with c++14 support (`gcc` ≥ 5, `clang` ≥ 3.5)

For recent Debian and Ubuntu derivatives (remove clang to only use
gcc):

    apt-get install build-essential cmake clang

## Installation

Default, baryonyx provides a shared library `libbaryonyx-0.2.so`, a static
library `libbaryonyx-0.2.a` and an executable `baryonyx-0.2`. To compile and
install in default cmake's install directory:

    cd baryonyx
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make install

Previous command line install baryonyx program and library into the
`/usr/local` prefix. If you install into another directory, you need to define
three environment variables (into your `.bashrc` or equivalent). If you install
into `$HOME/usr` for example, you need to define in your `.bashrc`:

    export PATH=$PATH:$HOME/usr
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/lib
    export PKG_CONFIG_PATH=$PKG_CONFIX_PATH:$HOME/usr/lib/pkgconfig

Then run the following commands:

    cd baryonyx
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/usr ..
    make install

To override the default build flags, remove the `CMAKE_BUILD_TYPE`
parameter and override the `CXXFLAGS` as follow:

    export CC=gcc-7
    export CXX=g++-7
    export CXXFLAGS='-Wall -Wextra -Werror -ftree-vectorize -mmmx -msse -msse2 -msse3 -O3'
    cmake -DCMAKE_INSTALL_PREFIX=$HOME/usr ..

## Usage

To run baryonyx in solver mode (ie. trying to valid all constraints):

    baryonyx-0.2 file.lp

To run baryonyx into the heuristic model (ie. trying to valid all constraints
and optimize the solution), add a `-o` or `--optimize` option to the command
line:

    baryonyx-0.2 file.lp

The baryonyx solver have many parameters. To assign it, use the `-p
[name]:value` syntax in the command line:

* `real` **time_limit**: time in second to stop the solver or the heuristic
* `real` **theta**:
* `real` **delta**:
* `real` **kappa_min**:
* `real` **kappa_step**:
* `real` **kappa_max**:
* `real` **alpha**:
* `real` **pushing_k_factor**:
* `integer` **pushes_limit**:
* `integer` **pushing_objective_amplifier**:
* `integer` **pushing_iteration_limit**:
* `integer` **limit**: number of loop to stop the solver or the heuristic
* `integer` **w**:
* `constraint_order` **order**: `none`, `reversing`, `random-sorting`, `infeasibility-decr`, `infeasibility-incr`.
* `string` **preprocessing**: `none`,  `variables-number`, `variables-weight`m `constraints-weight`, `implied`.
* `string` **norm**: `none`, `l1`, `l2`, `rng`, `infinity`
* `integer` **print-level**: show debug information if greater than 0.

For example:

    baryonyx -p limit:1000000 lib/test/prevl1.lp
    baryonyx -p limit: -p kappa-min:0.2 lib/test/prevl1.lp

## Upgrade

To upgrade to the latest version of baryonyx, under bash (or equivalent):

    cd baryonyx
    git pull -r
    cd build
    make install

# R

To use rbaryonyx, you must compile and install the baryonyx library. Follow
the previous section and install R.

## Installation

To install rbaryonyx, first, remove old installation under bash (or
equivalent):

    cd baryonyx/rbaryonyx
    # Removes previous installed version of rbaryonyx
    R CMD REMOVE rbaryonyx

The R rbaryonyx package requires several packages. Then, under a R terminal:

    # Install the dependencies of rbaryonyx
    install.packages("roxygen2")
    install.packages("Rcpp")
    install.packages("devtools")

    library(Rcpp)
    compileAttributes(".")
    library(devtools)
    devtools::document()
    devtools::build()
    devtools::install()

    library(rbaryonyx)
    ?rbaryonyx

## API

Two functions are provided to solve or optimize 01 linear programming problem. Parameters are the same as `C++ API`. These function returns a scalar:

- If a solution is found:
  - if the problem is a minimization: the value of the solution found.
  - if the problem is a maximization: the inverse of the solution found.
- If no solution is found, we use the limits of the objective
  function (minimal and maximal value possible.
  - if the problem is a minimization: the maximal value possible + the remaining constraints.
  - if the problem is a maximization: the inverse of the minimal value possible + the remaining constraints.
- If a error occurred (not enough memory, problem error etc.):
  - if the problem is a minimization: the maximal value possible + the number of constraints .
  - if the problem is a maximization: the inverse of the minimal value possible + the number of constraints.

```R
solve_01lp_problem <- function(file_path, limit = 1000L, theta = 0.5,
  delta = 1e-4, constraint_order = 0L, kappa_min = 0.1, kappa_step = 1e-4,
  kappa_max = 1.0, alpha = 1.0, w = 500L, time_limit = 10.0, seed = -1L,
  thread = 1L, norm = 4L, pushing_k_factor = 0.9,
  pushing_objective_amplifier = 5.0, pushes_limit = 10L,
  pushing_iteration_limit = 20L, float_type = 1L, verbose = TRUE)

optimize_01lp_problem <- function(file_path, limit = 1000L, theta = 0.5,
  delta = 1e-4, constraint_order = 0L, kappa_min = 0.1, kappa_step = 1e-4,
  kappa_max = 1.0, alpha = 1.0, w = 500L, time_limit = 10.0, seed = -1L,
  thread = 1L, norm = 4L, pushing_k_factor = 0.9,
  pushing_objective_amplifier = 5.0, pushes_limit = 10L,
  pushing_iteration_limit = 20L, float_type = 1L, verbose = TRUE)

```

## Usage

Apply morris method to found useful parameters:

```R
library(rbaryonyx)
library(sensitivity)

factors = c("theta", "delta", "constraint_order", "kappa_min", "kappa_step",
  "kappa_max", "alpha", "w", "norm", "pushing_k_factor",
  "pushing_objective_amplifier", "pushes_limit", "pushing_iteration_limit",
  "float_type")

bounds = data.frame(
  min=c(
    0,     # theta
    0,     # delta
    0,     # constraint_order
    0,     # kappa_min
    1e-16, # kappa_step
    1.0,   # kappa_max
    0.0,   # alpha
    50,    # w
    0,     # norm
    0.1,   # pushing_k_factor
    1.0,   # pushing_objective_amplifier
    10,    # pushes_limit
    20,    # pushing_iteration_limit
    0
    ),    # float_type
max=c(
    1,     # theta
    0,     # delta
    4,     # constraint_order
    0.1,   # kappa_min
    1e-1,  # kappa_step
    1.0,   # kappa_max
    2.0,   # alpha
    500,   # w
    4,     # norm
    1,     # pushing_k_factor
    10.0,  # pushing_objective_amplifier
    100,   # pushes_limit
    200,   # pushing_iteration_limit
    2))    # float_type

rownames(bounds) <- factors

morrisDesign <- morris(model = NULL,
                factors = factors,
                r = 10,
                design=list(type="oat", levels=10, grid.jump=5),
                binf = bounds$min,
                bsup = bounds$max,
                scale=TRUE)

solve_lp <- function(x, file_path, limit=10000, time_limit=10, seed=123456789, thread=1) {
  r <- rbaryonyx::solve_01lp_problem(file_path = file_path,
                   limit = limit,
                   theta = x["theta"],
                   delta = x["delta"],
                   constraint_order = x["constraint_order"],
                   kappa_min = x["kappa_min"],
                   kappa_step = x["kappa_step"],
                   kappa_max = x["kappa_max"],
                   alpha = x["alpha"],
                   w = x["w"],
                   time_limit = time_limit,
                   seed = seed,
                   thread = thread,
                   norm = x["norm"],
                   pushing_k_factor = x["pushing_k_factor"],
                   pushing_objective_amplifier = x["pushing_objective_amplifier,"],
                   pushes_limit = x["pushes_limit"],
                   pushing_iteration_limit = x["pushing_iteration_limit"],
                   float_type = x["float_type"])

  return(r)
}

r = apply(morrisDesign$X, 1, solve_lp, file_path="verger_5_5.lp", thread=1, limit=10000, time_limit=10, seed=123456789)

morrisDesign$Y <- r
mu <- apply(morrisDesign$X,2,mean)
mu.star <- apply(morrisDesign$X, 2, function(x) mean(abs(x)))
sigma <- apply(morrisDesign$ee, 2, sd)

apply(morrisDesign$X, 2, function(v) plot(factor(v), r))
```

# Upgrade

To upgrade to the latest version of rbaryonyx, under bash (or equivalent):

```bash
cd baryonyx
git pull -r
cd build
make install
R CMD REMOVE rbaryonyx
cd rbaryonyx
Rscript -e 'library(Rcpp); compileAttributes(".")'
Rscript -e 'library(devtools); devtools::document()'
cd ..
R CMD build rbaryonyx
R CMD INSTALL rbaryonyx_1.0.tar.gz
```
