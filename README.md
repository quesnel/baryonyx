# Baryonyx-0.2

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
* `bool` **serialize**: true to store for each loop, constraint and variable states.

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

# Usage

Apply morris method to found parameters:

````
library(rbaryonyx)
library(sensitivity)

factors=c("theta", "delta", "constraint_order", "kappa_min",
                   "kappa_step", "kappa_max", "alpha", "w")
bounds = data.frame(min=c(0,   0, 0, 0.0,    0,  1.0, 1.0, 50),
                    max=c(1, 0.1, 4, 0.25, 0.1,  1.0, 1.0, 50))
rownames(bounds) <- factors

morrisDesign <- morris(model = NULL,
                factors = factors,
                r = 100,
                design=list(type="oat", levels=20, grid.jump=5),
                binf = bounds$min,
                bsup = bounds$max,
                scale=TRUE)

# morrisDesign$X: # 20 * (8 + 1)

solve_lp <- function(x, file_path, thread=1, limit=10000, time_limit=10) {
  theta <- x["theta"]
  delta <- x["delta"]
  constraint_order <- x["constraint_order"]
  kappa_min <- x["kappa_min"]
  kappa_step <- x["kappa_step"]
  kappa_max <- x["kappa_max"]
  alpha <- x["alpha"]
  w <- x["w"]

  r <- rbaryonyx::optimize_01lp_problem(file_path = file_path,
           limit = limit,
           theta = theta,
           delta = delta,
           constraint_order = constraint_order,
           kappa_min = kappa_min,
           kappa_step = kappa_step,
           kappa_max = kappa_max,
           alpha = alpha,
           w = w,
           time_limit = time_limit,
           seed = 123654785,
           thread = 1L,
           verbose = FALSE)

  return(r)
}

r = apply(morrisDesign$X, 1, solve_lp, file_path="../lib/test/prevl1.lp")
mr = matrix(unlist(r), ncol=2, byrow = TRUE)
mr1 <- mr[,-3]

tell(morrisDesign, mr1)
plot(morrisDesign)
````

Apply rgenoud to found best parameters:

````
library(rbaryonyx)
library(rgenoud)
library(rilp)
library(parallel)

optim_gen_lp <- function(x) {
  r <- rilp::optimize_01lp_problem(
           file_path="../lib/test/prevl1.lp",
           limit = 4000,
           theta = x[1],
           delta = x[2],
           constraint_order = 1,
           kappa_min = x[3],
           kappa_step = x[4],
           kappa_max = 1.0,
           alpha = 1.0,
           w = 60,
           time_limit = 30,
           seed = 123654785,
           thread = 1L,
           verbose = FALSE)

  return(r)
}

d = matrix(c(0.0, 1e-10, 0.0, 1e-9,
             0.9,  1e-5, 0.9, 1e-7),
             nrow=4, ncol=2)

s = c(0.1, 1e-5, 0.2, 1e-8)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, outfile="debug.txt")

claw1 <- genoud(optim_gen_lp, nvars=4,
                Domains=d,
                starting.values=s,
                cluster=cl,
                max=FALSE, pop.size=10)
````


# Upgrade

To upgrade to the latest version of rbaryonyx, under bash (or equivalent):

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

