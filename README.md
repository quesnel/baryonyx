# Baryonyx-0.1

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

## Compilation

Default, `baryonyx` provides a shared library `libbaryonyx-0.1.so`, a
static library `libbaryonyx-0.1.a` and an executable
`baryonyx-0.1`. To compile and install:

    cd baryonyx
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
    make -j 4
    make install

For example, to use `clang-4.0`, overrides the `CC` and `CXX`
variables as follow:

    cd baryonyx
    mkdir build
    cd build
    export CC=clang-4.0
    export CXX='clang++-4.0 -stdlib=libc++`
    cmake -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
    make -j 4
    make install

To override the default build flags, remove the `CMAKE_BUILD_TYPE`
parameter and override the `CXXFLAGS` as follow:

    export CC=gcc-7
    export CXX=g++-7
    export CXXFLAGS='-Wall -Wextra -Werror -ftree-vectorize -mmmx -msse -msse2 -msse3 -O3'
    cmake -DCMAKE_INSTALL_PREFIX=/usr ..

## Usage

    baryonyx-0.1 file.lp

Previous command line install `baryonyx` program and library into the
`/usr` prefix. If you install into another directory, you need to
define three environment variables (into your `.bashrc` or
equivalent). If you install into `$HOME/usr` for example, you need to
define:

    export PATH=$PATH:$HOME/usr
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/usr/lib
    export PKG_CONFIG_PATH=$PKG_CONFIX_PATH:$HOME/usr/lib/pkgconfig

## R

To use `rbaryonyx`, you must compile and install the `baryonyx`
library. Follow the previous section. The R `rbaryonyx` package
requires several packages. Under the R terminal:

    cd baryonyx/rbaryonyx

    # Removes previous installed version of rbaryonyx
    R CMD REMOVE rbaryonyx

    # Install the dependencies of rbaryonyx
    install.packages("roxygen2")
    install.packages("Rcpp")
    install.packages("rtools")

    library(Rcpp)
    compileAttributes(".")
    library(devtools)
    devtools::document()
    devtools::build()
    devtools::install()

    library(rbaryonyx)
    ?rbaryonyx
