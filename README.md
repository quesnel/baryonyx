baryonyx-solver
===============

Copyright © 2014-2016 INRA

The software is released under the MIT license. See the COPYING file.

## Requirements

* cmake (≥ 3)
* c++ compiler with c++14 support (gcc ≥ 4.9, clang ≥ 3.5)

For recent Debian and Ubuntu derivatives (remove clang to only use
gcc):

    apt-get install build-essential cmake clang

## Compilation

Compiling and installing:

    cd baryonyx
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
    make -j 4
    make install

To use clang replace the previous `./configure ...` command with the following:

    cd baryonyx
    mkdir build
    cd build
    export CC=clang
    export CXX=clang++-libc++
    cmake -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
    make -j 4
    make install

# Usage

    $ baryonyx file.lp
