lp
==

Copyright © 2014-2016 INRA

The software is released under the MIT license. See the COPYING file.

## Requirements

* boost (≥ 1.50)
* eigen3 (≥ 3)
* cmake (≥ 3)
* c++ compiler with c++14 support (gcc ≥ 4.9, clang ≥ 3.5)

For recent Debian and Ubuntu derivatives (remove clang to only use
gcc):

    apt-get install build-essential cmake clang \
                    libboost-dev libeigen3-dev

## Compilation

Compiling and installing:

    cd lp
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
    make -j 4
    make install

To use clang replace the previous `./configure ...` command with the following:

    cd lp
    mkdir build
    cd build
    export CC=clang
    export CXX=clang++-libc++
    cmake -DCMAKE_INSTALL_PREFIX=/usr -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
    make -j 4
    make install

# Usage

    $ lp file.lp
