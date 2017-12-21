#! /bin/bash

# build directory
mkdir -p build
cd build

# CMake
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Make
make -j

# Tests
make -j test
