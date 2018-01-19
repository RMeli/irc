#! /bin/bash

# build directory
mkdir -p build
cd build

# CMake
cmake ..

# Make
make -j

# Tests
make -j test
