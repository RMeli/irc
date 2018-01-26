#! /bin/bash

mkdir -p build

cd build

cmake -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX ..

make -j

make -j test