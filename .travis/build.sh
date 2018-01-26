#! /bin/bash

mkdir -p build

cd build

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    cmake -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DWITH_ARMA:BOOLEAN=TRUE ..
elif [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    #cmake -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DWITH_EIGEN:BOOLEAN=TRUE ..
fi

make -j

make -j test