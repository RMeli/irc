#! /bin/bash

mkdir -p build

cd build

echo "Building with ${LIBLA}..."

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    if [[ "$LIBLA" == "arma" ]]; then
        cmake -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DWITH_ARMA:BOOLEAN=TRUE -DCOVERAGE:BOOLEAN=TRUE ..
    elif [[ "$LIBLA" == "eigen" ]]; then
        cmake -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DWITH_EIGEN:BOOLEAN=TRUE -DCOVERAGE:BOOLEAN=TRUE ..
    fi
fi

make -j

make -j test

cd ..