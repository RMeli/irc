#! /bin/bash

echo -n "Creating build directory..."
mkdir -p build && cd build
echo ""done

echo "Current directory: $PWD"

echo "Linear algebra library: ${LIBLA}"

if [[ "$LIBLA" == "arma" ]]; then
        cmake -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DWITH_ARMA:BOOLEAN=TRUE -DCOVERAGE:BOOLEAN=TRUE ..
    elif [[ "$LIBLA" == "eigen" ]]; then
        cmake -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DWITH_EIGEN:BOOLEAN=TRUE -DCOVERAGE:BOOLEAN=TRUE ..
fi

make -j

env CTEST_OUTPUT_ON_FAILURE=1 make -j test