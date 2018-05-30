#!/bin/bash

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    echo "Installing ${LIBLA}..."
    if [[ "$LIBLA" == "arma" ]]; then
        brew install armadillo
    elif [[ "$LIBLA" == "eigen" ]]; then
        brew install eigen
    fi
fi
