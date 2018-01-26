# Internal Redundant Coordinates

[![Build Status](https://travis-ci.org/RMeli/irc.svg?branch=master)](https://travis-ci.org/RMeli/irc)
[![codecov](https://codecov.io/gh/RMeli/irc/branch/master/graph/badge.svg)](https://codecov.io/gh/RMeli/irc)
[![GitHub top language](https://img.shields.io/github/languages/top/RMeli/irc.svg)](https://isocpp.org/)
[![GitHub license](https://img.shields.io/github/license/RMeli/irc.svg)](https://github.com/RMeli/irc/blob/master/LICENSE)


IRC is a modern C++ library allowing the determination of internal redundant coordinates (and transformation to and from Cartesian coordinates) for geometry optimization of molecules. The aim of this library is to provide an easy-to-use and portable implementation of internal redundant coordinates for modern electronic structure codes.

- [Compilation and Installation](https://github.com/RMeli/irc#compilation-and-installation)
- [Usage](https://github.com/RMeli/irc#usage)
- [Test Suite](https://github.com/RMeli/irc#test-suite)
- [Contributions](https://github.com/RMeli/irc#contributions)
- [Sources](https://github.com/RMeli/irc#sources)

## Compilation and Installation

### Dependencies

IRC uses Boost Graph Library ([BGL](http://www.boost.org/doc/libs/1_66_0/libs/graph/doc/index.html)) to determine the connectivity of the molecule. To simplify inclusion of IRC in other projects, BGL is provided in `inlcude/boost`. If you want to use your local Boost library, just remove `include/boost`.

For tests, IRC needs a linear algebra library. For the time being, the use of [Armadillo](http://arma.sourceforge.net/) is required but support for [Eigen](http://eigen.tuxfamily.org) is also provided.

### Build
Debug:
```
  mkdir build && cd build
  cmake -CMAKE_BUILD_TYPE=Debug ..
  make -j
```

Release:
```
  mkdir build && cd build
  cmake ..
  make -j
```

### Test

```
  make -j test
```

### Build and test

```
  bash build.sh
```

### Install

### Include in a CMake project

## Usage

## Test suite

### Catch2
Tests are written using the multi-paradigm test framework [Catch2](https://github.com/catchorg/Catch2). Catch2 is included as a single header file in `include/catch`.

### CTest
Tests are run using the CTest testing tool distributed as a part of CMake.

Run tests:
```
  make -j test
```
  
### Travis-CI

Continuous integration (CI) is implemented using [Travis-CI](https://travis-ci.org/). The test suite is run for every commit on all the branches and at least once a day for the `master` branch.

The supported operating systems and compilers are the following:

OS | Compiler
---|---------
macOS 10.13 | Apple Clang 9.0
macOS 10.13 | GCC 7.3
macOS 10.12 | Apple Clang 8.1
macOS 10.12 | GCC 6.4
macOS 10.11 | Apple Clang 7.3
macOS 10.11 | GCC 5.5


The current build status for every branch is the following:

Branch | Status
-------|-------
`master` | [![Build Status](https://travis-ci.org/RMeli/irc.svg?branch=master)](https://travis-ci.org/RMeli/irc)
`travis-ci` | [![Build Status](https://travis-ci.org/RMeli/irc.svg?branch=travis-ci)](https://travis-ci.org/RMeli/irc)

## Contributions

Any contribution to this open-source project is very welcome. If you are considering contributing don't hesitate to contact the main constributors. 

You may find beneficial to have a look at the [Open Source Guides](https://opensource.guide/).

### List of contributors
<center>
  
Contributor | Affiliation
------------|-------------
Rocco Meli | University of Bristol

</center>

## Sources

### Papers

- P. Puly and G. Fogarasi, *Geometry optimization in redundant internal coordinates*, J. Chem. Phys. **96** 2856 (1992).

- C. Peng, P. Y. Ayala and H. B. Schlegel, *Using Redundant Internal Coordinates to Optimize Equilibrium Geometries and Transition States*, J. Comp. Chem. **17**, 49-56 (1996).

- V. Bakken and T. Helgaker, *The efficient optimization of molecular geometries using redundant internal coordinates*, J. Chem. Phys. **117**, 9160 (2002).

### Books

- E. Bright Wilson Jr., J. C. Decius and P. C. Cross, *Molecular Vibrations: The Theory of Infrared and Raman Vibrational Spectra*, Dover Publications Inc. (2003).
