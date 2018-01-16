# Internal Redundant Coordinates

[![Build Status](https://travis-ci.org/RMeli/irc.svg?branch=master)](https://travis-ci.org/RMeli/irc)
[![codecov](https://codecov.io/gh/RMeli/irc/branch/master/graph/badge.svg)](https://codecov.io/gh/RMeli/irc)
[![License: MIT](https://img.shields.io/packagist/l/doctrine/orm.svg)](https://opensource.org/licenses/MIT)

## Usage

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

<center>
  
Branch | Status
-------|-------
`master` | [![Build Status](https://travis-ci.org/RMeli/irc.svg?branch=master)](https://travis-ci.org/RMeli/irc)
`travis-ci` | [![Build Status](https://travis-ci.org/RMeli/irc.svg?branch=travis-ci)](https://travis-ci.org/RMeli/irc)

</center>

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
