# Internal Redundant Coordinates

[![Build Status](https://travis-ci.org/RMeli/irc.svg?branch=master)](https://travis-ci.org/RMeli/irc)
[![codecov](https://codecov.io/gh/RMeli/irc/branch/master/graph/badge.svg)](https://codecov.io/gh/RMeli/irc)
[![GitHub top language](https://img.shields.io/github/languages/top/RMeli/irc.svg)](https://isocpp.org/)
[![GitHub license](https://img.shields.io/github/license/RMeli/irc.svg)](https://github.com/RMeli/irc/blob/master/LICENSE)

**DISCLAMIER:** *IRC is currently under developement. At the current sage is therefore unfinished, unpolished, not extensively tested and not used in any production code. It is open-sourced so that it is easy to contribute to its development. When the library is ready and safe to use, this disclaimer will be removed and a release will be created.*

IRC is a modern C++ library allowing the determination of internal redundant coordinates (and transformation to and from Cartesian coordinates) for geometry optimization of molecules. The aim of this library is to provide an easy-to-use, flexible and portable implementation of internal redundant coordinates for modern electronic structure codes.

- [Compilation and Installation](https://github.com/RMeli/irc#compilation-and-installation)
- [Usage](https://github.com/RMeli/irc#usage)
- [Test Suite](https://github.com/RMeli/irc#test-suite)
- [Contributions](https://github.com/RMeli/irc#contributions)
- [Sources](https://github.com/RMeli/irc#sources)

## Compilation and Installation

### Dependencies

IRC uses the Boost Graph Library ([BGL](http://www.boost.org/doc/libs/1_66_0/libs/graph/doc/index.html)) to determine the connectivity of the molecule. [Boost](http://www.boost.org/) is therefore a requirement for the library to work.

IRC needs a linear algebra library. Support for [Armadillo](http://arma.sourceforge.net/) and [Eigen](http://eigen.tuxfamily.org) is provided, but other linear algebra libraries can be easily added (see `include/linalg.h`).

### Build IRC
Debug with Armadillo:
```CMake
  mkdir build && cd build
  cmake -CMAKE_BUILD_TYPE=Debug -DWITH_ARMA:BOOLEAN=TRUE ..
  make -j
```

Debug with Eigen:
```CMake
  mkdir build && cd build
  cmake -CMAKE_BUILD_TYPE=Debug -DWITH_EIGEN:BOOLEAN=TRUE ..
  make -j
```

Release with Armadillo:
```CMake
  mkdir build && cd build
  cmake -DWITH_ARMA:BOOLEAN=TRUE ..
  make -j
```

Release with Eigen:
```CMake
  mkdir build && cd build
  cmake -DWITH_EIGEN:BOOLEAN=TRUE ..
  make -j
```

### Test

```Shell
  make -j test
```

### Include in a CMake project

## Usage

IRC is mainly based on free functions, therefore every function can be used out-of-the-box. However, a wrapper class `IRC` with the main functionalities needed in geomerty optimization is provided.

The first step to use IRC is to build a molecule. In IRC a molecule is represented as an `std::vector` of atoms (`irc::atom::Atom<Vector3>`):
```C++
template <typename Vector3>
using irc::molecule::Molecule<Vector3> = std::vector<irc::atoms::Atom<Vector3>>;
```
where an atom is simply defined by an atomic number and its position (stored in a three-dimensional vector of type `Vector3`). In order to use IRC, you will need to implement a function that converts your molecule to an `irc::molecule::Molecule`:
```C++
template<typename Vector3>
irc::molecule::Molecule<Vector3> convert_my_molecule_to_irc_molecule(const MyMolecule&);
```

Once you have an IRC molecule `irc::molecule::Molecule<Vector3>`, you can build an `irc::IRC<Vector3,Vector,Matrix>` object:
```C++
irc::IRC<Vector3, Vector, Matrix>(const irc::molecule::Molecule&);
```
The `irc::IRC` class provides an initial guess of the Hessian in redundant internal coordinates (`Matrix projected_initial_hessian_inv()`), a projector for the update Hessian (`Matrix projected_hessian_inv(const Matrix &)`), the tranformation of the gradient from Cartesian to redundant internal coordinates (`Vector grad_cartesian_to_projected_irc(const Vector &)`) and the transformation from the updated redundant internal coordinates to Cartesian coordinates (`Vector irc_to_cartesian(const Vector &, const Vector &, const Vector &x_c_old, size_t, double)`). In addition, the transformation from Cartesian to redundant internal coordinates is also provided (`Vector cartesian_to_irc(const Vector &)`).

### User-defined internal coordinates

## Tests and code coverage

### Catch2
Tests are written using the multi-paradigm test framework [Catch2](https://github.com/catchorg/Catch2). Catch2 is included as a single header file in `include/catch`.

### CTest
Tests are run using the CTest testing tool distributed as a part of CMake.

Run tests:
```Shell
  make -j test
```
  
### Travis-CI

Continuous integration (CI) is implemented using [Travis-CI](https://travis-ci.org/). The test suite is run for every commit on all the branches and at least once a day for the `master` branch.

The supported operating systems and compilers are the following:

OS | Compiler | Libraries
---|----------|----------
macOS `10.13` | Apple Clang `9.0` | Armadillo
macOS `10.13` | GCC `7.3` | Eigen
macOS `10.12` | Apple Clang `8.1` | Armadillo
macOS `10.12` | GCC `6.4` | Eigen
macOS `10.11` | Apple Clang `7.3` | Armadillo
macOS `10.11` | GCC `5.5` | Eigen


The current build status for every branch is the following:

Branch | Status
-------|-------
`master` | [![Build Status](https://travis-ci.org/RMeli/irc.svg?branch=master)](https://travis-ci.org/RMeli/irc)
`bugfix` | [![Build Status](https://travis-ci.org/RMeli/irc.svg?branch=bugfix)](https://travis-ci.org/RMeli/irc)
`travis-ci` | [![Build Status](https://travis-ci.org/RMeli/irc.svg?branch=travis-ci)](https://travis-ci.org/RMeli/irc)

### Code coverage
Code coverage is obtained using `Apple LLVM` by using the option `--coverage` (for both compiling and linking). The coverage reports are then processed and uploaded to [CodeCov](https://codecov.io) by Travis-CI. The directories `include/boost/` and `/include/catch/` are excluded from the code coverage report.

To compile with code coverage, use the following CMake option:
```CMake
-DCOVERAGE:BOOLEAN=TRUE
```

Manually upload a coverage report:
```Shell
bash <(curl -s https://codecov.io/bash)
```

The current coverage status for every branch is the following:

Branch | Status
-------|-------
`master` | [![codecov](https://codecov.io/gh/RMeli/irc/branch/master/graph/badge.svg)](https://codecov.io/gh/RMeli/irc)
`bugfix` | [![codecov](https://codecov.io/gh/RMeli/irc/branch/bugfix/graph/badge.svg)](https://codecov.io/gh/RMeli/irc)
`travis-ci` | [![codecov](https://codecov.io/gh/RMeli/irc/branch/travis-ci/graph/badge.svg)](https://codecov.io/gh/RMeli/irc)

Sunburst coverage graph for the branch `master`:

![codecov-graph](https://codecov.io/gh/RMeli/irc/branch/master/graphs/sunburst.svg)

### Code format
The code is formatted using `clang-format`.  The style configuration is based on `LLVM` style and saved in the file `tools/.clang-format`.

Format the code:
```Shell
bash tool/clang-format.sh
```

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
