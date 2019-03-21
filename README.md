# Internal Redundant Coordinates

[![Build Status](https://travis-ci.org/RMeli/irc.svg?branch=master)](https://travis-ci.org/RMeli/irc)
[![codecov](https://codecov.io/gh/RMeli/irc/branch/master/graph/badge.svg)](https://codecov.io/gh/RMeli/irc)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/6cd91c20ae064924bc5c4e0181710bf3)](https://www.codacy.com/app/RMeli/irc?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=RMeli/irc&amp;utm_campaign=Badge_Grade)
[![GitHub license](https://img.shields.io/github/license/RMeli/irc.svg)](https://github.com/RMeli/irc/blob/master/LICENSE)
[![GitHub top language](https://img.shields.io/github/languages/top/RMeli/irc.svg)](https://isocpp.org/)

**DISCLAMIER:** *IRC is currently under developement. It is integrated in the new quantum chemistry code [entos](https://www.doi.org/10.26434/chemrxiv.7762646.v2), soon to be released.*

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

### Include IRC in your project

IRC is an header-only library and its inclusion is quite straightforward.

### With CMake

### Without CMake

If you are using Armadillo as linear algebra library you just need to include the header files in `include/` and define the variable `HAVE_ARMA`.


If you are using Eigen3 as linear algebra library, you need to include the header files in `include/` as well as the extension to Eigen3's matrix constructors to support `std::initializer_list`s located in `external/eigen/plugins/`. The content of the file `Matrix_initializer_list.h` must be included Eigen3's matrix class as described [here](http://eigen.tuxfamily.org/dox-3.2/TopicCustomizingEigen.html), therefore you need to set the variable `EIGEN_MATRIXBASE_PLUGIN` to `"path_to/Matrix_initializer_list.h"`. Finally, you have to define the variable `HAVE_EIGEN3`.

If your are using a custom linear algebra library supporting the initialization of vectors and matrices from `std::initializer_list`s, you have to implement the linear algebra functions in `include/libirc/linalg.h` and you are good to go.

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

The current build status for every branch is the following:

Branch | Status
-------|-------
`master` | [![Build Status](https://travis-ci.org/RMeli/irc.svg?branch=master)](https://travis-ci.org/RMeli/irc)

### Code coverage
Code coverage is obtained using `Apple LLVM` with the option `--coverage` (for both compiling and linking). The coverage reports are then processed and uploaded to [CodeCov](https://codecov.io) by Travis-CI. The directories `include/boost/` and `/include/catch/` are excluded from the code coverage report.

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
Peter Bygrave | University of Bristol

</center>

## Sources

### Papers

- P. Puly and G. Fogarasi, *Geometry optimization in redundant internal coordinates*, J. Chem. Phys. **96** 2856 (1992).

- C. Peng, P. Y. Ayala and H. B. Schlegel, *Using Redundant Internal Coordinates to Optimize Equilibrium Geometries and Transition States*, J. Comp. Chem. **17**, 49-56 (1996).

- V. Bakken and T. Helgaker, *The efficient optimization of molecular geometries using redundant internal coordinates*, J. Chem. Phys. **117**, 9160 (2002).

### Books

- E. Bright Wilson Jr., J. C. Decius and P. C. Cross, *Molecular Vibrations: The Theory of Infrared and Raman Vibrational Spectra*, Dover Publications Inc. (2003).
