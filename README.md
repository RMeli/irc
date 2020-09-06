# Internal Redundant Coordinates

[![Build Status](https://travis-ci.org/RMeli/irc.svg?branch=master)](https://travis-ci.org/RMeli/irc)
[![codecov](https://codecov.io/gh/RMeli/irc/branch/master/graph/badge.svg)](https://codecov.io/gh/RMeli/irc)
[![GitHub license](https://img.shields.io/github/license/RMeli/irc.svg)](https://github.com/RMeli/irc/blob/master/LICENSE)
[![GitHub top language](https://img.shields.io/github/languages/top/RMeli/irc.svg)](https://isocpp.org/)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3928502.svg)](https://doi.org/10.5281/zenodo.3928502)

**DISCLAMIER:** *IRC is still in beta. It is integrated in the quantum chemistry code [entos](https://www.entos.info/).*

IRC is a modern C++ library for the determination of internal redundant coordinates (and transformation to and from Cartesian coordinates) for geometry optimization of molecules. The aim of this library is to provide an easy-to-use, flexible and portable implementation of internal redundant coordinates for modern electronic structure codes.

- [Compilation and Installation](https://github.com/RMeli/irc#compilation-and-installation)
- [Usage](https://github.com/RMeli/irc#usage)
- [Test Suite](https://github.com/RMeli/irc#test-suite)
- [Contributions](https://github.com/RMeli/irc#contributions)
- [Sources](https://github.com/RMeli/irc#sources)

## Compilation and Installation

### Dependencies

IRC has the following dependencies:

- Boost Graph Library ([BGL](http://www.boost.org/doc/libs/1_66_0/libs/graph/doc/index.html))
- A linear algebra library ([Armadillo](http://arma.sourceforge.net/) or [Eigen](http://eigen.tuxfamily.org))

[BGL](http://www.boost.org/doc/libs/1_66_0/libs/graph/doc/index.html) functionality is used determine the connectivity of the molecule and therefore the [Boost](http://www.boost.org/) is a requirement for the library to work.

IRC also needs a linear algebra library. [Armadillo](http://arma.sourceforge.net/) and [Eigen](http://eigen.tuxfamily.org) are supported out-of-the-box, but other linear algebra libraries can be easily added.

For standalone installation and testing, IRC also requires [CMake](https://cmake.org/).

### Build IRC

#### Armadillo

```cmake
  mkdir build && cd build
  cmake -DWITH_ARMA=TRUE ..
  make -j
```

#### Eigen

```cmake
  mkdir build && cd build
  cmake -DWITH_EIGEN=TRUE ..
  make -j
```

### Test

```bash
  ctest
```

### Include IRC in your project

IRC is an header-only library and its inclusion is quite straightforward.

### With CMake

Libirc installs CMake configuration files that enable use of the library with
`find_package`. The irc CMake target is named `irc`. Armadillo and Eigen
dependencies are resolved at configuration time of this library and the
relevant compile definitions and includes are set in the configuration files.

### Without CMake

#### Armadillo

If you are using Armadillo as linear algebra library you just need to include the header files in `include/` and define the variable `HAVE_ARMA`.

#### Eigen

If you are using Eigen3 as linear algebra library, you need to include the header files in `include/` as well as the extension to Eigen3's matrix constructors to support `std::initializer_list`s located in `external/eigen/plugins/`. The content of the file `Matrix_initializer_list.h` must be included Eigen3's matrix class as described [here](http://eigen.tuxfamily.org/dox-3.2/TopicCustomizingEigen.html), therefore you need to set the variable `EIGEN_MATRIXBASE_PLUGIN` to `"path_to/Matrix_initializer_list.h"`. Finally, you have to define the variable `HAVE_EIGEN3`.

#### Custom Linear Algebra Library

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

The `irc::IRC` class provides an initial guess of the Hessian in redundant internal coordinates (`Matrix projected_initial_hessian_inv()`), a projector for the update Hessian (`Matrix projected_hessian_inv(const Matrix &)`), the transformation of the gradient from Cartesian to redundant internal coordinates (`Vector grad_cartesian_to_projected_irc(const Vector &)`) and the transformation from the updated redundant internal coordinates to Cartesian coordinates (`Vector irc_to_cartesian(const Vector &, const Vector &, const Vector &x_c_old, size_t, double)`). In addition, the transformation from Cartesian to redundant internal coordinates is also provided (`Vector cartesian_to_irc(const Vector &)`).

### User-defined internal coordinates

## Tests and Code Coverage

### Catch2

Tests are written using the multi-paradigm test framework [Catch2](https://github.com/catchorg/Catch2). Catch2 is included as a single header file in `include/catch`.

### CTest

Tests are run using the CTest testing tool distributed as a part of CMake.

To run the test use:

```bash
  ctest
```
  
### Travis-CI

Continuous integration (CI) is implemented using [Travis-CI](https://travis-ci.org/). The test suite is run for every commit on all the branches and at least once a day for the `master` branch.

### Code Coverage

Code coverage is obtained using `Apple LLVM`. The coverage reports are then processed and uploaded to [CodeCov](https://codecov.io) by Travis-CI.

### Code Formatting

The code is formatted using `clang-format`.  The style configuration is based on `LLVM` style.

To format your code before opening a PR use:

```bash
bash tool/clang-format.sh
```

## Contributions

Any contribution to this open-source project is very welcome. If you are considering contributing don't hesitate to contact the main contributors.

You may find beneficial to have a look at the [Open Source Guides](https://opensource.guide/).

## Sources

### Papers

- P. Puly and G. Fogarasi, *Geometry optimization in redundant internal coordinates*, J. Chem. Phys. **96** 2856 (1992).

- C. Peng, P. Y. Ayala and H. B. Schlegel, *Using Redundant Internal Coordinates to Optimize Equilibrium Geometries and Transition States*, J. Comp. Chem. **17**, 49-56 (1996).

- V. Bakken and T. Helgaker, *The efficient optimization of molecular geometries using redundant internal coordinates*, J. Chem. Phys. **117**, 9160 (2002).

### Books

- E. Bright Wilson Jr., J. C. Decius and P. C. Cross, *Molecular Vibrations: The Theory of Infrared and Raman Vibrational Spectra*, Dover Publications Inc. (2003).
