#include "catch.hpp"

#include "libirc/molecule.h"

#include "libirc/periodic_table.h"

#include <cassert>
#include <iostream>

#ifdef HAVE_ARMA
#include <armadillo>
using arma::vec3;
#elif HAVE_EIGEN3
#include <Eigen3/Eigen/Dense>
using vec3 = Eigen::Vector3d;
#else
#error
#endif

using namespace irc;

TEST_CASE("Molecule") {

  using namespace std;
  using namespace molecule;
  using namespace periodic_table;

  const auto molecule = Molecule<vec3>{
      {{1, {0.0, 1.1, 2.2}}, {2, {0.0, 1.1, 2.2}}, {3, {0.0, 1.1, 2.2}}}};

  CHECK(mass(molecule) == Approx(masses[1] + masses[2] + masses[3]));

  SECTION("Position multiplier") {
    auto multiplier = 2.;

    auto scaled_molecule = molecule;
    multiply_positions(scaled_molecule, multiplier);

    for (std::size_t i{0}; i < scaled_molecule.size(); i++) {
      for (std::size_t j{0}; j < 3; j++) {
        CHECK(scaled_molecule[i].position(j) ==
              Approx(multiplier * molecule[i].position(j)));
      }
    }

    scaled_molecule = molecule * multiplier;
    for (std::size_t i{0}; i < scaled_molecule.size(); i++) {
      for (std::size_t j{0}; j < 3; j++) {
        CHECK(scaled_molecule[i].position(j) ==
              Approx(multiplier * molecule[i].position(j)));
      }
    }

    scaled_molecule = multiplier * molecule;
    for (std::size_t i{0}; i < scaled_molecule.size(); i++) {
      for (std::size_t j{0}; j < 3; j++) {
        CHECK(scaled_molecule[i].position(j) ==
              Approx(multiplier * molecule[i].position(j)));
      }
    }
  }
}
