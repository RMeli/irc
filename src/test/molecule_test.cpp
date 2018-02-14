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
    auto scaled_molecule = molecule;
    multiply_positions(scaled_molecule, 2.);

    const auto pos = vec3{0.0, 2.2, 4.4};

    for (const auto& atom : scaled_molecule) {
      for (size_t i{0}; i < 3; i++) {
        CAPTURE(i);
        CHECK(atom.position(i) == Approx(pos(i)));
      }
    }
  }
}
