#include "catch.hpp"

#include "libirc/atom.h"

#include "libirc/periodic_table.h"

#include <iostream>
#include <stdexcept>

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

TEST_CASE("Test atom and periodic table lookup functions", "[atom]") {

  using namespace atom;

  SECTION("invalid atom") {
    CHECK_THROWS_AS(Atom<vec3>(periodic_table::pt_size + 1), std::logic_error);
  }

  SECTION("atom from atomic number") {
    for (size_t i{1}; i < periodic_table::pt_size; i++) {
      REQUIRE(periodic_table::valid_atomic_number(i));

      const auto a = Atom<vec3>{i};

      REQUIRE(symbol(a.atomic_number) == periodic_table::pt_symbols[i]);
      REQUIRE(symbol(i) == periodic_table::pt_symbols[i]);

      REQUIRE(mass(a.atomic_number) == Approx(periodic_table::pt_masses[i]));
      REQUIRE(mass(i) == Approx(periodic_table::pt_masses[i]));

      REQUIRE(covalent_radius(a.atomic_number) == Approx(periodic_table::pt_covalent_radii[i]));
      REQUIRE(covalent_radius(i) == Approx(periodic_table::pt_covalent_radii[i]));

      REQUIRE(vdw_radius(a.atomic_number) == Approx(periodic_table::pt_vdv_radii[i]));
      REQUIRE(vdw_radius(i) == Approx(periodic_table::pt_vdv_radii[i]));

      if (i == 1) {
        REQUIRE(is_H(a.atomic_number));
        REQUIRE(is_H(i));
      } else if (i == 7 or i == 8 or i == 9 or i == 15 or i == 16 or
                 i == 17) {
        REQUIRE(is_NOFPSCl(a.atomic_number));
        REQUIRE(is_NOFPSCl(i));
      }
    }
  }

  SECTION("atom from atomic symbol") {
    for (size_t i{1}; i < periodic_table::pt_size; i++) {

      REQUIRE(periodic_table::valid_atomic_number(i));

      const auto a = Atom<vec3>{periodic_table::pt_symbols[i]};

      REQUIRE(symbol(a.atomic_number) == periodic_table::pt_symbols[i]);
      REQUIRE(symbol(AtomicNumber(periodic_table::pt_symbols[i])) == periodic_table::pt_symbols[i]);

      REQUIRE(mass(a.atomic_number) == Approx(periodic_table::pt_masses[i]));
      REQUIRE(mass(AtomicNumber(periodic_table::pt_symbols[i])) == Approx(periodic_table::pt_masses[i]));

      REQUIRE(covalent_radius(a.atomic_number) == Approx(periodic_table::pt_covalent_radii[i]));
      REQUIRE(covalent_radius(AtomicNumber(periodic_table::pt_symbols[i])) == Approx(periodic_table::pt_covalent_radii[i]));

      REQUIRE(vdw_radius(a.atomic_number) == Approx(periodic_table::pt_vdv_radii[i]));
      REQUIRE(vdw_radius(AtomicNumber(periodic_table::pt_symbols[i])) == Approx(periodic_table::pt_vdv_radii[i]));

      // Atom in H-bond
      if (i == 1) {
        REQUIRE(is_H(a.atomic_number));
      } else if (i == 7 or i == 8 or i == 9 or i == 15 or i == 16 or
                 i == 17) {
        REQUIRE(is_NOFPSCl(a.atomic_number));
      }
    }
  }
}
