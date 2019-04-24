#include "catch.hpp"

#include "config.h"

#include "libirc/irc.h"

#include "libirc/conversion.h"
#include "libirc/io.h"
#include "libirc/molecule.h"

#ifdef HAVE_ARMA
#include <armadillo>
using vec3 = arma::vec3;
using vec = arma::vec;
using mat = arma::mat;

template<typename T>
using Mat = arma::Mat<T>;
#elif HAVE_EIGEN3
#include <Eigen3/Eigen/Dense>
using vec3 = Eigen::Vector3d;
using vec = Eigen::VectorXd;
using mat = Eigen::MatrixXd;

template<typename T>
using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
#else
#error
#endif

using namespace irc;

TEST_CASE("Internal Redundant Coordinates") {
  using namespace std;

  using namespace connectivity;
  using namespace molecule;
  using namespace tools::conversion;

  SECTION("User-defined Coordinates") {

    // Define methanol molecule
    Molecule<vec3> molecule{{"O", {0.0000000000, 0.0000000000, 0.0000000000}},
                            {"H", {-0.9658081475, 0.0000000000, 0.0000000000}},
                            {"C", {0.4371881141, 1.3664973705, 0.0000000000}},
                            {"H", {0.0965595997, 1.9062256626, 0.8987823813}},
                            {"H", {1.5337031672, 1.3396209523, 0.0000000000}},
                            {"H", {0.0965596000, 1.9062256619, -0.8987823831}}};

    // Transform molecular coordinates from angstrom to bohr
    multiply_positions(molecule, angstrom_to_bohr);

    // Build internal reaction coordinates
    // Add H-H bonds: (3,4), (4,5), (5,3)
    // Add H-O-H bonds: (1,0,3), (1,0,4), (1,0,5)
    // Add H-H-H-H dihedral: (1,3,4,5)
    // Add O-H-H-H dihedral: (0,3,4,5)
    IRC<vec3, vec, mat> irc(molecule,
                            {{3, 4}, {4, 5}, {5, 3}},
                            {{1, 0, 3}, {1, 0, 4}, {1, 0, 5}},
                            {{1, 3, 4, 5}, {0, 3, 4, 5}});

    // Compute internal coordinates
    vec q_irc{
        irc.cartesian_to_irc(molecule::to_cartesian<vec3, vec>(molecule))};

    // Check size (5+3 bonds, 7+3 angles, 3+2 dihedrals)
    REQUIRE(linalg::size(q_irc) == 23);

    // Check manually added H-H bond
    SECTION("Manually added H-H bond (1)") {
      Approx target(q_irc(5));
      target.margin(1e-6);

      // Compute O-H distance
      double d{distance(molecule[3].position, molecule[4].position)};

      CHECK(d == target);
    }

    // Check manually added H-H bond
    SECTION("Manually added H-H bond (2)") {
      Approx target(q_irc(6));
      target.margin(1e-6);

      // Compute O-H distance
      double d{distance(molecule[4].position, molecule[5].position)};

      CHECK(d == target);
    }

    // Check manually added H-H bond
    SECTION("Manually added H-H bond (3)") {
      Approx target(q_irc(7));
      target.margin(1e-6);

      // Compute O-H distance
      double d{distance(molecule[5].position, molecule[3].position)};

      CHECK(d == target);
    }

    // Check manually added H-O-H angle (1)
    SECTION("Manually added H-O-H angle (1)") {
      Approx target(q_irc(15));
      target.margin(1e-6);

      // Compute H-O-H angle
      double a{angle(
          molecule[1].position, molecule[0].position, molecule[3].position)};

      CHECK(a == target);
    }

    // Check manually added H-O-H angle (2)
    SECTION("Manually added H-O-H angle (2)") {
      Approx target(q_irc(16));
      target.margin(1e-6);

      // Compute H-O-H angle
      double a{angle(
          molecule[1].position, molecule[0].position, molecule[4].position)};

      CHECK(a == target);
    }

    // Check manually added H-O-H angle (3)
    SECTION("Manually added H-O-H angle (3)") {
      Approx target(q_irc(17));
      target.margin(1e-6);

      // Compute H-O-H angle
      double a{angle(
          molecule[1].position, molecule[0].position, molecule[5].position)};

      CHECK(a == target);
    }

    // Check manually added H-H-H-H dihedral angle
    SECTION("Manually added H-H-H-H dihedral angle") {
      Approx target(q_irc(21));
      target.margin(1e-6);

      // Compute H-O-C-H dihedral angle
      double d{dihedral(molecule[1].position,
                        molecule[3].position,
                        molecule[4].position,
                        molecule[5].position)};

      CHECK(d == target);
    }

    // Check manually added O-H-H-H dihedral angle
    SECTION("Manually added O-H-H-H dihedral angle") {
      Approx target(q_irc(22));
      target.margin(1e-6);

      // Compute H-O-C-H dihedral angle
      double d{dihedral(molecule[0].position,
                        molecule[3].position,
                        molecule[4].position,
                        molecule[5].position)};

      CHECK(d == target);
    }

    // Check manually added H-H-O-C dihedral angle
    SECTION("Manually added H-H-O-C out of plane angle") {
      Approx target(q_irc(10));
      target.margin(1e-6);

      // Compute H-H-O-C dihedral angle
      double d{out_of_plane_angle(molecule[0].position,
                                  molecule[1].position,
                                  molecule[2].position,
                                  molecule[3].position)};

      CHECK(d == target);
      CHECK(d == Approx(0).margin(1e-6));
    }
  }
  /*
  SECTION("Constraints") {
    // Define formaldehyde molecule (CH2O)
    Molecule<vec3> molecule{{"C", {0.000000, 0.000000, -0.537500}},
                            {"O", {0.000000, 0.000000, 0.662500}},
                            {"H", {0.000000, 0.866025, -1.037500}},
                            {"H", {0.000000, -0.866025, -1.037500}}};

    // Transform molecular coordinates from angstrom to bohr
    multiply_positions(molecule, angstrom_to_bohr);

    // Build internal reaction coordinates
    // Add C-O bond constraint: (0,1)
    // Add H-H bond constraint: (2,3)
    // Add H-C-H angle constraint: (2,0,3)
    // Add H-H-O-C dihedral constraint: (3,2,1,0)
    IRC<vec3, vec, mat> irc(
        molecule, {{0,1,Constraint::constrained},
  {2,3,Constraint::constrained}}, {{2,0,3,Constraint::constrained}},
  {{3,2,1,0,Constraint::constrained}});

    // Compute internal coordinates
    vec q_irc{
        irc.cartesian_to_irc(molecule::to_cartesian<vec3, vec>(molecule))};

    // Check size (3+1 bonds, 3+0 angles, 0+1 dihedrals)
    REQUIRE(linalg::size(q_irc) == 8);

    // Get bonds
    auto B = irc.get_bonds();

    // Check bonds
    REQUIRE(B.size() == 4 );
    CHECK(B[0].constraint == Constraint::constrained);
    CHECK(B[3].constraint == Constraint::constrained);

    // Get angles
    auto A = irc.get_angles();

    // Check bonds
    REQUIRE(A.size() == 3 );
    CHECK(A[2].constraint == Constraint::constrained);

    // Get angles
    //auto D = irc.get_dihedrals();

    // Check bonds
    //REQUIRE(D.size() == 1 );
    CHECK(D[0].constraint == Constraint::constrained);
  }
  */

  SECTION("Hessian projection") {

    // Define methanol molecule
    Molecule<vec3> molecule{{"O", {0.0000000000, 0.0000000000, 0.0000000000}},
                            {"H", {-0.9658081475, 0.0000000000, 0.0000000000}},
                            {"C", {0.4371881141, 1.3664973705, 0.0000000000}},
                            {"H", {0.0965595997, 1.9062256626, 0.8987823813}},
                            {"H", {1.5337031672, 1.3396209523, 0.0000000000}},
                            {"H", {0.0965596000, 1.9062256619, -0.8987823831}}};

    // Transform molecular coordinates from angstrom to bohr
    multiply_positions(molecule, angstrom_to_bohr);

    // Build internal reaction coordinates
    // (Manually added dihedral to increase coverage)
    IRC<vec3, vec, mat> irc(molecule, {}, {}, {{0, 1, 2, 3}});

    // Compute initial hessian
    mat iH0{irc.projected_initial_hessian_inv()};

    // Project Hessian again
    mat iH{irc.projected_hessian_inv(iH0)};

    // Check sizes
    REQUIRE(linalg::size(iH0) == linalg::size(iH));

    // Check that second projection has no effect
    std::size_t n{linalg::size(iH0)};
    for (std::size_t i{0}; i < n; i++) {
      Approx target(iH0(i));
      target.margin(1e-6);

      CHECK(iH(i) == target);
    }
  }

  SECTION("IRC to Cartesian") {

    // Define formaldehyde molecule (CH2O)
    const auto molecule =
        io::load_xyz<vec3>(config::molecules_dir + "ethanol.xyz");

    // Build internal reaction coordinates
    IRC<vec3, vec, mat> irc(molecule);

    // Get cartesian coordinates
    vec x_c{to_cartesian<vec3, vec>(molecule)};

    // Compute internal redundant coordinates
    vec q_irc{irc.cartesian_to_irc(x_c)};

    // Define no displacement in IRC
    vec dq{linalg::zeros<vec>(linalg::size(q_irc))};

    // Compute cartesian coordinate from IRC
    vec x_c_from_irc{irc.irc_to_cartesian(q_irc, dq, x_c)};

    std::size_t n{linalg::size(x_c)};
    for (std::size_t i{0}; i < n; i++) {
      Approx target(x_c(i));
      target.margin(1e-6);

      REQUIRE(x_c_from_irc(i) == target);
    }
  }
}
