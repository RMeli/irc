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
    vec x_c_from_irc{irc.irc_to_cartesian(q_irc, dq, x_c).x_c};

    std::size_t n{linalg::size(x_c)};
    for (std::size_t i{0}; i < n; i++) {
      Approx target(x_c(i));
      target.margin(1e-6);

      REQUIRE(x_c_from_irc(i) == target);
    }
  }
}

TEST_CASE("Issues") {
  using namespace molecule;
  using namespace connectivity;
  using namespace transformation;
  using namespace tools::conversion;
  using namespace io;
  using namespace std;

  SECTION("Issue #41") {
    // Define molecule in Issue #41
    const auto molecule = load_xyz<vec3>(config::molecules_dir + "issue41.xyz");

    // Build internal reaction coordinates
    IRC<vec3, vec, mat> ircs(molecule);

    vec q_irc_next = {
        2.69873,    2.63216,     2.75051,     2.7331,      2.6465,
        2.66307,    2.83097,     2.74455,     2.58185,     2.84907,
        2.82958,    2.75401,     2.58345,     2.33196,     2.4169,
        2.08872,    2.10057,     2.08658,     2.09745,     2.08976,
        2.1038,     2.09314,     2.09084,     2.07154,     2.09745,
        2.10667,    2.09451,     2.13115,     2.06417,     2.04088,
        2.11113,    2.13647,     2.08623,     2.04094,     2.10441,
        1.95211,    2.0311,      1.96448,     2.02538,     2.09322,
        2.13935,    2.08234,     2.10632,     2.14114,     2.11665,
        2.08578,    2.09212,     2.17156,     2.07068,     2.10512,
        2.0892,     2.12405,     2.06167,     2.10058,     2.11018,
        2.10038,    2.07596,     2.08923,     2.12978,     0.0192478,
        0.0107708,  0.0258413,   0.0163308,   -0.0410394,  -0.0314206,
        3.09261,    -3.1075,     1.31824,     -1.83119,    -1.95765,
        0.0562704,  0.00600911,  1.9557,      -0.638266,   -0.683332,
        -1.21737,   1.91638,     0.598524,    0.645881,    2.0076,
        -2.03129,   -3.12078,    3.10563,     -3.14202,    -2.48917,
        1.11684,    -1.14069,    2.54853,     3.10222,     -0.04581,
        -3.08914,   2.46185,     -1.13041,    1.26207,     3.12046,
        -2.46494,   -0.0249528,  -3.12882,    -3.11901,    0.0270607,
        3.11603,    -3.13808,    -0.00442377, -0.0320436,  -3.1364,
        -3.11046,   0.0134492,   -3.13311,    3.11241,     -0.0337284,
        0.0288444,  -0.00699303, 0.0393071,   -0.00540334, 0.00393414,
        0.00201568, -0.0175075,  0.0673447,   -0.00320022, 0.00578351,
        0.0390592,
    };

    // Get cartesian coordinates
    vec x_c{to_cartesian<vec3, vec>(molecule)};

    // Compute internal redundant coordinates
    vec q_irc{ircs.cartesian_to_irc(x_c)};

    cout << q_irc << endl;

    // Define no displacement in IRC
    vec dq{q_irc_next - q_irc};

    // Compute cartesian coordinate from IRC
    auto result = ircs.irc_to_cartesian(q_irc_next, dq, x_c);

    vec x_c_from_irc{result.x_c};

    CAPTURE(ircs.get_bonds().size());
    CAPTURE(ircs.get_angles().size());
    CAPTURE(ircs.get_dihedrals().size());
    CAPTURE(ircs.get_linear_angles().size());
    CAPTURE(ircs.get_out_of_plane_bends().size());

    CHECK(result.converged);
  }
}
