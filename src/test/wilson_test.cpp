#include "catch.hpp"

#include "libirc/wilson.h"

#include "libirc/atom.h"
#include "libirc/connectivity.h"
#include "libirc/conversion.h"
#include "libirc/io.h"
#include "libirc/molecule.h"

#include "config.h"

#include <cmath>

#ifdef HAVE_ARMA
#include <armadillo>
using vec3 = arma::vec3;
using vec = arma::vec;
using mat = arma::mat;
#elif HAVE_EIGEN3
#include <Eigen3/Eigen/Dense>
using vec3 = Eigen::Vector3d;
using vec = Eigen::VectorXd;
using mat = Eigen::MatrixXd;
#else
#error
#endif

using namespace irc;

TEST_CASE("Wilson B matrix for single fragments", "[wilson]") {
  using namespace connectivity;
  using namespace molecule;
  using namespace wilson;
  using namespace tools;

  SECTION("H2 stretching") {
    // Define molecule
    Molecule<vec3> mol{{"H", {0., 0., 0.}}, {"H", {1., 0., 0.}}};

    // Molecular connectivity
    mat dd{distances<vec3, mat>(mol)};
    UGraph adj{adjacency_matrix(dd, mol)};
    mat dist{distance_matrix<mat>(adj)};

    // Compute bonds
    std::vector<Bond> B{bonds(dist, mol)};
    REQUIRE(B.size() == 1);

    // Compute Wilson B matrix for H2 analytically
    mat Bwilson =
        wilson_matrix<vec3, vec, mat>(to_cartesian<vec3, vec>(mol), B);

    // Check Wilson B matrix size
    INFO("Wilson B matrix (analytical):\n" << Bwilson);
    REQUIRE(linalg::size(Bwilson) == 6);

    // Compute Wilson B matrix for H2 numerically
    mat BwilsonN = wilson_matrix_numerical<vec3, vec, mat>(
        to_cartesian<vec3, vec>(mol), B);

    // Check Wilson B matrix size
    INFO("Wilson B matrix (numerical):\n" << BwilsonN);
    REQUIRE(linalg::size(BwilsonN) == 6);

    // Check analytical and numerical Wilson matrices are the same
    INFO("Analytical vs Numerical");
    for (std::size_t i{0}; i < 6; i++)
      REQUIRE(Bwilson(i) == Approx(BwilsonN(i)).margin(1e-6));

    INFO("Transformation with bond stretch");
    const double d{0.01};
    vec dx{-d, 0.00, 0.00, d, 0.00, 0.00};

    INFO("Analytical transformation");
    vec analytical_transformation = Bwilson * dx;
    CAPTURE(analytical_transformation);
    REQUIRE(analytical_transformation(0) == Approx(2 * d).margin(1e-5));

    INFO("Numerical transformation");
    vec numerical_transformation = BwilsonN * dx;
    CAPTURE(numerical_transformation);
    REQUIRE(numerical_transformation(0) == Approx(2 * d).margin(1e-5));
  } // H2 stretching

  SECTION("H2O bending") {

    double angle(0.5);
    double angle_rad(angle / 180. * constants::pi);

    const std::vector<mat> R{
        {// Rotation for H1
         {cos(angle_rad), -sin(angle_rad), 0},
         {sin(angle_rad), cos(angle_rad), 0},
         {0, 0, 1}},
        {// Rotation for O
         {1, 0, 0},
         {0, 1, 0},
         {0, 0, 1}},
        {// Rotation for H2
         {cos(-angle_rad), -sin(-angle_rad), 0},
         {sin(-angle_rad), cos(-angle_rad), 0},
         {0, 0, 1}},
    };

    const Molecule<vec3> mol{
        {"H", {1.43, -1.10, 0.00}}, // H1
        {"O", {0.00, 0.00, 0.00}},  // O
        {"H", {-1.43, -1.10, 0.00}} // H2
    };

    // Allocate displacements in cartesian coordinates
    vec dx{linalg::zeros<vec>(3 * mol.size())};

    // Compute displacements
    for (std::size_t i{0}; i < 3; i++) {
      vec3 v{R[i] * mol[i].position - mol[i].position};

      dx(3 * i + 0) = v(0);
      dx(3 * i + 1) = v(1);
      dx(3 * i + 2) = v(2);
    }

    // Molecular connectivity
    mat dd{distances<vec3, mat>(mol)};
    UGraph adj{adjacency_matrix(dd, mol)};
    mat dist{distance_matrix<mat>(adj)};

    // Compute bonds
    const auto B = bonds(dist, mol);
    REQUIRE(B.size() == 2);

    // Compute angles
    const auto A = angles(dist, mol);
    REQUIRE(A.size() == 1);

    // Compute Wilson B matrix for H2O analytically
    mat Bwilson =
        wilson_matrix<vec3, vec, mat>(to_cartesian<vec3, vec>(mol), B, A);
    REQUIRE(linalg::size(Bwilson) == 27);
    INFO("Wilson B matrix (analytical):\n" << Bwilson);

    // Compute Wilson B matrix for H2O numerically
    const mat BwilsonN = wilson_matrix_numerical<vec3, vec, mat>(
        to_cartesian<vec3, vec>(mol), B, A);
    REQUIRE(linalg::size(BwilsonN) == 27);
    INFO("Wilson B matrix (numerical):\n" << BwilsonN);

    // Check analytical and numerical Wilson matrices are the same
    for (std::size_t i{0}; i < 27; i++)
      REQUIRE(Bwilson(i) == Approx(BwilsonN(i)).margin(1e-6));

    INFO("Compute displacements in internal coordinates");
    const vec displacement = Bwilson * dx;
    CAPTURE(displacement);

    INFO("Check bond change");
    REQUIRE(displacement(0) == Approx(0).margin(1e-4));
    REQUIRE(displacement(1) == Approx(0).margin(1e-4));

    INFO("Check angle change");
    REQUIRE(displacement(2) == Approx(2 * angle_rad).margin(1e-3));
  }

  SECTION("H2O2 torsion") {

    const double angle(1.0);
    const double angle_rad(angle / 180. * constants::pi);

    const mat R{{cos(angle_rad), -sin(angle_rad), 0},
                {sin(angle_rad), cos(angle_rad), 0},
                {0, 0, 1}};

    molecule::Molecule<vec3> molecule{
        {"H", {0.000, 0.947, -0.079}}, // H1
        {"O", {0.000, 0.000, 0.000}},  // O1
        {"O", {0.000, 0.000, 1.474}},  // O2
        {"H", {-0.854, -0.407, 1.553}} // H2
    };
    molecule::multiply_positions(molecule, conversion::angstrom_to_bohr);

    // Allocate displacements in cartesian coordinates
    vec dx{linalg::zeros<vec>(3 * molecule.size())};

    // Compute transformation for H1 rotation
    vec3 v{R * molecule[0].position - molecule[0].position};
    dx(0) = v(0);
    dx(1) = v(1);
    dx(2) = v(2);

    // Compute old dihedral (before rotation)
    double d_old{dihedral<vec3>(molecule[0].position,
                                molecule[1].position,
                                molecule[2].position,
                                molecule[3].position)};

    // Compute new dihedral angle (after rotation)
    double d_new{dihedral<vec3>(R * molecule[0].position,
                                molecule[1].position,
                                molecule[2].position,
                                molecule[3].position)};

    // Compute dihedral variation
    double d_diff{d_new - d_old};

    // Compute interatomic distances
    mat dd{distances<vec3, mat>(molecule)};
    UGraph adj{adjacency_matrix(dd, molecule)};
    mat dist{distance_matrix<mat>(adj)};

    // Compute bonds, angles and dihedrals
    const std::vector<Bond> B{bonds(dist, molecule)};
    CHECK(B.size() == 3);
    const std::vector<Angle> A{angles(dist, molecule)};
    CHECK(A.size() == 2);
    const std::vector<Dihedral> D{dihedrals(dist, molecule)};
    CHECK(D.size() == 1);

    // Compute Wilson's B matrix
    const mat Bwilson = wilson_matrix<vec3, vec, mat>(
        to_cartesian<vec3, vec>(molecule), B, A, D);
    REQUIRE(linalg::size(Bwilson) == 72);
    INFO("Wilson B matrix (analytical):\n" << Bwilson);

    // Compute Wilson's B matrix
    const mat BwilsonN = wilson_matrix_numerical<vec3, vec, mat>(
        to_cartesian<vec3, vec>(molecule), B, A, D);
    REQUIRE(linalg::size(BwilsonN) == 72);
    INFO("Wilson B matrix (numerical):\n" << BwilsonN);

    // Check analytical and numerical Wilson matrices are the same
    INFO("Check Analytical vs Numerical Wilson B matrix");
    for (std::size_t i{0}; i < 72; i++) {
      REQUIRE(Bwilson(i) == Approx(BwilsonN(i)).margin(1e-6));
    }

    SECTION("Analytical displacements") {
      // Compute displacement in internal coordinates
      const vec displacement{Bwilson * dx};
      INFO("Displacement (analytical):\n" << displacement);

      INFO("Bonds do not change");
      REQUIRE(displacement(0) == Approx(0).margin(1e-3));
      REQUIRE(displacement(1) == Approx(0).margin(1e-3));
      REQUIRE(displacement(2) == Approx(0).margin(1e-3));

      INFO("Angles do not change");
      REQUIRE(displacement(3) == Approx(0).margin(1e-4));
      REQUIRE(displacement(4) == Approx(0).margin(1e-4));

      INFO("Dihedral should change");
      REQUIRE(displacement(5) == Approx(d_diff).margin(1e-4));
    }

    SECTION("Numerical displacements") {
      // Compute displacement in internal coordinates
      const vec displacement{BwilsonN * dx};
      INFO("Displacement (numerical):\n" << displacement);

      INFO("Bonds do not change");
      REQUIRE(displacement(0) == Approx(0).margin(1e-3));
      REQUIRE(displacement(1) == Approx(0).margin(1e-3));
      REQUIRE(displacement(2) == Approx(0).margin(1e-3));

      INFO("Angles do not change");
      REQUIRE(displacement(3) == Approx(0).margin(1e-4));
      REQUIRE(displacement(4) == Approx(0).margin(1e-4));

      INFO("Dihedral should change");
      REQUIRE(displacement(5) == Approx(d_diff).margin(1e-4));
    }
  }
}

TEST_CASE("Wilson B matrix for formalydehyde", "[wilson]") {
  using namespace connectivity;
  using namespace molecule;
  using namespace tools;
  using namespace wilson;
  using namespace io;

  const auto mol = load_xyz<vec3>(config::molecules_dir + "formaldehyde.xyz");

  // Compute interatomic distances
  mat dd{distances<vec3, mat>(mol)};
  UGraph adj{adjacency_matrix(dd, mol)};
  mat dist{distance_matrix<mat>(adj)};

  // Compute bonds
  std::vector<Bond> B{bonds(dist, mol)};
  std::vector<Angle> A{angles(dist, mol)};
  std::vector<Dihedral> D{dihedrals(dist, mol)};
  std::vector<LinearAngle<vec3>> LA{linear_angles(dist, mol)};
  std::vector<OutOfPlaneBend> OOPB{out_of_plane_bends(dist, mol)};

  const auto q = connectivity::cartesian_to_irc<vec3, vec>(
      to_cartesian<vec3, vec>(mol), B, A, D, LA, OOPB);
  CAPTURE(q);

  const mat wilson_b_analytical =
      wilson_matrix<vec3, vec, mat>(to_cartesian<vec3, vec>(mol), B, A, D, LA, OOPB);
  CAPTURE(wilson_b_analytical);

  const mat wilson_b_numerical = wilson_matrix_numerical<vec3, vec, mat>(
      to_cartesian<vec3, vec>(mol), B, A, D, LA, OOPB);
  CAPTURE(wilson_b_numerical);

  REQUIRE(linalg::size(wilson_b_analytical) ==
          linalg::size(wilson_b_numerical));


  const std::size_t n = linalg::size(wilson_b_analytical);
  for (std::size_t i{0}; i < n; i++) {
    CAPTURE(i);
    REQUIRE(wilson_b_analytical(i) ==
            Approx(wilson_b_numerical(i)).margin(1e-5));
  }
}

TEST_CASE("Wilson B matrix for water dimer", "[wilson]") {
  using namespace connectivity;
  using namespace molecule;
  using namespace tools;
  using namespace wilson;
  using namespace io;

  const auto mol = load_xyz<vec3>(config::molecules_dir + "water_dimer_2.xyz");

  // Compute interatomic distances
  mat dd{distances<vec3, mat>(mol)};
  UGraph adj{adjacency_matrix(dd, mol)};
  mat dist{distance_matrix<mat>(adj)};

  // Compute bonds
  std::vector<Bond> B{bonds(dist, mol)};
  std::vector<Angle> A{angles(dist, mol)};
  std::vector<Dihedral> D{dihedrals(dist, mol)};
  std::vector<LinearAngle<vec3>> LA{linear_angles(dist, mol)};
  std::vector<OutOfPlaneBend> OOPB{out_of_plane_bends(dist, mol)};

  const auto q = connectivity::cartesian_to_irc<vec3, vec>(
      to_cartesian<vec3, vec>(mol), B, A, D, LA, OOPB);
  CAPTURE(q);

  const mat wilson_b_analytical =
      wilson_matrix<vec3, vec, mat>(to_cartesian<vec3, vec>(mol), B, A, D, LA, OOPB);
  CAPTURE(wilson_b_analytical);

  const mat wilson_b_numerical = wilson_matrix_numerical<vec3, vec, mat>(
      to_cartesian<vec3, vec>(mol), B, A, D, LA, OOPB);
  CAPTURE(wilson_b_numerical);

  REQUIRE(linalg::size(wilson_b_analytical) ==
          linalg::size(wilson_b_numerical));

  const std::size_t n = linalg::size(wilson_b_analytical);
  for (std::size_t i{0}; i < n; i++) {
    REQUIRE(wilson_b_analytical(i) ==
            Approx(wilson_b_numerical(i)).margin(1e-5));
  }
}


TEST_CASE("Linear angle gradient", "[wilson]") {
    using namespace connectivity;
    using namespace molecule;
    using namespace tools;
    using namespace wilson;
    using namespace io;

    const auto mol = load_xyz<vec3>(config::molecules_dir + "hcn.xyz");
    REQUIRE(mol.size() == 3);

    // Compute interatomic distances
    mat dd{distances<vec3, mat>(mol)};
    UGraph adj{adjacency_matrix(dd, mol)};
    mat dist{distance_matrix<mat>(adj)};

    // Compute bonds
    std::vector<Bond> B{bonds(dist, mol)};
    std::vector<Angle> A{angles(dist, mol)};
    std::vector<Dihedral> D{dihedrals(dist, mol)};
    std::vector<LinearAngle<vec3>> LA{linear_angles(dist, mol)};
    std::vector<OutOfPlaneBend> OOPB{out_of_plane_bends(dist, mol)};
    REQUIRE(LA.size() == 2);

    const auto q = connectivity::cartesian_to_irc<vec3, vec>(
        to_cartesian<vec3, vec>(mol), B, A, D, LA, OOPB);
    CAPTURE(q);

    const mat wilson_b_analytical =
        wilson_matrix<vec3, vec, mat>(to_cartesian<vec3, vec>(mol), B, A, D, LA, OOPB);
    CAPTURE(wilson_b_analytical);

    const mat wilson_b_numerical = wilson_matrix_numerical<vec3, vec, mat>(
        to_cartesian<vec3, vec>(mol), B, A, D, LA, OOPB);
    CAPTURE(wilson_b_numerical);

    REQUIRE(linalg::size(wilson_b_analytical) ==
            linalg::size(wilson_b_numerical));

    const std::size_t n = linalg::size(wilson_b_analytical);
    for (std::size_t i{0}; i < n; i++) {
      REQUIRE(wilson_b_analytical(i) ==
              Approx(wilson_b_numerical(i)).margin(1e-5));
    }
  }
