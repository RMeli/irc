#include "../../include/catch/catch.hpp"

#include "libirc/wilson.h"

#include "libirc/atom.h"
#include "libirc/connectivity.h"
#include "libirc/conversion.h"
#include "libirc/io.h"
#include "libirc/molecule.h"

#include "config.h"

#include <cmath>
#include <iostream>

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

TEST_CASE("Wilson B matrix", "[wilson]") {

  bool verbose{true};

  using namespace std;

  using namespace connectivity;
  using namespace molecule;
  using namespace wilson;
  using namespace tools;

  SECTION("H2 stretching") {

    // Define molecule
    Molecule<vec3> mol{{"H", {0., 0., 0.}}, {"H", {1., 0., 0.}}};

    // Compute interatomic distances
    mat dd{distances<vec3, mat>(mol)};

    // Compute adjacency matrix (graph)
    UGraph adj{adjacency_matrix(dd, mol)};

    // Compute distance matrix and predecessor matrix
    mat dist, pred;
    tie(dist, pred) = distance_matrix<mat>(adj);

    // Compute bonds
    vector<Bond> B{bonds(dist, mol)};

    REQUIRE(B.size() == 1);

    // Compute Wilson B matrix for H2 analytically
    mat Bwilson =
        wilson_matrix<vec3, vec, mat>(to_cartesian<vec3, vec>(mol), B);

    // Check Wilson B matrix size
    REQUIRE(linalg::size(Bwilson) == 6);

    if (verbose) {
      cout << "Wilson B matrix (analytical):" << endl;
      cout << Bwilson << endl;
    }

    // Compute Wilson B matrix for H2 numerically
    mat BwilsonN = wilson_matrix_numerical<vec3, vec, mat>(
        to_cartesian<vec3, vec>(mol), B);

    // Check Wilson B matrix size
    REQUIRE(linalg::size(BwilsonN) == 6);

    if (verbose) {
      cout << "Wilson B matrix (numerical):" << endl;
      cout << BwilsonN << endl;
    }

    // Check analytical and numerical Wilson matrices are the same
    SECTION("Analytical vs Numerical") {
      for (size_t i{0}; i < 6; i++) {
        Approx target{BwilsonN(i)};
        target.margin(1e-6);

        REQUIRE(Bwilson(i) == target);
      }
    }

    // Set bond stretching
    double d{0.01};
    vec dx{-d, 0.00, 0.00, d, 0.00, 0.00};

    SECTION("Bond change (analytical)") {
      // Compute transformation in IRC
      vec transformation{Bwilson * dx};

      if (verbose) {
        cout << "\nTransformation (analytical)" << endl;
        cout << Bwilson * dx << endl;
      }

      Approx target{2 * d};
      target.margin(1e-5);

      REQUIRE(transformation(0) == target);
    }

    SECTION("Bond change (numerical)") {
      // Compute transformation in IRC
      vec transformation{BwilsonN * dx};

      if (verbose) {
        cout << "\nTransformation (numerical)" << endl;
        cout << Bwilson * dx << endl;
      }

      Approx target{2 * d};
      target.margin(1e-5);

      REQUIRE(transformation(0) == target);
    }
  } // H2 stretching

  SECTION("H2O bending") {

    double angle(0.5);
    double angle_rad(angle / 180. * constants::pi);

    std::vector<mat> R{
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

    Molecule<vec3> mol{
        {"H", {1.43, -1.10, 0.00}}, // H1
        {"O", {0.00, 0.00, 0.00}},  // O
        {"H", {-1.43, -1.10, 0.00}} // H2
    };

    // Allocate displacements in cartesian coordinates
    vec dx{linalg::zeros<vec>(3 * mol.size())};

    // Compute displacements
    for (size_t i{0}; i < 3; i++) {
      vec3 v{R[i] * mol[i].position - mol[i].position};

      dx(3 * i + 0) = v(0);
      dx(3 * i + 1) = v(1);
      dx(3 * i + 2) = v(2);
    }

    // Compute interatomic distances
    mat dd{distances<vec3, mat>(mol)};

    // Compute adjacency matrix (graph)
    UGraph adj{adjacency_matrix(dd, mol)};

    // Compute distance matrix and predecessor matrix
    mat dist, pred;
    tie(dist, pred) = distance_matrix<mat>(adj);

    // Compute bonds
    vector<Bond> B{bonds(dist, mol)};

    // Check number of bonds
    REQUIRE(B.size() == 2);

    // Compute bonds
    vector<Angle> A{angles(dist, pred, mol)};

    REQUIRE(A.size() == 1);

    // Compute Wilson B matrix for H2O analytically
    mat Bwilson =
        wilson_matrix<vec3, vec, mat>(to_cartesian<vec3, vec>(mol), B, A);

    // Check Wilson B matrix size
    REQUIRE(linalg::size(Bwilson) == 27);

    if (verbose) {
      // Print Wilson B matrix
      cout << "\nWilson B matrix (analytical):" << endl;
      cout << Bwilson << endl;
    }

    // Compute Wilson B matrix for H2O numerically
    mat BwilsonN = wilson_matrix_numerical<vec3, vec, mat>(
        to_cartesian<vec3, vec>(mol), B, A);

    // Check Wilson B matrix size
    REQUIRE(linalg::size(BwilsonN) == 27);

    if (verbose) {
      // Print Wilson B matrix
      cout << "\nWilson B matrix (numerical):" << endl;
      cout << BwilsonN << endl;
    }

    // Check analytical and numerical Wilson matrices are the same
    SECTION("Analytical vs Numerical") {
      for (size_t i{0}; i < 27; i++) {
        Approx target{BwilsonN(i)};
        target.margin(1e-6);

        REQUIRE(Bwilson(i) == target);
      }
    }

    SECTION("Analytical displacement") {
      // Compute displacement in internal coordinates
      vec displacement{Bwilson * dx};
      cout << "\nDisplacement:" << endl;
      cout << displacement << endl;

      // Check that bonds do not change
      SECTION("Bond change") {
        Approx target{0};

        target.margin(1e-4);

        REQUIRE(displacement(0) == target);
        REQUIRE(displacement(1) == target);
      }

      // Check change in angle
      SECTION("Angle change") {
        Approx target{2 * angle};

        target.margin(1e-3);

        REQUIRE(displacement(2) * 180 / constants::pi == target);
      }
    }
  }

  SECTION("H2O2 torsion") {

    double angle(1.0);
    double angle_rad(angle / 180. * constants::pi);

    mat R{{cos(angle_rad), -sin(angle_rad), 0},
          {sin(angle_rad), cos(angle_rad), 0},
          {0, 0, 1}};

    molecule::Molecule<vec3> molecule{
        {"H", {0.000, 0.947, -0.079}}, // H1
        {"O", {0.000, 0.000, 0.000}},  // O1
        {"O", {0.000, 0.000, 1.474}},  // O2
        {"H", {-0.854, -0.407, 1.553}} // H2
    };

    // Transform molecular coordinates from angstrom to bohr
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

    // Compute adjacency matrix (graph)
    UGraph adj{adjacency_matrix(dd, molecule)};

    // Compute distance matrix and predecessor matrix
    mat dist, predecessors;
    tie(dist, predecessors) = distance_matrix<mat>(adj);

    // Compute bonds
    vector<Bond> B{bonds(dist, molecule)};

    // Check number of bonds
    REQUIRE(B.size() == 3);

    // Compute angles
    vector<Angle> A{angles(dist, predecessors, molecule)};

    // Check number of angles
    REQUIRE(A.size() == 2);

    // Compute dihedrals
    vector<Dihedral> D{dihedrals(dist, predecessors, molecule)};

    // Check number of dihedrals
    REQUIRE(D.size() == 1);

    // Compute Wilson's B matrix
    mat Bwilson = wilson_matrix<vec3, vec, mat>(
        to_cartesian<vec3, vec>(molecule), B, A, D);

    REQUIRE(linalg::size(Bwilson) == 72);

    if (verbose) {
      // Print Wilson B matrix
      cout << "\nWilson B matrix (analytical):" << endl;
      cout << Bwilson << endl;
    }

    // Compute Wilson's B matrix
    mat BwilsonN = wilson_matrix_numerical<vec3, vec, mat>(
        to_cartesian<vec3, vec>(molecule), B, A, D);

    REQUIRE(linalg::size(BwilsonN) == 72);

    if (verbose) {
      // Print Wilson B matrix
      cout << "\nWilson B matrix (numerical):" << endl;
      cout << BwilsonN << endl;
    }

    // Check analytical and numerical Wilson matrices are the same
    SECTION("Analytical vs Numerical") {
      for (size_t i{0}; i < 72; i++) {
        Approx target{BwilsonN(i)};
        target.margin(1e-6);

        REQUIRE(Bwilson(i) == target);
      }
    }

    SECTION("Analytical displacements") {
      // Compute displacement in internal coordinates
      vec displacement{Bwilson * dx};

      // Print displacements
      if (verbose) {
        cout << "\nDisplacement (analytical):" << endl;
        cout << displacement << endl;
      }

      // Check that bonds do not change
      SECTION("Bonds change") {
        Approx target{0};

        target.margin(1e-3);

        REQUIRE(displacement(0) == target);
        REQUIRE(displacement(1) == target);
        REQUIRE(displacement(2) == target);
      }

      // Check change in angle
      SECTION("Angles change") {
        Approx target{0};

        target.margin(1e-4);

        REQUIRE(displacement(3) == target);
        REQUIRE(displacement(4) == target);
      }

      // Check change in angle
      SECTION("Dihedral change") {
        Approx target{d_diff};

        target.margin(1e-4);

        REQUIRE(displacement(5) == target);
      }
    }

    SECTION("Numerical displacements") {
      // Compute displacement in internal coordinates
      vec displacement{BwilsonN * dx};

      // Print displacements
      if (verbose) {
        cout << "\nDisplacement (numerical):" << endl;
        cout << displacement << endl;
      }

      // Check that bonds do not change
      SECTION("Bonds change") {
        Approx target{0};

        target.margin(1e-3);

        REQUIRE(displacement(0) == target);
        REQUIRE(displacement(1) == target);
        REQUIRE(displacement(2) == target);
      }

      // Check change in angle
      SECTION("Angles change") {
        Approx target{0};

        target.margin(1e-4);

        REQUIRE(displacement(3) == target);
        REQUIRE(displacement(4) == target);
      }

      // Check change in angle
      SECTION("Dihedral change") {
        Approx target{d_diff};

        target.margin(1e-4);

        REQUIRE(displacement(5) == target);
      }
    }
  }
}

TEST_CASE("Wilson") {
  using namespace std;

  using namespace connectivity;
  using namespace molecule;
  using namespace tools;
  using namespace wilson;

  bool verbose{true};

  Molecule<vec3> mol{
      io::load_xyz<vec3>(config::molecules_dir + "water_dimer_2.xyz")};

  // Transform molecule to bohr
  multiply_positions(mol, conversion::angstrom_to_bohr);

  // Compute interatomic distances
  mat dd{distances<vec3, mat>(mol)};

  // Compute adjacency matrix (graph)
  UGraph adj{adjacency_matrix(dd, mol)};

  // Compute distance matrix and predecessor matrix
  mat dist, predecessors;
  tie(dist, predecessors) = distance_matrix<mat>(adj);

  // Compute bonds
  vector<Bond> B{bonds(dist, mol)};

  // Compute angles
  vector<Angle> A{angles(dist, predecessors, mol)};

  // Compute dihedrals
  vector<Dihedral> D{dihedrals(dist, predecessors, mol)};

  cout << "q =\n"
       << connectivity::cartesian_to_irc<vec3, vec>(to_cartesian<vec3, vec>(mol), B, A, D)
       << endl;

  // Return Wilson's B matrix (analytical)
  mat WB{wilson_matrix<vec3, vec, mat>(to_cartesian<vec3, vec>(mol), B, A, D)};

  if (verbose) {
    cout << "B (analytical) =\n" << WB << endl;
  }

  // Return Wilson's B matrix (numerical)
  mat WBN{wilson_matrix_numerical<vec3, vec, mat>(
      to_cartesian<vec3, vec>(mol), B, A, D)};

  if (verbose) {
    cout << "B (numerical) =\n" << WBN << endl;
  }

  size_t n{linalg::size(WB)};
  for (size_t i{0}; i < n; i++) {
    Approx target(WBN(i));
    target.margin(1e-5);

    REQUIRE(WB(i) == target);
  }
}