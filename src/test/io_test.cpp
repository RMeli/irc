#include "../../include/catch/catch.hpp"

#include "libirc/io.h"

#include "config.h"
#include "libirc/conversion.h"
#include "libirc/molecule.h"

#include <iostream>
#include <stdexcept>

#ifdef HAVE_ARMA
#include <armadillo>
using arma::mat;
using arma::vec;
using arma::vec3;
#elif HAVE_EIGEN3
#include <Eigen3/Eigen/Dense>
using vec3 = Eigen::Vector3d;
using vec = Eigen::VectorXd;
using mat = Eigen::MatrixXd;

#else
#error
#endif

using namespace irc;

TEST_CASE("File not found") {
  using namespace io;
  using namespace molecule;

  bool exception{false};

  try {
    Molecule<vec3> mol{load_xyz<vec3>(config::molecules_dir + "ABC.xyz")};
  } catch (const std::runtime_error &e) {
    exception = true;
  }

  REQUIRE(exception == true);
}

TEST_CASE("Loading XYZ file") {

  using namespace io;
  using namespace molecule;

  Molecule<vec3> mol{load_xyz<vec3>(config::molecules_dir + "caffeine.xyz")};

  std::cout << mol << std::endl;

  try {
    Molecule<vec3> mol{load_xyz<vec3>(config::molecules_dir + "foo.xyz")};
  } catch (const std::runtime_error &e) {
  }
}

TEST_CASE("Print toluene") {
  using namespace io;

  using namespace connectivity;
  using namespace molecule;
  using namespace tools;

  // Load toluene molecule
  Molecule<vec3> mol{load_xyz<vec3>(config::molecules_dir + "toluene.xyz")};

  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(mol, conversion::angstrom_to_bohr);

  // Compute interatomic distance for formaldehyde molecule
  mat dd{distances<vec3, mat>(mol)};

  // Build graph based on the adjacency matrix
  UGraph adj{adjacency_matrix(dd, mol)};

  // Compute distance matrix and predecessor matrix
  mat dist, predecessors;
  std::tie(dist, predecessors) = distance_matrix<mat>(adj);

  // Compute bonds
  std::vector<Bond> B{bonds(dist, mol)};

  // Chek number of bonds
  REQUIRE( B.size() == 15 );

  // Print bonds to std::cout
  print_bonds<vec3, vec>(to_cartesian<vec3, vec>(mol), B);

  // Compute angles
  std::vector<Angle> A{angles(dist, predecessors, mol)};

  // Check number of angles
  REQUIRE( A.size() == 24 );

  // Print angles to std::cout
  print_angles<vec3, vec>(to_cartesian<vec3, vec>(mol), A);

  // Compute dihedral angles
  std::vector<Dihedral> D{dihedrals(dist, predecessors, mol)};

  // Check number of dihedral angles
  REQUIRE( D.size() == 30 );

  // Print dihedrals to std::cout
  print_dihedrals<vec3, vec>(to_cartesian<vec3, vec>(mol), D);
}


TEST_CASE("Print caffeine") {
  using namespace io;

  using namespace connectivity;
  using namespace molecule;
  using namespace tools;

  // Load toluene molecule
  Molecule<vec3> mol{load_xyz<vec3>(config::molecules_dir + "caffeine.xyz")};

  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(mol, conversion::angstrom_to_bohr);

  // Compute interatomic distance for formaldehyde molecule
  mat dd{distances<vec3, mat>(mol)};

  // Build graph based on the adjacency matrix
  UGraph adj{adjacency_matrix(dd, mol)};

  // Compute distance matrix and predecessor matrix
  mat dist, predecessors;
  std::tie(dist, predecessors) = distance_matrix<mat>(adj);

  // Compute bonds
  std::vector<Bond> B{bonds(dist, mol)};

  // Print bonds to std::cout
  print_bonds<vec3, vec>(to_cartesian<vec3, vec>(mol), B);

  // Chek number of bonds
  REQUIRE(B.size() == 25);

  // Compute angles
  std::vector<Angle> A{angles(dist, predecessors, mol)};

  // Print angles to std::cout
  print_angles<vec3, vec>(to_cartesian<vec3, vec>(mol), A);

  // Chek number of angles
  REQUIRE(A.size() == 43);

  // Compute dihedral angles
  std::vector<Dihedral> D{dihedrals(dist, predecessors, mol)};

  // Print dihedrals to std::cout
  print_dihedrals<vec3, vec>(to_cartesian<vec3, vec>(mol), D);

  // Chek number of dihedral angles
  REQUIRE(D.size() == 54);
}