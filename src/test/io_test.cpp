#include "../../include/catch/catch.hpp"

#include "libirc/io.h"

#include "config.h"
#include "libirc/molecule.h"
#include "libirc/conversion.h"

#include <iostream>
#include <stdexcept>

#ifdef HAVE_ARMA
#include <armadillo>
using arma::vec3;
using arma::vec;
using arma::mat;
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
}

TEST_CASE("Print"){
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

  // Print bonds to std::cout
  print_bonds<vec3,vec>( to_cartesian<vec3,vec>(mol), B );

  // Compute angles
  std::vector<Angle> A{angles(dist, predecessors, mol)};

  // Print angles to std::cout
  print_angles<vec3,vec>( to_cartesian<vec3,vec>(mol), A );

  // Compute dihedral angles
  std::vector<Dihedral> D{dihedrals(dist, predecessors, mol)};

  // Print dihedrals to std::cout
  print_dihedrals<vec3,vec>( to_cartesian<vec3,vec>(mol), D );
}