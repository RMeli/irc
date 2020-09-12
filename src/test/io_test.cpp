#include "catch.hpp"

#include "libirc/io.h"

#include "config.h"

#ifdef HAVE_ARMA
#include <armadillo>
using arma::mat;
using arma::vec;
using arma::vec3;
#elif HAVE_EIGEN3
#include <Eigen/Dense>
using vec3 = Eigen::Vector3d;
using vec = Eigen::VectorXd;
using mat = Eigen::MatrixXd;

#else
#error
#endif

TEST_CASE("File not found") {
  using irc::io::load_xyz;

  REQUIRE_THROWS_AS(load_xyz<vec3>(irc::config::molecules_dir + "ABC.xyz"),
                    std::runtime_error);

  REQUIRE_NOTHROW(load_xyz<vec3>(irc::config::molecules_dir + "caffeine.xyz"));
}

TEST_CASE("Print molecule") {
  using namespace irc;
  using namespace io;
  using namespace connectivity;
  using namespace molecule;
  using namespace tools;

  // Load toluene molecule
  const auto mol = load_xyz<vec3>(config::molecules_dir + "benzene_dimer.xyz");

  // Compute interatomic distance for formaldehyde molecule
  const mat dd{distances<vec3, mat>(mol)};

  // Build graph based on the adjacency matrix
  const UGraph adj{adjacency_matrix(dd, mol)};

  // Compute distance matrix and predecessor matrix
  mat dist{distance_matrix<mat>(adj)};

  // Compute bonds
  const std::vector<Bond> B{bonds(dist, mol)};

  // Print bonds to std::cout
  print_bonds<vec3, vec>(to_cartesian<vec3, vec>(mol), B);

  // Compute angles
  const std::vector<Angle> A{angles(dist, mol)};

  // Print angles to std::cout
  print_angles<vec3, vec>(to_cartesian<vec3, vec>(mol), A);

  // Compute dihedral angles
  const std::vector<Dihedral> D{dihedrals(dist, mol)};

  // Print dihedrals to std::cout
  print_dihedrals<vec3, vec>(to_cartesian<vec3, vec>(mol), D);

  // Compute out of plane bends
  const std::vector<OutOfPlaneBend> OOPB{out_of_plane_bends(dist, mol)};

  // Print dihedrals to std::cout
  print_out_of_plane_bends<vec3, vec>(to_cartesian<vec3, vec>(mol), OOPB);
}
