#include "../../include/catch/catch.hpp"

#include "config.h"

#include "libirc/irc.h"

#include "libirc/conversion.h"
#include "libirc/molecule.h"
#include "libirc/io.h"

#include <iostream>

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

    // Define formaldehyde molecule (CH2O)
    Molecule<vec3> molecule{{"C", {0.000000, 0.000000, -0.537500}},
                            {"O", {0.000000, 0.000000, 0.662500}},
                            {"H", {0.000000, 0.866025, -1.037500}},
                            {"H", {0.000000, -0.866025, -1.037500}}};

    // Transform molecular coordinates from angstrom to bohr
    multiply_positions(molecule, angstrom_to_bohr);

    // Build internal reaction coordinates
    // Add O-H bonds: (1,2)
    // Add H-H bond: (2,3)
    // Add H-O-H angle: (2,1,3)
    // Add H-H-O-C dihedral: (3,2,1,0)
    IRC<vec3, vec, mat> irc(
        molecule, {{1, 2}, {2, 3}}, {{2, 1, 3}}, {{3, 2, 1, 0}});

    // Compute internal coordinates
    vec q_irc{
        irc.cartesian_to_irc(molecule::to_cartesian<vec3, vec>(molecule))};

    // Check size (3+2 bonds, 3+1 angles, 0+1 dihedrals)
    REQUIRE(linalg::size(q_irc) == 10);

    // Check manually added O-H bond
    SECTION("Manually added O-H bond") {
      Approx target(q_irc(3));
      target.margin(1e-6);

      // Compute O-H distance
      double d{distance(molecule[1].position, molecule[2].position)};

      REQUIRE(d == target);
    }

    // Check manually added H-H bond
    SECTION("Manually added H-H bond") {
      Approx target(q_irc(4));
      target.margin(1e-6);

      // Compute H-H disance
      double b{distance(molecule[2].position, molecule[3].position)};

      REQUIRE(b == target);
    }

    // Check manually added H-O-H angle
    SECTION("Manually added H-O-H angle") {
      Approx target(q_irc(8));
      target.margin(1e-6);

      // Compute H-O-H angle
      double a{angle(
          molecule[2].position, molecule[1].position, molecule[3].position)};

      REQUIRE(a == target);
    }
    
    // Check manually added H-H-O-C dihedral angle
    SECTION("Manually added H-H-O-C dihedral angle"){
      Approx target(q_irc(9));
      target.margin(1e-6);
  
      // Compute H-H-O-C dihedral angle
      double d{dihedral(
          molecule[3].position, molecule[2].position, molecule[1].position, molecule[0].position)};
  
      REQUIRE(d == target);
    }
  }

  SECTION("Initial hessian") {

    // Define formaldehyde molecule (CH2O)
    Molecule<vec3> molecule{{"C", {0.000000, 0.000000, -0.537500}},
                            {"O", {0.000000, 0.000000, 0.662500}},
                            {"H", {0.000000, 0.866025, -1.037500}},
                            {"H", {0.000000, -0.866025, -1.037500}}};

    // Transform molecular coordinates from angstrom to bohr
    multiply_positions(molecule, angstrom_to_bohr);

    // Build internal reaction coordinates
    // (Manually added dihedral to increase to increase coverage)
    IRC<vec3, vec, mat> irc(molecule,{},{},{{0,1,2,3}});

    // Compute initial hessian
    mat iH0{irc.projected_initial_hessian_inv()};

    // Project Hessian again
    mat iH{irc.projected_hessian_inv(iH0)};

    // Check sizes
    REQUIRE(linalg::size(iH0) == linalg::size(iH));

    // Chech that second projection has no effect
    size_t n{linalg::size(iH0)};
    for (size_t i{0}; i < n; i++) {
      Approx target(iH0(i));
      target.margin(1e-6);

      REQUIRE(iH(i) == target);
    }
  }
  
  SECTION("IRC to Cartesian"){
    
    // Define formaldehyde molecule (CH2O)
    Molecule<vec3> molecule{io::load_xyz<vec3>(config::molecules_dir + "ethanol.xyz")};
  
    // Transform molecular coordinates from angstrom to bohr
    multiply_positions(molecule, angstrom_to_bohr);
  
    // Build internal reaction coordinates
    IRC<vec3, vec, mat> irc(molecule);
    
    // Get cartesian coordinates
    vec x_c{ to_cartesian<vec3,vec>(molecule) };
    
    // Compute internal redundant coordinates
    vec q_irc{ irc.cartesian_to_irc(x_c) };
    
    // Define no displacement in IRC
    vec dq{linalg::zeros<vec>(linalg::size(q_irc))};
    
    // Compute cartesian coordinate from IRC
    vec x_c_from_irc{ irc.irc_to_cartesian(q_irc, dq, x_c) };
    
    size_t n{linalg::size(x_c)};
    for(size_t i{0}; i < n; i++){
      Approx target( x_c(i) );
      target.margin(1e-6);
      
      REQUIRE(x_c_from_irc(i) == target);
    }
  }
}