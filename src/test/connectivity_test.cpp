#include "../../include/catch/catch.hpp"

#include "libirc/connectivity.h"

#include "libirc/atom.h"
#include "libirc/molecule.h"
#include "config.h"
#include "libirc/io.h"
#include "libirc/conversion.h"


#include <iostream>

#ifdef HAVE_ARMA
#include <armadillo>
using vec3 = arma::vec3;
using vec = arma::vec;
using mat = arma::mat;

template <typename T>
using Mat = arma::Mat<T>;
#elif HAVE_EIGEN3
#include <Eigen3/Eigen/Dense>
using vec3 = Eigen::Vector3d;
using vec = Eigen::VectorXd;
using mat = Eigen::MatrixXd;

template <typename T>
using Mat = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;
#else
#error
#endif

using namespace irc;

TEST_CASE("Distance, angle and dihedral"){
  using namespace connectivity;
  
  // Define three points in 3D space
  vec3 p1{ 0.00,  0.00, -0.25};
  vec3 p2{ 0.00,  0.00,  1.50};
  vec3 p3{ 0.00,  1.00,  1.50};
  vec3 p4{ 5.00,  1.00,  1.50};
  
  SECTION("Distance"){
    Approx target{1.75};
    
    target.margin(1e-12);
    
    REQUIRE( distance(p1, p2) == target);
  }
  
  SECTION("Angle"){
    Approx target{90};
    
    target.margin(1e-12);
    
    REQUIRE( angle(p1, p2, p3) == target);
  }
  
  SECTION("Dihedral"){
    Approx target{-90};
  
    target.margin(1e-12);
  
    REQUIRE( dihedral(p1, p2, p3, p4) == target);
  }
  
}

// Regular bond
TEST_CASE("Connectivity for compressed H2"){
  using namespace std;
  
  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;
  
  double d{0.5};
  
  // Define compressed H2 molecule
  Molecule<vec3> molecule{
      {"H", {0.0, 0.0, d  }},
      {"H", {0.0, 0.0, 0.0}}
  };
  
  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);
  
  // Compute interatomic distance for compressed H2
  mat dd{ distances<vec3, mat>(molecule) };
  
  // Compute adjacency graph for compressed H2
  UGraph adj{ adjacency_matrix(dd, molecule) };
  
  // Compute distance and predecessor matrices for compressed H2
  Mat<int> dist, predecessors;
  std::tie(dist, predecessors) = distance_matrix<Mat<int>>(adj);
  
  SECTION("Bond"){
    // Compute bond
    std::vector<Bond> B{ bonds(dist, molecule) };
    
    // Check number of bonds
    REQUIRE( B.size() == 1);
    
    // Compute IRC
    vec q{irc_from_bad<vec3,vec>(to_cartesian<vec3,vec>(molecule), B, {}, {})};
    
    // Check number of IRC
    REQUIRE( linalg::size<vec>(q) == 1);
    
    // Check bond length
    Approx bb(d * angstrom_to_bohr);
    bb.margin(1e-12);
    REQUIRE( q(0) == bb );
  }
  
}

// Interfragment bond
TEST_CASE("Connectivity for stretched H2"){
  using namespace std;
  
  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;
  
  double d{2.5};
  
  // Define compressed H2 molecule
  Molecule<vec3> molecule{
      {"H", {0.0, 0.0, d  }},
      {"H", {0.0, 0.0, 0.0}}
  };
  
  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);
  
  // Compute interatomic distance for compressed H2
  mat dd{ distances<vec3, mat>(molecule) };
  
  // Compute adjacency graph for compressed H2
  UGraph adj{ adjacency_matrix(dd, molecule) };
  
  // Compute distance and predecessor matrices for compressed H2
  Mat<int> dist, predecessors;
  std::tie(dist, predecessors) = distance_matrix<Mat<int>>(adj);
  
  SECTION("Bond"){
    // Compute bond
    std::vector<Bond> B{ bonds(dist, molecule) };
    
    // Check number of bonds
    REQUIRE( B.size() == 1);
    
    // Compute IRC
    vec q{irc_from_bad<vec3,vec>(to_cartesian<vec3,vec>(molecule), B, {}, {})};
    
    // Check number of IRC
    REQUIRE( linalg::size<vec>(q) == 1);
    
    // Check bond length
    Approx bb(d * angstrom_to_bohr);
    bb.margin(1e-12);
    REQUIRE( q(0) == bb );
  }
}

// Angle
TEST_CASE("Connectivity for compressed H2O"){
  using namespace std;
  
  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;
  
  double d1{0.8};
  double d2{0.6};
  double angle{134};
  
  double a{(180. - angle) / 2.};
  
  // Define compressed H2 molecule
  Molecule<vec3> molecule{
      {"O",{0, 0, 0}},
      {"H",{ d1 * std::cos(a * deg_to_rad), -d1 * std::sin(a * deg_to_rad), 0}},
      {"H",{-d2 * std::cos(a * deg_to_rad), -d2 * std::sin(a * deg_to_rad), 0}}
  };
  
  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);
  
  // Compute interatomic distance for compressed H2
  mat dd{ distances<vec3, mat>(molecule) };
  
  // Compute adjacency graph for compressed H2
  UGraph adj{ adjacency_matrix(dd, molecule) };
  
  // Compute distance and predecessor matrices for compressed H2
  Mat<int> dist, predecessors;
  std::tie(dist, predecessors) = distance_matrix<Mat<int>>(adj);
  
  // Compute bonds
  std::vector<Bond> B{ bonds(dist, molecule) };
  
  // Check number of bonds
  REQUIRE( B.size() == 2);
  
  // Compute angles
  std::vector<Angle> A{ angles(dist, predecessors, molecule) };
  
  // Check number of angles
  REQUIRE( A.size() == 1);
  
  // Compute IRC
  vec q{irc_from_bad<vec3,vec>(to_cartesian<vec3,vec>(molecule), B, A, {})};
  
  // Check number of IRC
  REQUIRE( linalg::size<vec>(q) == 3);
  
  SECTION("Bonds"){
    Approx bb1(d1 * angstrom_to_bohr);
    bb1.margin(1e-12);
    REQUIRE( q(0) == bb1 );
  
    Approx bb2(d2 * angstrom_to_bohr);
    bb2.margin(1e-12);
    REQUIRE( q(1) == bb2 );
  }
  
  SECTION("Angle"){
    Approx aa(angle);
    aa.margin(1e-12);
    REQUIRE( q(2) == aa );
  }
}

// Angle from interfragment bonds
TEST_CASE("Connectivity for stretched H2O"){
  using namespace std;
  
  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;
  
  double d1{1.3};
  double d2{1.4};
  double angle{102.03};
  
  double a{(180. - angle) * deg_to_rad};
  
  double cos_a{ std::cos(a)};
  double sin_a{ std::sin(a) };
  
  // Define compressed H2 molecule
  Molecule<vec3> molecule{
      {"O",{0., 0., 0.}},
      {"H",{-d1, 0., 0.}},
      {"H",{d2 * cos_a, -d2 * sin_a, 0.}}
  };
  
  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);
  
  // Compute interatomic distance for compressed H2
  mat dd{ distances<vec3, mat>(molecule) };
  
  // Compute adjacency graph for compressed H2
  UGraph adj{ adjacency_matrix(dd, molecule) };
  
  // Compute distance and predecessor matrices for compressed H2
  Mat<int> dist, predecessors;
  std::tie(dist, predecessors) = distance_matrix<Mat<int>>(adj);
  
  // Compute bonds
  std::vector<Bond> B{ bonds(dist, molecule) };
  
  // Check number of bonds
  REQUIRE( B.size() == 3); // Three interfragment bonds
  
  // Compute angles
  std::vector<Angle> A{ angles(dist, predecessors, molecule) };
  
  // Check number of angles
  REQUIRE( A.size() == 0);  // No angles for this bonding structure
  
  // Compute IRC
  vec q{irc_from_bad<vec3,vec>(to_cartesian<vec3,vec>(molecule), B, A, {})};
  
  // Because the three atoms belong to three different fragments, there
  // are three bonds (and no angles) for this structure.
  // TODO: Check this properly with other codes!!!
  
  // Check number of IRC
  REQUIRE( linalg::size<vec>(q) == 3);
  
  SECTION("Bonds"){
    // O-H1
    Approx bb1(d1 * angstrom_to_bohr);
    bb1.margin(1e-12);
    REQUIRE( q(0) == bb1 );
    
    // O-H2
    Approx bb2(d2 * angstrom_to_bohr);
    bb2.margin(1e-12);
    REQUIRE( q(1) == bb2 );
  
    // H1-H2
    Approx bb3( distance<vec3>(molecule[1].position, molecule[2].position) );
    bb3.margin(1e-12);
    REQUIRE( q(2) == bb3 );
  }
}

// Dihedral
TEST_CASE("Connectivity for compressed H2O2"){
  // TODO
}

// Dihedral from interfragment bonds
TEST_CASE("Connectivity for stretched H2O2"){
  // TODO
}

// Hydrogen bond (without quasi-linear angles)
TEST_CASE("Connectivity for bent water dimer"){
  using namespace std;
  
  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;
  
  double d1{1.3};
  double d2{1.4};
  double angle{102.03};
  
  double a{(180. - angle) * deg_to_rad};
  
  double cos_a{ std::cos(a)};
  double sin_a{ std::sin(a) };
  
  // Define compressed H2 molecule
  Molecule<vec3> molecule{
      {"O",{-1.464,   0.099,  -0.300}},
      {"H",{-1.956,   0.624,  -0.340}},
      {"H",{-1.797,  -0.799,   0.206}},
      {"O",{ 1.369,   0.146,  -0.395}},
      {"H",{ 1.894,   0.486,   0.335}},
      {"H",{ 0.451,   0.165,  -0.083}}
  };

  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);
  
  // Compute interatomic distance for compressed H2
  mat dd{ distances<vec3, mat>(molecule) };
  
  // Compute adjacency graph for compressed H2
  UGraph adj{ adjacency_matrix(dd, molecule) };
  
  // Compute distance and predecessor matrices for compressed H2
  Mat<int> dist, predecessors;
  std::tie(dist, predecessors) = distance_matrix<Mat<int>>(adj);
  
  // Compute bonds
  std::vector<Bond> B{ bonds(dist, molecule) };
  
  // Check number of bonds
  REQUIRE( B.size() == 5); // 4 regular bonds plus H-bond
  
  // Compute angles
  std::vector<Angle> A{ angles(dist, predecessors, molecule) };
  
  // Check number of angles
  REQUIRE( A.size() == 5);
  
  // Compute dihedral angles
  std::vector<Dihedral> D{ dihedrals(dist, predecessors, molecule) };
  
  // TODO: Check wit other codes (where dihedral 2-1-3-6 is added).
  
  // Check number of angles
  REQUIRE( D.size() == 3);
  
  // Compute IRC
  vec q{irc_from_bad<vec3,vec>(to_cartesian<vec3,vec>(molecule), B, A, D)};
  
  // Check number of IRC
  REQUIRE( linalg::size<vec>(q) == 13);
  
  SECTION("Bonds"){
    // TODO
  }
  
  SECTION("Angles"){
    // TODO
  }
  
  SECTION("Dihedrals"){
    // TODO
  }
}

// Quasi-linear angles
TEST_CASE("Connectivity for water dimer"){
  // TODO
}

// Interfragment bonds
TEST_CASE("Connectivity for benzene dimer"){

}

TEST_CASE("Connectivity test for CH2O"){
  using namespace std;
  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;

  // Define formaldehyde molecule (CH2O)
  Molecule<vec3> molecule{
      {"C", {0.000000,  0.000000,  -0.537500}},
      {"O", {0.000000,  0.000000,   0.662500}},
      {"H", {0.000000,  0.866025,  -1.037500}},
      {"H", {0.000000, -0.866025,  -1.037500}}
  };
  
  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);
  
  // Compute interatomic distance for formaldehyde molecule
  mat dd{ distances<vec3, mat>(molecule) };
  
  // Print interatomic distances for formaldehyde molecule
  cout << "Distances" << endl;
  cout << dd * bohr_to_angstrom << endl;
  
  UGraph adj{ adjacency_matrix(dd, molecule) };
  
  Mat<int> dist, predecessors;
  std::tie(dist, predecessors) = distance_matrix<Mat<int>>(adj);
  
  // Distance matrix
  cout << "Distance matrix:" << endl;
  cout << dist << endl;
  
  // Predecessors matrix
  cout << "Predecessors matrix:" << endl;
  cout << predecessors << endl;
  
  SECTION("Bonds"){
    // Compute bonds
    std::vector<Bond> B{ bonds(dist, molecule) };
  
    // Check number of bonds
    REQUIRE( B.size() == 3);
  
    // Define correct bond lengths
    std::vector<double> bb{ 1.2, 1., 1. };
    
    cout << "\nBonds:" << endl;
    for(size_t i{0}; i < B.size(); i++){
      cout << bond(B[i], molecule) * bohr_to_angstrom << endl;
      
      Approx target{ bb[i] };
      
      target.margin(1e-6);
      
      REQUIRE( bond(B[i], molecule) * bohr_to_angstrom == target );
    }
  }
  
  std::vector<Angle> A{angles(dist, predecessors, molecule)};
  cout << "\nAngles:" << endl;
  for(const auto& a : A){
    cout << angle(a, molecule) << endl;
  }
  
  REQUIRE( A.size() == 3);
}

TEST_CASE("Connectivity test with molecule from input"){
  using namespace std;
  
  using namespace io;
  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;
  
  // Load molecule from file
  Molecule<vec3> molecule{ load_xyz<vec3>(config::molecules_dir + "caffeine.xyz") };
  
  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);
  
  // Compute interatomic distance for formaldehyde molecule
  mat dd{ distances<vec3, mat>(molecule) };

  UGraph adj{ adjacency_matrix(dd, molecule) };
  
  Mat<int> dist, predecessors;
  std::tie(dist, predecessors) = distance_matrix<Mat<int>>(adj);
  
  std::cout << "\nPredecessor matrix:\n" << predecessors << std::endl;
  
  std::cout << "\nDistance matrix:\n" << dist << std::endl;
  
  // Compute bonds
  std::vector<Bond> B{ bonds(dist, molecule) };
  
  // Print bonds
  cout << '\n' << B.size() << " bonds:" << endl;
  for(const auto& b : B){
    cout << '(' << b.i + 1 << ',' << b.j + 1 << ") "
         << bond(b, molecule) * bohr_to_angstrom << endl;
  }
  
  std::vector<Angle> A{angles(dist, predecessors, molecule)};
  cout << '\n' << A.size() << " angles:" << endl;
  for(const auto& a : A){
    cout << '(' << a.i + 1 << ',' << a.j + 1 << ',' << a.k + 1 << ") "
         << angle(a, molecule) << endl;
  }
  
  std::vector<Dihedral> D{dihedrals(dist, predecessors, molecule)};
  cout << '\n' << D.size() << " dihedrals:" << endl;
  for(const auto& d : D){
    cout << '(' << d.i + 1 << ',' << d.j + 1 << ','
         << d.k + 1 << ',' << d.l + 1 << ") "
         << dihedral(d, molecule) << endl;
  }
}

TEST_CASE("Fragment recognition"){
  using namespace std;
  
  using namespace io;
  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;
  
  // Load molecule from file
  Molecule<vec3> molecule{ load_xyz<vec3>(config::molecules_dir + "benzene_dimer.xyz") };
  
  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);
  
  // Compute interatomic distance for formaldehyde molecule
  mat dd{ distances<vec3, mat>(molecule) };
  
  std::cout << "\nDistance matrix:\n" << dd << std::endl;
  
  try{
    UGraph adj{ adjacency_matrix(dd, molecule) };
  }
  catch(const std::logic_error& le){
    cerr << le.what() << endl;
  }
  
}