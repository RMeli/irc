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
  Molecule<vec3> molecule{ load_xyz<vec3>(config::molecules_dir + "water_dimer_2.xyz") };
  
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