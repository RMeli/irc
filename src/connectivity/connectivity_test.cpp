#define CATCH_CONFIG_MAIN
#include "../catch/catch.hpp"

#include "connectivity.h"

#include "../atoms/atom.h"
#include "../atoms/molecule.h"
#include "../io/io.h"
#include "../tools/conversion.h"


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

TEST_CASE("Distance, angle and dihedral"){
  using namespace connectivity;
  
  // Define three points in 3D space
  vec3 p1{ 0.00,  0.00, -0.25};
  vec3 p2{ 0.00,  0.00,  1.50};
  vec3 p3{ 0.00,  1.00,  1.50};
  
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
  
  mat dist, predecessors;
  std::tie(dist, predecessors) = distance_matrix<mat>(adj) ;
  
  // Distance matrix
  cout << "Distance matrix:" << endl;
  cout << dist << endl;
  
  // Predecessors matrix
  cout << "Predecessors matrix:" << endl;
  cout << predecessors << endl;
  
  SECTION("Bonds"){
    // Compute bonds
    std::vector<Bond<vec3>> B{ bonds(dist, molecule) };
  
    // Check number of bonds
    REQUIRE( B.size() == 3);
  
    // Define correct bond lengths
    std::vector<double> bb{ 1.2, 1., 1. };
    
    cout << "\nBonds:" << endl;
    for(size_t i{0}; i < B.size(); i++){
      cout << B[i].bond * bohr_to_angstrom << endl;
      
      Approx target{ bb[i] };
      
      target.margin(1e-6);
      
      REQUIRE( B[i].bond * bohr_to_angstrom == target );
    }
  }
  
  std::vector<Angle<vec3>> A{angles(dist, predecessors, molecule)};
  cout << "\nAngles:" << endl;
  for(const auto& a : A){
    cout << a.angle << endl;
  }
  
  REQUIRE( A.size() == 3);
}

TEST_CASE("Connectivity test for caffeine"){
  using namespace std;
  using namespace io;
  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;
  
  // Load molecule from file
  Molecule<vec3> molecule{ load_xyz<vec3>("../test/caffeine.xyz") };
  
  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);
  
  // Compute interatomic distance for formaldehyde molecule
  mat dd{ distances<vec3, mat>(molecule) };

  UGraph adj{ adjacency_matrix(dd, molecule) };
  
  mat dist, predecessors;
  std::tie(dist, predecessors) = distance_matrix<mat>(adj);
  
  // Compute bonds
  std::vector<Bond<vec3>> B{ bonds(dist, molecule) };
  
  // Print bonds
  cout << '\n' << B.size() << " bonds:" << endl;
  for(const auto& b : B){
    cout << '(' << b.i + 1 << ',' << b.j + 1 << ") "
         << b.bond * bohr_to_angstrom << endl;
  }
  
  std::vector<Angle<vec3>> A{angles(dist, predecessors, molecule)};
  cout << '\n' << A.size() << " angles:" << endl;
  for(const auto& a : A){
    cout << '(' << a.i + 1 << ',' << a.j + 1 << ',' << a.k + 1 << ") "
         << a.angle << endl;
  }
}