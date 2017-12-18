#define CATCH_CONFIG_MAIN
#include "../catch/catch.hpp"

#include "transformation.h"

#include "../atoms/molecule.h"
#include "../io/io.h"
#include "../tools/conversion.h"
#include "../connectivity/wilson.h"

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

TEST_CASE("Transformation"){
  SECTION("Cartesian to internal for ethanol"){
    using namespace molecule;
    using namespace connectivity;
    using namespace transformation;
    using namespace tools::conversion;
    using namespace io;
    using namespace std;
    
    // Load molecule from file
    Molecule<vec3> molecule{ load_xyz<vec3>("../test/ethanol.xyz") };
    
    // Transform molecular coordinates from angstrom to bohr
    multiply_positions(molecule, angstrom_to_bohr);
    
    // Compute interatomic distance for formaldehyde molecule
    mat dd{ distances<vec3, mat>(molecule) };
    
    // Build graph based on the adjacency matrix
    UGraph adj{ adjacency_matrix(dd, molecule) };
    
    // Compute distance matrix and predecessor matrix
    mat dist, predecessors;
    std::tie(dist, predecessors) = distance_matrix<mat>(adj);
    
    // Compute bonds
    std::vector<Bond<vec3>> B{ bonds(dist, molecule) };
    
    // Print bonds
    cout << '\n' << B.size() << " bonds (a.u.):" << endl;
    for(const auto& b : B){
      cout << '(' << b.i + 1 << ',' << b.j + 1 << ") "
           << b.bond << endl;
    }
    
    // Compute angles
    std::vector<Angle<vec3>> A{angles(dist, predecessors, molecule)};
    
    // Print angles
    cout << '\n' << A.size() << " angles (deg):" << endl;
    for(const auto& a : A){
      cout << '(' << a.i + 1 << ',' << a.j + 1 << ',' << a.k + 1 << ") "
           << a.angle << endl;
    }
    
    // Compute dihedral angles
    std::vector<Dihedral<vec3>> D{dihedrals(dist, predecessors, molecule)};
    
    // Print dihedral angles
    cout << '\n' << D.size() << " dihedrals (deg):" << endl;
    for(const auto& d : D){
      cout << '(' << d.i + 1 << ',' << d.j + 1 << ','
           << d.k + 1 << ',' << d.l + 1 << ") "
           << d.dihedral << endl;
    }
    
    // Compute number of cartesian coordinates
    size_t n_c{ 3 * molecule.size() };
    
    // Allocate vector for cartesian positions
    vec x_c{ linalg::zeros<vec>(n_c) };
    
    // Fill vector with cartesian positions
    for(size_t i{0}; i < molecule.size(); i++){
      x_c(3 * i + 0) = molecule[i].position(0) ;
      x_c(3 * i + 1) = molecule[i].position(1) ;
      x_c(3 * i + 2) = molecule[i].position(2) ;
    }
    
    // Print cartesian coordinates
    cout << "\nCartesian coordinates (a.u.):\n " << x_c << endl;
    
    // Compute and print internal redundant coordinates
    cout << "Internal redundant coordinates (a.u.):\n"
         << cartesian_to_irc(x_c, B, A, D) << endl;
  }
  
  SECTION("Internal to cartesian for H2"){
    using namespace molecule;
    using namespace connectivity;
    using namespace transformation;
    using namespace tools::conversion;
    using namespace io;
    using namespace std;
    using namespace wilson;
  
    // Load molecule from file
    Molecule<vec3> molecule{{ {"H", {0.,0.,0.}}, {"H", {1.,0.,0.}}  }};
  
    // Compute interatomic distance for formaldehyde molecule
    mat dd{ distances<vec3, mat>(molecule) };
  
    // Build graph based on the adjacency matrix
    UGraph adj{ adjacency_matrix(dd, molecule) };
  
    // Compute distance matrix and predecessor matrix
    mat dist, predecessors;
    std::tie(dist, predecessors) = distance_matrix<mat>(adj);
  
    // Compute bonds
    std::vector<Bond<vec3>> B{ bonds(dist, molecule) };
  
    // Print bonds
    cout << '\n' << B.size() << " bonds (a.u.):" << endl;
    for(const auto& b : B){
      cout << '(' << b.i + 1 << ',' << b.j + 1 << ") "
           << b.bond << endl;
    }
  
    // Check number of bonds
    REQUIRE( B.size() == 1);
  
    // Wilson B matrix
    mat W = wilson_matrix<vec3,mat>(molecule.size(), B);
    
    mat G, iG;
    std::tie(G, iG) = G_matirces(W);
  
    // Allocate vector for internal reaction coordinates
    vec q_irc{ B[0].bond };
    
    // Displacement in internal coordinates
    vec dq_irc{ 0.1 };
  
    // Compute number of cartesian coordinates
    size_t n_c{ 3 * molecule.size() };
  
    // Allocate vector for cartesian positions
    vec x_c_old{ linalg::zeros<vec>(n_c) };
  
    // Fill vector with cartesian positions
    for(size_t i{0}; i < molecule.size(); i++){
      x_c_old(3 * i + 0) = molecule[i].position(0) ;
      x_c_old(3 * i + 1) = molecule[i].position(1) ;
      x_c_old(3 * i + 2) = molecule[i].position(2) ;
    }
    
    // Compute new cartesian coordinates
    vec x_c{ irc_to_cartesian(q_irc, dq_irc, x_c_old, B, {}, {}, W, iG) };
    
    // Print cartesian coordinates
    cout << "\nNew cartesian coordinates (a.u.):\n " << x_c << endl;
    
    vec3 p1{x_c(0), x_c(1), x_c(2)};
    vec3 p2{x_c(3), x_c(4), x_c(5)};
    
    Approx target{ q_irc(0) + dq_irc(0) };
    target.margin(1e-6);
    
    REQUIRE( distance(p1,p2) == target);
  }
}