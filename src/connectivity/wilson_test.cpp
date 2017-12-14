#define CATCH_CONFIG_MAIN
#include "../catch/catch.hpp"

#include "wilson.h"

#include "../atoms/atom.h"
#include "../atoms/molecule.h"
#include "connectivity.h"
#include "../tools/conversion.h"

#include <cmath>
#include <iostream>

#ifdef HAVE_ARMA
#include <armadillo>
using vec3 = arma::vec3;
using vec = arma::vec;
using mat = arma::mat;
#else
#error
#endif

TEST_CASE("Wilson B matrix","[wilson]"){
  
  using namespace std;
  using namespace wilson;
  
  SECTION("H2 stretching"){
    // Define dihydrogen molecule (H2)
    molecule::Molecule<vec3> molecule{
        {"H", {0.00,  0.00,  0.0}},
        {"H", {1.00,  0.00,  0.0}}
    };
    
    double d{0.01};
    vec dx{-d, 0.00, 0.00, d, 0.00, 0.00};
  
    mat Bwilson = wilson_matrix<vec3, mat>(molecule);
  
    cout << "Wilson B matrix:" << endl;
    cout << Bwilson << endl;
  
    vec transformation{Bwilson * dx};
    cout << "\nTransformation" << endl;
    cout << Bwilson * dx << endl;
    
    SECTION("Bond change"){
      Approx target{2 * d};
      
      target.margin(1e-5);
      
      REQUIRE( transformation(0) == target );
    }
  }
  
  SECTION("H2O bending"){
    
    double angle( 0.5 );
    double angle_rad( angle / 180. * tools::constants::pi );
    
    std::vector<mat> R{
        {
            // Rotation for H1
            {cos(angle_rad), -sin(angle_rad), 0},
            {sin(angle_rad),  cos(angle_rad), 0},
            {             0,               0, 1}
        },
        {
            // Rotation for O
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
        },
        {
            // Rotation for H2
            {cos(-angle_rad), -sin(-angle_rad), 0},
            {sin(-angle_rad),  cos(-angle_rad), 0},
            {              0,                0, 1}
        },
    };
    
    molecule::Molecule<vec3> molecule{
        {"H", { 1.43,  -1.10,  0.00}}, // H1
        {"O", { 0.00,   0.00,  0.00}}, // O
        {"H", {-1.43,  -1.10,  0.00}}  // H2
    };
    
    // Allocate displacements in cartesian coordinates
    vec dx{ linalg::zeros<vec>(3 * molecule.size()) };
    
    // Compute displacements
    for(size_t i{0}; i < 3; i++){
      vec3 v{ R[i] * molecule[i].position - molecule[i].position};
      
      dx(3*i + 0) = v(0);
      dx(3*i + 1) = v(1);
      dx(3*i + 2) = v(2);
    }
  
    // Compute interatomic distances for water
    mat dd{ connectivity::distances<vec3, mat>(molecule) };
  
    // Compute adjacency matrix (graph)
    connectivity::UGraph adj{ connectivity::adjacency_matrix(dd, molecule) };
  
    // Compute distance matrix and predecessor matrix
    mat dist, predecessors;
    std::tie(dist, predecessors) = connectivity::distance_matrix<mat>(adj) ;
  
    // Compute bonds
    vector<connectivity::Bond<vec3>> bonds{ connectivity::bonds(dist, molecule)};
    
    // Print bonds
    cout << "\nBonds:" << endl;
    for(const auto& b : bonds){
      cout << b.bond << endl;
    }
    
    // Check number of bonds
    REQUIRE( bonds.size() == 2 );
  
    // Compute angles
    vector<connectivity::Angle<vec3>> angles{
        connectivity::angles(dist, predecessors, molecule)};
    
    // Print angles
    cout << "\nAngles:" << endl;
    for(const auto& a : angles){
      cout << a.angle << endl;
    }
    
    // Check number of angles
    REQUIRE( angles.size() == 1 );
  
    // Compute Wilson's B matrix
    mat Bwilson = wilson_matrix<vec3, mat>(molecule.size(), bonds, angles);
  
    // Print Wilson B matrix
    cout << "\nWilson B matrix:" << endl;
    cout << Bwilson << endl;
  
    // Compute displacement in internal coordinates
    vec displacement{Bwilson * dx};
    cout << "\nDisplacement:" << endl;
    cout << displacement << endl;
    
    // Check that bonds do not change
    SECTION("Bond change"){
      Approx target{0};
      
      target.margin(1e-4);
      
      REQUIRE( displacement(0) == target );
      REQUIRE( displacement(1) == target );
    }
    
    // Check change in angle
    SECTION("Angle change"){
      Approx target{2 * angle};
    
      target.margin(1e-3);
    
      REQUIRE( displacement(2) * 180 / tools::constants::pi == target );
    }
  }
  
  SECTION("H2O2 torsion"){
  
    double angle( 1.0 );
    double angle_rad( angle / 180. * tools::constants::pi );
  
    mat R{
        {cos(angle_rad), -sin(angle_rad), 0},
        {sin(angle_rad), cos(angle_rad),  0},
        {0,              0,               1}
    };
    
    molecule::Molecule<vec3> molecule{
        {"H", {  0.000,  0.947, -0.079}}, // H1
        {"O", {  0.000,  0.000,  0.000}}, // O1
        {"O", {  0.000,  0.000,  1.474}}, // O2
        {"H", { -0.854, -0.407,  1.553}}  // H2
    };

    // Transform molecular coordinates from angstrom to bohr
    molecule::multiply_positions(molecule,
                                 tools::conversion::angstrom_to_bohr);
  
    // Allocate displacements in cartesian coordinates
    vec dx{ linalg::zeros<vec>(3 * molecule.size()) };
    
    // Apply rotation to H1
    vec3 v{ R * molecule[0].position - molecule[0].position};
    dx(0) = v(0);
    dx(1) = v(1);
    dx(2) = v(2);
    
    // Compute dihedral variation
    double d_old{ connectivity::dihedral<vec3>( molecule[0].position,
                                                molecule[1].position,
                                                molecule[2].position,
                                                molecule[3].position) };
    double d_new{ connectivity::dihedral<vec3>( R * molecule[0].position,
                                                molecule[1].position,
                                                molecule[2].position,
                                                molecule[3].position) };
    double d_diff{ d_new - d_old };
    
    // Compute interatomic distances
    mat dd{ connectivity::distances<vec3, mat>(molecule) };
  
    // Compute adjacency matrix (graph)
    connectivity::UGraph adj{ connectivity::adjacency_matrix(dd, molecule) };
  
    // Compute distance matrix and predecessor matrix
    mat dist, predecessors;
    std::tie(dist, predecessors) = connectivity::distance_matrix<mat>(adj) ;
  
    // Compute bonds
    vector<connectivity::Bond<vec3>> bonds{ connectivity::bonds(dist, molecule)};
  
    // Print bonds
    cout << "\nBonds:" << endl;
    for(const auto& b : bonds){
      cout << b.bond * tools::conversion::bohr_to_angstrom << endl;
    }
  
    // Check number of bonds
    REQUIRE( bonds.size() == 3 );
  
    // Compute angles
    vector<connectivity::Angle<vec3>> angles{
        connectivity::angles(dist, predecessors, molecule)};
  
    // Print angles
    cout << "\nAngles:" << endl;
    for(const auto& a : angles){
      cout << a.angle << endl;
    }
  
    // Check number of angles
    REQUIRE( angles.size() == 2 );
  
    // Compute dihedrals
    vector<connectivity::Dihedral<vec3>> dihedrals{
        connectivity::dihedrals(dist, predecessors, molecule)};
  
    // Print angles
    cout << "\nDihedrals:" << endl;
    for(const auto& d : dihedrals){
      cout << d.dihedral << endl;
    }
  
    // Check number of dihedrals
    REQUIRE( dihedrals.size() == 1 );
  
    // Compute Wilson's B matrix
    mat Bwilson = wilson_matrix<vec3, mat>(molecule.size(),
                                           bonds, angles, dihedrals);
  
    // Print Wilson B matrix
    cout << "\nWilson B matrix:" << endl;
    cout << Bwilson << endl;
    
    // Compute displacement in internal coordinates
    vec displacement{Bwilson * dx};
    cout << "\nDisplacement:" << endl;
    cout << displacement << endl;
  
    // Check that bonds do not change
    SECTION("Bond change"){
      Approx target{0};
    
      target.margin(1e-3);
    
      REQUIRE( displacement(0) == target );
      REQUIRE( displacement(1) == target );
      REQUIRE( displacement(2) == target );
    }
  
    // Check change in angle
    SECTION("Angle change"){
      Approx target{0};
    
      target.margin(1e-3);
  
      REQUIRE( displacement(3) == target );
      REQUIRE( displacement(4) == target );
    }
  
    // Check change in angle
    SECTION("Dihedral change"){
      Approx target{d_diff};
      
      target.margin(1e-4);
    
      REQUIRE( displacement(5) * 180 / tools::constants::pi == target );
    }
  }
}