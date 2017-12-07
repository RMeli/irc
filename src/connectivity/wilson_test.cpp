#define CATCH_CONFIG_MAIN
#include "../catch/catch.hpp"

#include "wilson.h"

#include "../atoms/atom.h"
#include "../atoms/molecule.h"
#include "connectivity.h"

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
  
  SECTION("H2 bond stretching"){
    // Define dihydrogen molecule (H2)
    molecule::Molecule<vec3> molecule{
        {"H", {0.00,  0.00,  0.0}},
        {"H", {1.00,  0.00,  0.0}}
    };
  
    
    double d{0.01};
    vec dx{-d, 0.00, 0.00, d, 0.00, 0.00};
  
    mat C{ connectivity::connectivity_matrix<vec3, mat>(molecule)};
  
    vector<connectivity::Bond<vec3>> bonds{ connectivity::bonds(molecule, C)};
  
    mat Bwilson = wilson_matrix<vec3, mat>(molecule.size(), bonds);
  
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
            {cos(angle_rad), -sin(angle_rad), 0},
            {sin(angle_rad),  cos(angle_rad), 0},
            {             0,               0, 1}
        },
        {
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
        },
        {
            {cos(-angle_rad), -sin(-angle_rad), 0},
            {sin(-angle_rad),  cos(-angle_rad), 0},
            {              0,                0, 1}
        },
    };
  
    
    molecule::Molecule<vec3> molecule{
        {"H", { 1.43,  -1.10,  0.00}},
        {"O", { 0.00,   0.00,  0.00}},
        {"H", {-1.43,  -1.10,  0.00}}
    };
    
    vec dx{ linalg::zeros<vec>(3 * molecule.size()) };
    for(size_t i{0}; i < 3; i++){
      vec3 v{ R[i] * molecule[i].position - molecule[i].position};
      
      dx(3*i + 0) = v(0);
      dx(3*i + 1) = v(1);
      dx(3*i + 2) = v(2);
    }
    
    /*
    molecule::Molecule<vec3> molecule{
        {"H", {0.00,   0.00, -0.96}},
        {"O", {0.00,   0.00,  0.00}},
        {"H", {0.00,   0.93,  0.24}},
    };
     */
    
    mat C{ connectivity::connectivity_matrix<vec3, mat>(molecule)};
    
    cout << "Connectivity matrix:" << endl;
    cout << C << endl;
  
    vector<connectivity::Bond<vec3>> bonds{ connectivity::bonds(molecule, C)};
    cout << "\nBonds:" << endl;
    for(const auto& b : bonds){
      cout << b.bond << endl;
    }
  
    vector<connectivity::Angle<vec3>> angles{ connectivity::angles(molecule, C)};
    cout << "\nAngles:" << endl;
    for(const auto& a : angles){
      cout << a.angle << endl;
    }
  
    mat Bwilson = wilson_matrix<vec3, mat>(molecule.size(), bonds, angles);
  
    cout << "\nWilson B matrix:" << endl;
    cout << Bwilson << endl;
  
    vec displacement{Bwilson * dx};
    cout << "\nDisplacement:" << endl;
    cout << displacement << endl;
    
    SECTION("Bond change"){
      Approx target{0};
      
      target.margin(1e-4);
      
      REQUIRE( displacement(0) == target );
      REQUIRE( displacement(1) == target );
    }
    SECTION("Angle change"){
      Approx target{2 * angle};
    
      target.margin(1e-3);
    
      REQUIRE( displacement(2) * 180 / tools::constants::pi == target );
    }
  }
}