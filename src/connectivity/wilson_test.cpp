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
  
    vec dx{0.01, 0.00, 0.00, -0.01, 0.00, 0.00};
  
    mat C{ connectivity::connectivity_matrix<vec3, mat>(molecule)};
  
    vector<connectivity::Bond<vec3>> bonds{ connectivity::bonds(molecule, C)};
  
    mat Bwilson = wilson_matrix<vec3, mat>(molecule.size(), bonds);
  
    cout << "Wilson B matrix:" << endl;
    cout << Bwilson << endl;
  
    cout << "\nTransformation" << endl;
    cout << Bwilson * dx << endl;
  }
  
  SECTION("H2O bending"){
    double angle( 1 / 180. * tools::constants::pi );
    
    std::vector<mat> R{
        {
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1}
        },
        {
            {cos(angle), -sin(angle), 0},
            {sin(angle),  cos(angle), 0},
            {         0,           0, 1}
        },
        {
            {cos(-angle), -sin(-angle), 0},
            {sin(-angle),  cos(-angle), 0},
            {          0,            0, 1}
        },
    };
  
    molecule::Molecule<vec3> molecule{
        {"O", { 0.00,   0.00,  0.00}},
        {"H", { 1.43,  -1.10,  0.00}},
        {"H", {-1.43,  -1.10,  0.00}}
    };
    
    vec dx{ linalg::zeros<vec>(3 * molecule.size()) };
    for(size_t i{0}; i < 3; i++){
      vec3 v{ R[i] * molecule[i].position - molecule[i].position};
      
      dx(3*i + 0) = v(0);
      dx(3*i + 1) = v(1);
      dx(3*i + 2) = v(2);
    }
  
    mat C{ connectivity::connectivity_matrix<vec3, mat>(molecule)};
  
    vector<connectivity::Bond<vec3>> bonds{ connectivity::bonds(molecule, C)};
    vector<connectivity::Angle<vec3>> angles{ connectivity::angles(molecule, C)};
  
    mat Bwilson = wilson_matrix<vec3, mat>(molecule.size(), bonds, angles);
  
    cout << "Wilson B matrix:" << endl;
    cout << Bwilson << endl;
  
    cout << "\nTransformation" << endl;
    cout << Bwilson * dx << endl;
  }
}