#define CATCH_CONFIG_MAIN
#include "../catch/catch.hpp"

#include "wilson.h"

#include "../atoms/atom.h"
#include "../atoms/molecule.h"
#include "connectivity.h"


#include <iostream>

#ifdef HAVE_ARMA
#include <armadillo>
using vec3 = arma::vec3;
using vec = arma::vec;
using mat = arma::mat;
#else
#error
#endif

TEST_CASE("Wildon B matrix","[wilson]"){
  
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
  
    std::vector<connectivity::Bond<vec3>> bonds{ connectivity::bonds(molecule, C)};
  
    mat Bwilson = wilson_matrix<vec3, mat>(molecule.size(), bonds);
  
    cout << "Wilson B matrix:" << endl;
    cout << Bwilson << endl;
  
    cout << "\nTransformation" << endl;
    cout << Bwilson * dx << endl;
  }
  
  SECTION("H2O bending"){
  
  }
}