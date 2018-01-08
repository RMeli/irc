#define CATCH_CONFIG_MAIN
#include "catch/catch.hpp"

#include "irc.h"

#include "atoms/molecule.h"
#include "tools/conversion.h"

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

TEST_CASE("Internal Redundant Coordinates"){
  using namespace std;
  using namespace tools::conversion;
  using namespace molecule;
  
  SECTION("Initial hessian"){
  
    // Define formaldehyde molecule (CH2O)
    Molecule<vec3> molecule{
        {"C", {0.000000,  0.000000,  -0.537500}},
        {"O", {0.000000,  0.000000,   0.662500}},
        {"H", {0.000000,  0.866025,  -1.037500}},
        {"H", {0.000000, -0.866025,  -1.037500}}
    };
  
    // Transform molecular coordinates from angstrom to bohr
    multiply_positions(molecule, angstrom_to_bohr);
    
    // Build internal reaction coordinates
    IRC<vec3, vec, mat> irc(molecule);
    
    // Compute initial hessian
    mat H0{ irc.projected_initial_hessian_inv() };
    
    // Print initial hessian
    cout << "H0 =\n" << H0 << endl;
  }
}