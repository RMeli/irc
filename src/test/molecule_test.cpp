#include "../../include/catch/catch.hpp"

#include "libirc/molecule.h"

#include "libirc/periodic_table.h"

#include <cassert>
#include <iostream>

#ifdef HAVE_ARMA
#include <armadillo>
using arma::vec3;
#elif HAVE_EIGEN3
#include <Eigen3/Eigen/Dense>
using vec3 = Eigen::Vector3d;
#else
#error
#endif

using namespace irc;

TEST_CASE("Molecule") {
  
  using namespace std;
  using namespace molecule;
  using namespace periodic_table;
  
  
  Molecule<vec3> molecule{{
                              {1, {0.0, 1.1, 2.2}},
                              {2, {0.0, 1.1, 2.2}},
                              {3, {0.0, 1.1, 2.2}}
                          }};
  
  SECTION("Mass") {
    Approx target{pt_masses[1] + pt_masses[2] + pt_masses[3]};
    
    target.margin(1e-12);
    
    REQUIRE( mass(molecule) == target );
  }
  
  SECTION("Positon multiplier") {
    multiply_positions(molecule, 2.);
    
    for (const auto &atom : molecule) {
      
      vec3 pos{{0.0, 2.2, 4.4}};
      
      for(unsigned i{0}; i < 3; i++){
        Approx target{pos(i)};
        
        target.margin(1e-12);
        
        REQUIRE( atom.position(i) == target );
      }
    }
  }
  
}