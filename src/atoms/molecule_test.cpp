#define CATCH_CONFIG_MAIN
#include "../catch/catch.hpp"

#include "molecule.h"

#include "../atoms/periodic_table.h"

#include <cassert>
#include <iostream>

#ifdef HAVE_ARMA
#include <armadillo>
using arma::vec3;
#elif HAVE_EIGEN3
# error // TODO: Solve problem with std::initializer_list
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