#include "../../include/catch/catch.hpp"

#include "../../include/atom.h"

#include "../../include/periodic_table.h"

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

TEST_CASE("Test atom and periodic table lookup functions","[atom]"){
  
  using namespace atom;
  
  SECTION("Atom from atomic number"){
    for(size_t i{1}; i < periodic_table::pt_size; i++) {
      REQUIRE( periodic_table::valid_atomic_number(i) );
  
      // Define atom from atomic number
      Atom<vec3> a{i};
      
      SECTION("symbol from atomic number"){
        REQUIRE( symbol(a.atomic_number) == periodic_table::pt_symbols[i] );
      }
      
      SECTION("mass from atomic number"){
        Approx target{periodic_table::pt_masses[i]};
        
        target.margin(1e-12);
        
        REQUIRE( mass(a.atomic_number) == target );
      }
      
      SECTION("covalent radius from atomic number"){
        Approx target{periodic_table::pt_covalent_radii[i]};
  
        target.margin(1e-12);
        
        REQUIRE( covalent_radius(a.atomic_number) == target );
      }
    }
  }
  
  SECTION("Atom from atomic symbol"){
    for(size_t i{1}; i < periodic_table::pt_size; i++) {
      REQUIRE( periodic_table::valid_atomic_number(i) );
    
      // Define atom from atomic number
      Atom<vec3> a{periodic_table::pt_symbols[i]};
    
      SECTION("symbol from atomic number"){
        REQUIRE( symbol(a.atomic_number) == periodic_table::pt_symbols[i] );
      }
    
      SECTION("mass from atomic number"){
        Approx target{periodic_table::pt_masses[i]};
  
        target.margin(1e-12);
        
        REQUIRE( mass(a.atomic_number) == target );
      }
    
      SECTION("covalent radius from atomic number"){
        Approx target{periodic_table::pt_covalent_radii[i]};
  
        target.margin(1e-12);
        
        REQUIRE( covalent_radius(a.atomic_number) == target );
      }
    }
  }
  
}