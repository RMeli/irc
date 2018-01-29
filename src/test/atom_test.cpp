#include "../../include/catch/catch.hpp"

#include "libirc/atom.h"

#include "libirc/periodic_table.h"

#include <iostream>
#include <stdexcept>

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
  
  SECTION("invalid atom"){
    
    bool exception{false};
    
    try{
      Atom<vec3> a{periodic_table::pt_size + 1};
    }
    catch(const std::logic_error& e){
      exception = true;
    }
    
    REQUIRE( exception == true );
  }
  
  SECTION("atom from atomic number"){
    for(size_t i{1}; i < periodic_table::pt_size; i++) {
      REQUIRE( periodic_table::valid_atomic_number(i) );
  
      // Define atom from atomic number
      Atom<vec3> a{i};
      
      REQUIRE( symbol(a.atomic_number) == periodic_table::pt_symbols[i] );
      
      {
        Approx target{periodic_table::pt_masses[i]};
        
        target.margin(1e-12);
        
        REQUIRE( mass(a.atomic_number) == target );
      }
      
      {
        Approx target{periodic_table::pt_covalent_radii[i]};
  
        target.margin(1e-12);
        
        REQUIRE( covalent_radius(a.atomic_number) == target );
      }

      {
        Approx target{periodic_table::pt_vdv_radii[i]};

        target.margin(1e-12);

        REQUIRE( vdw_radius(a.atomic_number) == target );
      }

      {
        if(i == 1){
          REQUIRE( is_H(a.atomic_number) );
        }
        else if(i == 7 or i == 8 or i == 9 or i == 15 or i == 16 or i == 17){
          REQUIRE( is_NOFPSCl(a.atomic_number) );
        }
      }

    }
  }
  
  SECTION("atom from atomic symbol"){
    for(size_t i{1}; i < periodic_table::pt_size; i++) {
      
      REQUIRE( periodic_table::valid_atomic_number(i) );
    
      // Define atom from atomic number
      Atom<vec3> a{periodic_table::pt_symbols[i]};
    
      REQUIRE( symbol(a.atomic_number) == periodic_table::pt_symbols[i] );
      
      // Mass from atomic number
      {
        Approx target{periodic_table::pt_masses[i]};
  
        target.margin(1e-12);
        
        REQUIRE( mass(a.atomic_number) == target );
      }
    
      // Covalent radius from atomic number
      {
        Approx target{periodic_table::pt_covalent_radii[i]};
  
        target.margin(1e-12);
        
        REQUIRE( covalent_radius(a.atomic_number) == target );
      }

      // Van der Waals radius from atomic number
      {
        Approx target{periodic_table::pt_vdv_radii[i]};

        target.margin(1e-12);

        REQUIRE( vdw_radius(a.atomic_number) == target );
      }

      // Atom in H-bond
      {
        if(i == 1){
          REQUIRE( is_H(a.atomic_number) );
        }
        else if(i == 7 or i == 8 or i == 9 or i == 15 or i == 16 or i == 17){
          REQUIRE( is_NOFPSCl(a.atomic_number) );
        }
      }
    }
  }
  
}