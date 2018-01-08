#include "catch/catch.hpp"

#include "periodic_table.h"

#include <iomanip>
#include <iostream>

TEST_CASE("Test periodic table","[periodic_table]"){
  using namespace std;
  
  using namespace irc;
  using namespace periodic_table;
  
  bool verbose{false};
  
  // Check periodic table size for supported elements
  SECTION("size"){
    REQUIRE( pt_size == 96 );
  }
  
  // Check valid atomic numbers
  SECTION("atomic numbers"){
    bool valid{false};
    for(size_t i{0}; i < pt_size + 1; i++){
      valid = valid_atomic_number(i);
    
      if(valid){
        if(verbose){
          cout << std::setw(3) << i;
          cout << " is a valid atomic number corresponding to element ";
          cout << pt_symbols[i] << ".\n";
        }
  
        REQUIRE( i > 0);
        REQUIRE( i < pt_size);
      }
      else{
        if(verbose){
          cout << std::setw(3) << i;
          cout << " is not a valid atomic number.\n";
        }
      
        REQUIRE( (i == 0 or i == pt_size) );
      }
    }
  }
  
  // Check conversion from symbol to atomic number
  SECTION("atomic number from symbol"){
    for(size_t i{1}; i < pt_size; i++){
      REQUIRE( atomic_number(pt_symbols[i]) == i );
    }
  }
  
}