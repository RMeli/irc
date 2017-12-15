#define CATCH_CONFIG_MAIN
#include "../catch/catch.hpp"

#include "math.h"

TEST_CASE("Angle in the interval ]-pi,pi]"){
  
  using tools::constants::pi;
  using namespace tools::math;
  
  SECTION("Angle in ]-pi,pi]"){
    double angle{ pi / 2.};
    
    Approx target{ pi / 2. };
    
    target.margin(1e-6);
    
    REQUIRE( pirange(angle) == target );
  }
  
  SECTION("Angle in ]pi,2*pi]"){
    double angle{pi / 2. + pi};
    
    Approx target{ -pi / 2. };
    
    target.margin(1e-6);
    
    REQUIRE( pirange(angle) == target );
  }
  
  SECTION("Angle in ]-2*pi,-pi]"){
    double angle{ -pi / 2. - pi};
    
    Approx target{ pi / 2. };
    
    target.margin(1e-6);
    
    REQUIRE( pirange(angle) == target );
  }
  
  SECTION("Angle"){
    double angle{ pi / 2. + 6 * pi };
    
    Approx target{ pi / 2. };
    
    target.margin(1e-6);
    
    REQUIRE( pirange(angle) == target );
  }
  
  SECTION("Angle"){
    double angle{ pi / 2. - 5 * pi };
    
    Approx target{ -pi / 2. };
    
    target.margin(1e-6);
    
    REQUIRE( pirange(angle) == target );
  }
  
}