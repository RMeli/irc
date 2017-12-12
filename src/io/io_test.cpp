#define CATCH_CONFIG_MAIN
#include "../catch/catch.hpp"

#include "io.h"

#include "../atoms/molecule.h"

#include <iostream>

#ifdef HAVE_ARMA
#include <armadillo>
using vec3 = arma::vec3;
#else
#error
#endif

TEST_CASE("Loading XYZ file"){
  
  using namespace io;
  using namespace molecule;
  
  // TODO: Filename independent of building directory
  Molecule<vec3> mol{ load_xyz<vec3>("../test/caffeine.xyz") };
  
  std::cout << mol << std::endl;
}