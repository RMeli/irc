#include "../../include/catch/catch.hpp"

#include "../../include/io.h"

#include "../../include/molecule.h"
#include "config.h"

#include <iostream>

#ifdef HAVE_ARMA
#include <armadillo>
using vec3 = arma::vec3;
#else
#error
#endif

using namespace irc;

TEST_CASE("Loading XYZ file"){
  
  using namespace io;
  using namespace molecule;
  
  // TODO: Filename independent of building directory
  Molecule<vec3> mol{ load_xyz<vec3>(config::molecules_dir + "caffeine.xyz") };
  
  std::cout << mol << std::endl;
}