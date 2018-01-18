#include "../../include/catch/catch.hpp"

#include "libirc/io.h"

#include "libirc/molecule.h"
#include "config.h"

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

TEST_CASE("Loading XYZ file"){
  
  using namespace io;
  using namespace molecule;
  
  Molecule<vec3> mol{ load_xyz<vec3>(config::molecules_dir + "caffeine.xyz") };
  
  std::cout << mol << std::endl;
}