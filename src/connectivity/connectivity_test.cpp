#undef NDEBUG

#include "connectivity.h"

#include "../atoms/atom.h"
#include "../atoms/molecule.h"
#include "../tools/comparison.h"
#include "../tools/conversion.h"


#include <iostream>
#include <cassert>

#ifdef HAVE_ARMA
#include <armadillo>
using vec3 = arma::vec3;
using vec = arma::vec;
using mat = arma::mat;
#else
#error
#endif

int main(){
  using namespace std;
  using namespace tools::comparison;
  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;
  
  // Define three points in 3D space
  vec3 p1{ 0.00,  0.00, -0.25};
  vec3 p2{ 0.00,  0.00,  1.50};
  vec3 p3{ 0.00,  1.00,  1.50};
  
  // Test distance function
  assert( nearly_equal(distance(p1, p2), 1.75) );
  
  // Test angle function
  assert( nearly_equal(connectivity::angle(p1, p2, p3), 90) );

  // Define formaldehyde molecule (CH2O)
  Molecule<vec3> molecule{
      {"C", {0.000000,  0.000000,  -0.537500}},
      {"O", {0.000000,  0.000000,   0.662500}},
      {"H", {0.000000,  0.866025,  -1.037500}},
      {"H", {0.000000, -0.866025,  -1.037500}}
  };
  
  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);
  
  // Compute connectivity matrix for formaldehyde molecule (in angstrom)
  mat connectivity{ connectivity_matrix<vec3, mat>(molecule) * bohr_to_angstrom };
  
  // Print formaldehyde molecule
  cout << molecule << endl;
  
  // Print connectivity matrix for formaldehyde molecule
  cout << "Connectivity matrix:" << endl;
  cout << connectivity << endl;
  
  // Define correct connectivity matrix for comparison
  mat C{
      {0.0, 1.2, 1.0, 1.0},
      {1.2, 0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0, 0.0}
  };
  
  // Check connectivity matrix
  for(size_t j{0}; j < molecule.size(); j++) {
    for (size_t i{0}; i < molecule.size(); i++) {
      assert( nearly_equal(connectivity(i,j), C(i,j), 1e-6) );
    }
  }
  
  
  // Print bonds
  vector<double> b{ bonds(connectivity) };
  cout << "Bonds:" << endl;
  for(const auto& bond : b){
    cout << bond << endl;
  }
  
  // Print angles
  vector<double> ang{ angles(molecule, connectivity) };
  for(const auto& angle : ang){
    cout << angle << endl;
  }
  
  return 0;
}