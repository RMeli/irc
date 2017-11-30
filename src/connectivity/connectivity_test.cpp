#undef NDEBUG

#include "connectivity.h"

#include "../atoms/atom.h"
#include "../atoms/molecule.h"
#include "../tools/comparison.h"
#include "../tools/conversion.h"


#include <iostream>
#include <cassert>

#include <armadillo>

int main(){
  using namespace std;
  using namespace tools::comparison;
  using namespace tools::conversion;
  using namespace molecule;
  using namespace connectivity;
  
  // Define three points in 3D space
  arma::vec p1{ 0.00,  0.00, -0.25};
  arma::vec p2{ 0.00,  0.00,  1.50};
  arma::vec p3{ 0.00,  1.00,  1.50};
  
  // Test distance function
  assert( nearly_equal(distance(p1, p2), 1.75) );
  
  // Test angle function
  assert( nearly_equal(angle(p1, p2, p3), 90) );

  // Define formaldehyde molecule (CH2O)
  Molecule molecule{
      {6, {0.000000,  0.000000,  -0.537500}},
      {8, {0.000000,  0.000000,   0.662500}},
      {1, {0.000000,  0.866025,  -1.037500}},
      {1, {0.000000, -0.866025,  -1.037500}}
  };
  
  // Transform molecular coordinates from angstrom to bohr
  multiply_positions(molecule, angstrom_to_bohr);
  
  // Compute connectivity matrix for formaldehyde molecule (in angstrom)
  arma::mat connectivity{ connectivity_matrix(molecule) * bohr_to_angstrom };
  
  // Print formaldehyde molecule
  cout << molecule << endl;
  
  // Print connectivity matrix for formaldehyde molecule
  cout << "Connectivity matrix:" << endl;
  cout << connectivity << endl;
  
  // Define correct connectivity matrix for comparison
  arma::mat C{
      {0.0, 1.2, 1.0, 1.0},
      {1.2, 0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0, 0.0}
  };
  
  // Check connectivity matrix
  assert( nearly_equal(connectivity, C, 1e-6) );
  
  // Print bonds
  vector<double> b{ bonds(connectivity) };
  cout << "Bonds:" << endl;
  for(const auto& bond : b){
    cout << bond << endl;
  }
  
  // Print angles
  vector<double> ang{ angles(molecule, connectivity) };
  cout << "Angles:" << endl;
  for(const auto& angle : ang){
    cout << angle << endl;
  }
}