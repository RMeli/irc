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

  Molecule molecule{
      {6, {0.000000,  0.000000,  -0.537500}},
      {8, {0.000000,  0.000000,   0.662500}},
      {1, {0.000000,  0.866025,  -1.037500}},
      {1, {0.000000, -0.866025,  -1.037500}}
  };
  
  multiply_positions(molecule, angstrom_to_bohr);
  
  arma::mat connectivity{ connectivity_matrix(molecule) * bohr_to_angstrom };
  
  cout << molecule << endl;
  
  cout << "Connectivity matrix:" << endl;
  cout << connectivity_matrix(molecule) * bohr_to_angstrom << endl;
  
  arma::mat C{
      {0.0, 1.2, 1.0, 1.0},
      {1.2, 0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0, 0.0},
      {1.0, 0.0, 0.0, 0.0}
  };
  
  assert( nearly_equal(connectivity, C, 1e-6) );
}