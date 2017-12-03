#undef NDEBUG

#include "molecule.h"

#include "../tools/comparison.h"

#include <cassert>
#include <iostream>

#ifdef HAVE_ARMA
#include <armadillo>
using vec3 = arma::vec3;
#else
#error
#endif

int main(){
  using namespace std;
  using namespace molecule;
  using namespace tools::comparison;
  
  Molecule<vec3> molecule{
      {1, {0.0, 1.1, 2.2}},
      {2, {0.0, 1.1, 2.2}},
      {3, {0.0, 1.1, 2.2}},
      {4, {0.0, 1.1, 2.2}}
  };
  
  cout << molecule << endl;
  
  double m{mass(molecule)};
  
  assert( nearly_equal(m, 1.0079 + 4.0026 + 6.941 + 9.0122) );
  
  multiply_positions(molecule, 2);
  
  for(const auto& atom : molecule){
    assert( nearly_equal(atom.position, vec3{0.0, 2.2, 4.4}) );
  }
  
  return 0;
}