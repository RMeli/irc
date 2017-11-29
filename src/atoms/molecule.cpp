#include "molecule.h"

#include "atom.h"

#include <algorithm>

namespace molecule{

double mass(const Molecule& molecule){
  double m{0};
  for(auto atom : molecule){
    m += atom::mass(atom.atomic_number);
  }
  
  return m;
}

}