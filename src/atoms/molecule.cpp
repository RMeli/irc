#include "molecule.h"

#include "atom.h"

#include <algorithm>

namespace molecule{

double mass(const Molecule& molecule){
  double m{0};
  for(const auto& atom : molecule){
    m += atom::mass(atom.atomic_number);
  }
  
  return m;
}

void multiply_positions(Molecule& molecule, double multiplier){
  for(auto& atom : molecule){
    atom.position *= multiplier;
  }
}

std::ostream& operator<<(std::ostream& out, const Molecule& molecule){
  for(const auto& atom : molecule){
    out << atom;
  }
  
  return out;
}

}