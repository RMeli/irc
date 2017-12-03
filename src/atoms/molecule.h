#ifndef IRC_MOLECULE_H
#define IRC_MOLECULE_H

#include "atom.h"

#include <ostream>
#include <vector>

namespace molecule{

template<typename T>
using Molecule = std::vector<atom::Atom<T>>;

template<typename T>
double mass(const Molecule<T>& molecule){
  double m{0};
  for(const auto& atom : molecule){
    m += atom::mass(atom.atomic_number);
  }
  
  return m;
}

template<typename T>
void multiply_positions(Molecule<T>& molecule, double multiplier){
  for(auto& atom : molecule){
    atom.position = atom.position * multiplier;
  }
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const Molecule<T>& molecule){
  for(const auto& atom : molecule){
    out << atom;
  }
  
  return out;
}

}


#endif //IRC_MOLECULE_H
