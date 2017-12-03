#include "atom.h"

#include <stdexcept>

namespace atom {

AtomicNumber::AtomicNumber(size_t an) {
  if( !periodic_table::valid_atomic_number(an) ){
    throw std::logic_error("");
  }
  
  atomic_number = an;
}

std::string symbol(const AtomicNumber& an){
  return periodic_table::pt_symbols[an.atomic_number];
}

double mass(const AtomicNumber& an){
  return periodic_table::pt_masses[an.atomic_number];
}

double covalent_radius(const AtomicNumber& an){
  return periodic_table::pt_covalent_radii[an.atomic_number];
}

std::ostream& operator<<(std::ostream& out, const AtomicNumber& an){
  out << an.atomic_number;
  
  return out;
}

}

