#include "atom.h"

#include <iomanip>
#include <ostream>
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

Atom::Atom(const AtomicNumber& an, const arma::vec3& pos)
: atomic_number(an), position(pos)
{}

std::ostream& operator<<(std::ostream& out, const AtomicNumber& an){
  out << an.atomic_number;
  
  return out;
}

std::ostream& operator<<(std::ostream& out, const Atom& a){
  // Print top line
  out << std::left << std::setw(15) << std::setfill('-') << '+';
  out  << ' ';
  out << std::left<< std::setw(2) << std::setfill(' ');
  out << symbol(a.atomic_number) << ' ';
  out << std::right << std::setw(15) << std::setfill('-') << '+' << std::endl;
  
  // Print atomic number
  out << std::setfill(' ');
  out << std::left << std::setw(20) << "| atomic number:";
  out << std::right << std::setw(12) << a.atomic_number;
  out << " |" << std::endl;
  
  // Print atomic mass
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "| atomic mass:";
  out << std::right << std::setw(12);
  out << mass(a.atomic_number) << " |" << std::endl;
  
  // Print covalent radius
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "| covalent radius:";
  out << std::right << std::setw(12);
  out << covalent_radius(a.atomic_number) << " |" << std::endl;
  
  // Print position
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "| position:";
  out << std::right << std::setw(14);
  out << " |" << std::endl;
  
  // x coordinate
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "|    x =";
  out << std::right << std::setw(12) << a.position(0);
  out << " |" << std::endl;
  
  // y coordinate
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "|    y =";
  out << std::right << std::setw(12) << a.position(1);
  out << " |" << std::endl;
  
  // z coordinate
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "|    x =";
  out << std::right << std::setw(12) << a.position(2);
  out << " |" << std::endl;
  
  // Print bottom line
  out << std::left << std::setw(17) << std::setfill('-') << '+';
  out << std::right << std::setw(17) << std::setfill('-') << '+' << std::endl;
  
  return out;
}

}

