#include "atom.h"

#include <iomanip>

namespace atoms {

Atom::Atom(const std::string &s, const arma::vec &pos)
    : symbol(s),
      position(pos) ,
      number(atomic_number(symbol)),
      mass(pt_masses[number]),
      covalent_radius(pt_covalent_radii[number])
{}

Atom::Atom(size_t idx, const arma::vec &pos)
    : symbol(pt_symbols[idx]),
      number(idx),
      mass(pt_masses[idx]),
      covalent_radius(pt_covalent_radii[idx])
{}

const std::string& Atom::get_symbol() const{
  return symbol;
}

size_t Atom::get_atomic_number() const{
  return number;
}

double Atom::get_mass() const{
  return mass;
}

double Atom::get_covalent_radius() const{
  return covalent_radius;
}

const arma::vec& Atom::get_position() const{
  return position;
}

std::ostream& operator<<(std::ostream& out, const Atom& a){
  // Print top line
  out << std::left << std::setw(15) << std::setfill('-') << '+';
  out  << ' ';
  out << std::left<< std::setw(2) << std::setfill(' ') << a.get_symbol();
  out << ' ' ;
  out << std::right << std::setw(15) << std::setfill('-') << '+' << std::endl;
  
  // Print atomic number
  out << std::setfill(' ');
  out << std::left << std::setw(20) << "| atomic number:";
  out << std::right << std::setw(12) << a.get_atomic_number();
  out << " |" << std::endl;
  
  // Print atomic mass
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "| atomic mass:";
  out << std::right << std::setw(12) << a.get_mass() << " |" << std::endl;
  
  // Print covalent radius
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "| covalent radius:";
  out << std::right << std::setw(12) << a.get_covalent_radius();
  out << " |" << std::endl;
  
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
  out << std::right << std::setw(12) << a.get_position()[0];
  out << " |" << std::endl;
  
  // y coordinate
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "|    y =";
  out << std::right << std::setw(12) << a.get_position()[1];
  out << " |" << std::endl;
  
  // z coordinate
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "|    x =";
  out << std::right << std::setw(12) << a.get_position()[2];
  out << " |" << std::endl;
  
  // Print bottom line
  out << std::left << std::setw(17) << std::setfill('-') << '+';
  out << std::right << std::setw(17) << std::setfill('-') << '+' << std::endl;
  
  return out;
}

}

