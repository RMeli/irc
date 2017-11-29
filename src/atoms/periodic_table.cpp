#include "periodic_table.h"

#include <iomanip>

namespace atoms{

size_t atomic_number(const std::string& symbol){
  size_t idx{0};
  for(size_t i{0}; i < pt_size; i++){
    if(symbol == pt_symbols[i]){
      idx = i;
      break;
    }
  }
  
  return idx;
}

void print_atom_info(const std::string& symbol, std::ostream& out){
  
  size_t idx{ atomic_number(symbol) };
  
  // Print top line
  out << std::left << std::setw(15) << std::setfill('-') << '+';
  out  << ' ';
  out << std::left<< std::setw(2) << std::setfill(' ') << pt_symbols[idx];
  out << ' ' ;
  out << std::right << std::setw(15) << std::setfill('-') << '+' << std::endl;
  
  // Print atomic number
  out << std::setfill(' ');
  out << std::left << std::setw(20) << "| atomic number:";
  out << std::right << std::setw(12) << idx << " |" << std::endl;
  
  // Print atomic mass
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "| atomic mass:";
  out << std::right << std::setw(12) << pt_masses[idx] << " |" << std::endl;
  
  // Print covalent radius
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "| covalent radius:";
  out << std::right << std::setw(12) << pt_covalent_radii[idx];
  out << " |" << std::endl;
  
  // Print bottom line
  out << std::left << std::setw(17) << std::setfill('-') << '+';
  out << std::right << std::setw(17) << std::setfill('-') << '+' << std::endl;
}

}