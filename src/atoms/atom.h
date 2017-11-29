#ifndef IRC_ATOM_H
#define IRC_ATOM_H

#include "periodic_table.h"

#include <ostream>
#include <string>

#include <armadillo>

namespace atoms{

class Atom {
 public:
  Atom(const std::string& s, const arma::vec& pos);
  Atom(size_t idx, const arma::vec& pos);
  
  const std::string& get_symbol() const;
  
  size_t get_atomic_number() const;
  
  double get_mass() const;
  
  double get_covalent_radius() const;
  
  const arma::vec& get_position() const;
 
 private:
  /// Atomic symbol
  std::string symbol;
  
  /// Atomic number
  size_t number;
  
  /// Atomic mass
  double mass;
  
  /// Covalent radius
  double covalent_radius;
  
  /// Position
  arma::vec position;
  
  /// Listo of atoms bonded with this
  std::vector<std::shared_ptr<Atom>> neighbours;
};

std::ostream& operator<<(std::ostream& out, const Atom& a);

}

#endif //IRC_ATOM_H