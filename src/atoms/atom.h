#ifndef IRC_ATOM_H
#define IRC_ATOM_H

#include "periodic_table.h"

#include <ostream>
#include <string>

#include <armadillo>

namespace atom{

struct AtomicNumber{
 public:
  AtomicNumber(size_t an);
  
  size_t atomic_number;
};

struct Atom {
  Atom(const AtomicNumber& an, const arma::vec& pos);
  
  AtomicNumber atomic_number;
  
  arma::vec position;
};

std::string symbol(const AtomicNumber& an);

double mass(const AtomicNumber& an);

double covalent_radius(const AtomicNumber& an);

std::ostream& operator<<(std::ostream& out, const Atom& a);

}

#endif //IRC_ATOM_H