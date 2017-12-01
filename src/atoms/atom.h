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

// No ARMADILLO as dependency
//#ifdef ARMADILLO
//class Matrix : arma::mat {};
//#else EIGEN
//class Matrix : eigen::mat {};
//#endif

struct Atom {
  Atom(const AtomicNumber& an, const arma::vec3& pos = {0., 0. , 0.});
  Atom(const std::string& symbol, const arma::vec3& pos = {0., 0. , 0.});
  
  AtomicNumber atomic_number;
  
  arma::vec3 position;
};

std::string symbol(const AtomicNumber& an);

double mass(const AtomicNumber& an);

double covalent_radius(const AtomicNumber& an);

std::ostream& operator<<(std::ostream& out, const Atom& a);

}

#endif //IRC_ATOM_H