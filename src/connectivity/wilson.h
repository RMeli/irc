#ifndef IRC_WILSON_H
#define IRC_WILSON_H

#include "connectivity.h"

#include <tuple>
#include <utility>

namespace wilson{

template <typename Vector3>
std::pair<Vector3> bond_gradient(const connectivity::Bond<Vector3>& b){
  Vector v{ (b.a1.position - b.a2.position) / b.bond };
  
  return {v, -v};
}

template <typename Vector3, typename Matrix>
Matrix wilson_matrix(cosnt std::vector<Bonds<Vector3>>& bonds, size_t n_atoms){
  Matrix B{ linalg::zeros(3 * n_atoms, bonds.size()) };
  
  size_t idx;
  Vector3 g1, g2;
  for(size_t i{0}; i < bonds.size(); i++){
    std::tie(g1, g2) = bond_gradient(bonds[i]);
  }
}

}

#endif //IRC_WILSON_H
