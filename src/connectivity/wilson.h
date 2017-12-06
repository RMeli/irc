#ifndef IRC_WILSON_H
#define IRC_WILSON_H

#include "connectivity.h"

#include "../tools/constants.h"

#include <cmath>
#include <iostream>
#include <tuple>
#include <utility>

namespace wilson{

template <typename Vector3>
std::pair<Vector3,Vector3> bond_gradient(const connectivity::Bond<Vector3>& b){
  Vector3 v{ (b.p1 - b.p2) / b.bond };
  
  return {v, -v};
}

template <typename Vector3>
std::tuple<Vector3, Vector3, Vector3> angle_gradient(const connectivity::Angle<Vector3>& a){
  double angle_rad{ a.angle / 180. * tools::constants::pi};
  
  // TODO: Check pyberny for more robust implementations
  
  double sin_angle{ std::sin(angle_rad) };
  double cos_angle{ std::cos(angle_rad) };
  
  Vector3 b21{ a.p1 - a.p2 };
  Vector3 b23{ a.p3 - a.p2 };
  
  double bond21{ linalg::norm(b21) };
  double bond23{ linalg::norm(b23) };
  
  b21 = b21 / bond21;
  b23 = b23 / bond23;
  
  Vector3 v1{ (cos_angle * b21 - b23) / (sin_angle * bond21)};
  Vector3 v3{ (cos_angle * b23 - b21) / (sin_angle * bond23)};
  Vector3 v2{ -v1 -v3 };
  
  return {v1, v2, v3};
}

/// Function computing Wilson's \f$\mathbf{B}\f$ matrix from a set of internal
/// redundant coordinates.
///
/// \tparam Vector3 3D vector type
/// \tparam Matrix Matrix type
/// \param n_atoms Total number of atoms
/// \param bonds Collection of bonds
/// \param angles Collection of angles between bonded atoms
/// \return Wilson's B matrix
///
/// This function resurns Wilson's \f$\mathbf{B}\f$ matrix given a collection
/// of bonds, angles and dihedral angles.
///
/// Wilson's \f$\mathbf{B}\f$ matrix
/// \f[
///   B_{ij} = \frac{\partial q_i}{\partial x_j}
/// \f]
/// defines the transfromation from Cartesia displacements
/// \f$\delta\mathbf{x}\f$ to redundant internal displacements
/// \f$\delta\mathbf{q}\f$:
/// \f[
///   \delta\mathbf{q} = \mathbf{B} \delta\mathbf{x}
/// \f]
///
/// More details can be found in Peng et al., J. Comp. Chem. 17, 49-56, 1996.
template <typename Vector3, typename Matrix>
Matrix wilson_matrix(size_t n_atoms,
                     const std::vector<connectivity::Bond<Vector3>>& bonds,
                     const std::vector<connectivity::Angle<Vector3>>& angles = {}){
  // Get the total number of internal redundant coordinates
  size_t n_irc{ bonds.size() + angles.size()};
  
  // Allocate Wilson's B matrix
  Matrix B{ linalg::zeros<Matrix>(n_irc, 3 * n_atoms) };
  
  // Utility vector for gradients storage
  Vector3 g1, g2, g3;
  
  // Populate B matrix's rows corresponding to bonds
  connectivity::Bond<Vector3> bond;
  for(size_t i{0}; i < bonds.size(); i++){
    bond = bonds[i];
    
    std::tie(g1, g2) = bond_gradient(bond);
    
    for(size_t idx{0}; idx < 3; idx++){
      B(i, 3 * bond.i + idx) = g1(idx);
      B(i, 3 * bond.j + idx) = g2(idx);
    }
  }
  
  // Populate B matrix's rows corresponding to angles
  connectivity::Angle<Vector3> angle;
  for(size_t i{0}; i < angles.size(); i++){
    angle = angles[i];
    
    std::tie(g1, g2, g3) = angle_gradient(angle);
  
    for(size_t idx{0}; idx < 3; idx++){
      B(i + bonds.size(), 3 * angle.i + idx) = g1(idx);
      B(i + bonds.size(), 3 * angle.j + idx) = g2(idx);
      B(i + bonds.size(), 3 * angle.k + idx) = g3(idx);
    }
  }
  
  // TODO: Populate B matrix's rows corresponding to dihedrals
  
  return B;
}

}

#endif //IRC_WILSON_H