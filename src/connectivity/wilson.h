#ifndef IRC_WILSON_H
#define IRC_WILSON_H

#include "connectivity.h"

#include "../atoms/molecule.h"
#include "../connectivity/connectivity.h"
#include "../tools/constants.h"

#include <cmath>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

namespace wilson{

template <typename Vector3>
std::pair<Vector3,Vector3> bond_gradient(const connectivity::Bond<Vector3>& b){
  Vector3 v{ (b.p1 - b.p2) / b.bond };
  
  return {v, -v};
}

template <typename Vector3>
std::tuple<Vector3, Vector3, Vector3> angle_gradient(
    const connectivity::Angle<Vector3>& a){
  double angle_rad{ a.angle / 180. * tools::constants::pi};
  
  // TODO: Check pyberny for more robust implementation
  // https://github.com/azag0/pyberny
  
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

template <typename Vector3>
std::tuple<Vector3, Vector3, Vector3, Vector3> dihedral_gradient(
    const connectivity::Dihedral<Vector3>& d){
  
  // TODO: Check pyberny for more robust implementation
  // https://github.com/azag0/pyberny
  
  double angle123{ connectivity::angle(d.p1, d.p2, d.p3)};
  double angle123_rad{ angle123 * tools::constants::pi / 180. };
  double sin_angle123{ std::sin(angle123_rad) };
  double cos_angle123{ std::cos(angle123_rad) };
  
  double angle234{ connectivity::angle(d.p2, d.p3, d.p4) };
  double angle234_rad{ angle234 * tools::constants::pi / 180. };
  double sin_angle234{ std::sin(angle234_rad) };
  double cos_angle234{ std::cos(angle234_rad) };
  
  Vector3 b12{ d.p2 - d.p1 };
  Vector3 b23{ d.p3 - d.p2 };
  Vector3 b34{ d.p4 - d.p3 };
  
  double bond12{ linalg::norm(b12) };
  double bond23{ linalg::norm(b23) };
  double bond34{ linalg::norm(b34) };
  
  b12 = b12 / bond12;
  b23 = b23 / bond23;
  b34 = b34 / bond34;
  
  Vector3 b32{ -b23 };
  Vector3 b43{ -b34 };
  
  Vector3 v1{-linalg::cross(b12, b23) / (bond12 * sin_angle123 * sin_angle123)};
  
  double vc1{0.}, vc2{0.};
  Vector3 vv1, vv2{0.,0.,0.};
  
  vc1 = (bond23 - bond12 * cos_angle123) / (bond12 * bond23 * sin_angle123);
  vc2 = cos_angle234 / (bond23 * sin_angle234);
  vv1 = linalg::cross(b12, b23) / sin_angle123;
  vv2 = linalg::cross(b43, b32) / sin_angle234;
  
  Vector3 v2{ vc1 * vv1 + vc2 * vv2 };
  
  vc1 = (bond23 - bond34 * cos_angle234) / (bond23 * bond34 * sin_angle234);
  vc2 = cos_angle123 / (bond23 * sin_angle123);
  vv1 = linalg::cross(b43, b32) / sin_angle234;
  vv2 = linalg::cross(b12, b23) / sin_angle123;
  
  Vector3 v3{ vc1 * vv1 + vc2 * vv2 };

  Vector3 v4{-linalg::cross(b43, b32) / (bond34 * sin_angle234 * sin_angle234)};
  
  return {v1, v2, v3, v4};
}

/// Function computing Wilson's \f$\mathbf{B}\f$ matrix from a set of internal
/// redundant coordinates, defined as a collection of bonds, angles and
/// dihedral angles.
///
/// \tparam Vector3 3D vector type
/// \tparam Matrix Matrix type
/// \param n_atoms Total number of atoms
/// \param bonds Collection of bonds
/// \param angles Collection of angles between bonded atoms
/// \return Wilson's B matrix
///
/// This function returns Wilson's \f$\mathbf{B}\f$ matrix given a collection
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
                     const std::vector<connectivity::Angle<Vector3>>&
                        angles = {},
                     const std::vector<connectivity::Dihedral<Vector3>>&
                        dihedrals = {} ){
  // Get the total number of internal redundant coordinates
  size_t n_irc{ bonds.size() + angles.size() + dihedrals.size() };
  
  // Allocate Wilson's B matrix
  Matrix B{ linalg::zeros<Matrix>(n_irc, 3 * n_atoms) };
  
  // Utility vector for gradients storage
  Vector3 g1, g2, g3, g4;
  
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
  
  // Populate B matrix's rows corresponding to dihedrals
  connectivity::Dihedral<Vector3> dihedral;
  for(size_t i{0}; i < dihedrals.size(); i++){
    dihedral = dihedrals[i];
    
    std::tie(g1, g2, g3, g4) = dihedral_gradient(dihedral);
    
    for(size_t idx{0}; idx < 3; idx++){
      B(i + bonds.size() + angles.size(), 3 * dihedral.i + idx) = g1(idx);
      B(i + bonds.size() + angles.size(), 3 * dihedral.j + idx) = g2(idx);
      B(i + bonds.size() + angles.size(), 3 * dihedral.k + idx) = g3(idx);
      B(i + bonds.size() + angles.size(), 3 * dihedral.l + idx) = g4(idx);
    }
  }
  
  return B;
}


template <typename Vector3, typename Matrix>
Matrix wilson_matrix(const molecule::Molecule<Vector3>& molecule){
  // Compute interatomic distances
  Matrix dd{ connectivity::distances<Vector3, Matrix>(molecule) };
  
  // Compute adjacency matrix (graph)
  connectivity::UGraph adj{ connectivity::adjacency_matrix(dd, molecule) };
  
  // Compute distance matrix and predecessor matrix
  Matrix dist, predecessors;
  std::tie(dist, predecessors) = connectivity::distance_matrix<Matrix>(adj) ;
  
  // Compute bonds
  std::vector<connectivity::Bond<Vector3>> bonds{
      connectivity::bonds(dist, molecule)};
  
  // Compute angles
  std::vector<connectivity::Angle<Vector3>> angles{
      connectivity::angles(dist, predecessors, molecule)};
  
  // Compute dihedrals
  std::vector<connectivity::Dihedral<Vector3>> dihedrals{
      connectivity::dihedrals(dist, predecessors, molecule)};
  
  // Return Wilson's B matrix
  return wilson_matrix<Vector3, Matrix>(molecule.size(),
                                        bonds, angles, dihedrals);
}

/*
template <typename Vector3, typename Matrix>
Matrix wilson_matrix(
    const std::vector<std::tuple<size_t, double, double, double>>& atoms) {
  molecule::Molecule<Vector3> molecule{
      molecule::make_molecule<Vector3>(atoms)
  };
  
  return wilson_matrix<Vector3, Matrix>(molecule);
}
*/

template <typename Matrix>
std::pair<Matrix, Matrix> G_matirces(const Matrix& B){
  Matrix G{ B * linalg::transpose(B) };
  
  return std::make_pair( G, linalg::pseudo_inverse(G)  );
}

template <typename Matrix>
Matrix projector(const Matrix& G, const Matrix& iG){
  return G * iG;
}

} // namespace wilson

#endif //IRC_WILSON_H