#ifndef IRC_WILSON_H
#define IRC_WILSON_H

#include "connectivity/connectivity.h"

#include "atoms/molecule.h"
#include "connectivity/connectivity.h"
#include "tools/constants.h"
#include "tools/math.h"

#include <cmath>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

namespace irc {

namespace wilson {

template<typename Vector3>
std::pair<Vector3, Vector3>
bond_gradient(const Vector3& p1, const Vector3& p2) {
  double bond{ connectivity::distance(p1, p2) };
  
  Vector3 v{(p1 - p2) / bond};
  
  return {v, -v};
}

template<typename Vector3>
std::tuple<Vector3, Vector3, Vector3> angle_gradient(const Vector3& p1,
                                                     const Vector3& p2,
                                                     const Vector3& p3) {
  double angle{ connectivity::angle(p1, p2, p3) };
  
  double angle_rad{ tools::math::deg_to_rad(angle) };
  
  // TODO: Check pyberny for more robust implementation
  // https://github.com/azag0/pyberny
  
  double sin_angle{std::sin(angle_rad)};
  double cos_angle{std::cos(angle_rad)};
  
  Vector3 b21{p1 - p2};
  Vector3 b23{p3 - p2};
  
  double bond21{linalg::norm(b21)};
  double bond23{linalg::norm(b23)};
  
  b21 = b21 / bond21;
  b23 = b23 / bond23;
  
  Vector3 v1{(cos_angle * b21 - b23) / (sin_angle * bond21)};
  Vector3 v3{(cos_angle * b23 - b21) / (sin_angle * bond23)};
  Vector3 v2{-v1 - v3};
  
  return std::make_tuple(v1, v2, v3);
}

template<typename Vector3>
std::tuple<Vector3, Vector3, Vector3, Vector3> dihedral_gradient(
    const Vector3& p1,
    const Vector3& p2,
    const Vector3& p3,
    const Vector3& p4) {
  
  // TODO: Check pyberny for more robust implementation
  // https://github.com/azag0/pyberny
  
  double angle123{ connectivity::angle(p1, p2, p3) };
  double angle123_rad{tools::math::deg_to_rad(angle123)};
  double sin_angle123{std::sin(angle123_rad)};
  double cos_angle123{std::cos(angle123_rad)};
  
  double angle234{ connectivity::angle(p2, p3, p4) };
  double angle234_rad{tools::math::deg_to_rad(angle234)};
  double sin_angle234{std::sin(angle234_rad)};
  double cos_angle234{std::cos(angle234_rad)};
  
  Vector3 b12{p2 - p1};
  Vector3 b23{p3 - p2};
  Vector3 b34{p4 - p3};
  
  double bond12{linalg::norm(b12)};
  double bond23{linalg::norm(b23)};
  double bond34{linalg::norm(b34)};
  
  b12 = b12 / bond12;
  b23 = b23 / bond23;
  b34 = b34 / bond34;
  
  Vector3 b32{-b23};
  Vector3 b43{-b34};
  
  Vector3 v1{-linalg::cross(b12, b23) / (bond12 * sin_angle123 * sin_angle123)};
  
  double vc1{0.}, vc2{0.};
  Vector3 vv1, vv2{0., 0., 0.};
  
  vc1 = (bond23 - bond12 * cos_angle123) / (bond12 * bond23 * sin_angle123);
  vc2 = cos_angle234 / (bond23 * sin_angle234);
  vv1 = linalg::cross(b12, b23) / sin_angle123;
  vv2 = linalg::cross(b43, b32) / sin_angle234;
  
  Vector3 v2{vc1 * vv1 + vc2 * vv2};
  
  vc1 = (bond23 - bond34 * cos_angle234) / (bond23 * bond34 * sin_angle234);
  vc2 = cos_angle123 / (bond23 * sin_angle123);
  vv1 = linalg::cross(b43, b32) / sin_angle234;
  vv2 = linalg::cross(b12, b23) / sin_angle123;
  
  Vector3 v3{vc1 * vv1 + vc2 * vv2};
  
  Vector3 v4{-linalg::cross(b43, b32) / (bond34 * sin_angle234 * sin_angle234)};
  
  return std::make_tuple(v1, v2, v3, v4);
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
/// \patam Atomic positions in cartesian coordinates
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
template<typename Vector3, typename Vector, typename Matrix>
Matrix wilson_matrix(
    const Vector& x_cartesian,
    const std::vector<connectivity::Bond>& bonds,
    const std::vector<connectivity::Angle>& angles = {},
    const std::vector<connectivity::Dihedral>& dihedrals = {}){
  // Get number of atoms
  size_t n_atoms{ linalg::size<Vector>(x_cartesian) / 3 };

  // Get the total number of internal redundant coordinates
  size_t n_irc{ bonds.size() + angles.size() + dihedrals.size() };

  // Allocate Wilson's B matrix
  Matrix B{ linalg::zeros<Matrix>(n_irc, 3 * n_atoms) };

  // Utility vectors for atomic positions
  Vector3 p1, p2, p3, p4;

  // Utility vector for gradients storage
  Vector3 g1, g2, g3, g4;

  // B-matrix rows offset
  size_t offset{0};

  // Populate B matrix's rows corresponding to bonds
  connectivity::Bond bond;
  for(size_t i{0}; i < bonds.size(); i++) {
    bond = bonds[i];

    for(size_t m{0}; m < 3; m++){
      p1(m) = x_cartesian(3 * bond.i + m);
      p2(m) = x_cartesian(3 * bond.j + m);
    }

    std::tie(g1, g2) = bond_gradient(p1, p2);

    for(size_t idx{0}; idx < 3; idx++) {
      B(i, 3 * bond.i + idx) = g1(idx);
      B(i, 3 * bond.j + idx) = g2(idx);
    }
  }

  // Populate B matrix's rows corresponding to angles
  offset = bonds.size();
  connectivity::Angle angle;
  for(size_t i{0}; i < angles.size(); i++){
    angle = angles[i];

    for(size_t m{0}; m < 3; m++){
      p1(m) = x_cartesian(3 * angle.i + m);
      p2(m) = x_cartesian(3 * angle.j + m);
      p3(m) = x_cartesian(3 * angle.k + m);
    }

    std::tie(g1, g2, g3) = angle_gradient(p1, p2, p3);

    for(size_t idx{0}; idx < 3; idx++){
      B(i + offset, 3 * angle.i + idx) = g1(idx);
      B(i + offset, 3 * angle.j + idx) = g2(idx);
      B(i + offset, 3 * angle.k + idx) = g3(idx);

    }
  }

  // Populate B matrix's rows corresponding to dihedrals
  offset = bonds.size() + angles.size();
  connectivity::Dihedral dihedral;
  for(size_t i{0}; i < dihedrals.size(); i++){
    dihedral = dihedrals[i];

    for(size_t m{0}; m < 3; m++){
      p1(m) = x_cartesian(3 * dihedral.i + m);
      p2(m) = x_cartesian(3 * dihedral.j + m);
      p3(m) = x_cartesian(3 * dihedral.k + m);
      p4(m) = x_cartesian(3 * dihedral.l + m);
    }

    std::tie(g1, g2, g3, g4) = dihedral_gradient(p1, p2, p3, p4);

    for(size_t idx{0}; idx < 3; idx++){
      B(i + offset, 3 * dihedral.i + idx) = g1(idx);
      B(i + offset, 3 * dihedral.j + idx) = g2(idx);
      B(i + offset, 3 * dihedral.k + idx) = g3(idx);
      B(i + offset, 3 * dihedral.l + idx) = g4(idx);
    }
  }

  return B;
}

// TODO: Remove this function
template<typename Vector3, typename Vector, typename Matrix>
Matrix wilson_matrix(const molecule::Molecule <Vector3> &molecule) {
  // Compute interatomic distances
  Matrix dd{connectivity::distances<Vector3, Matrix>(molecule)};
  
  // Compute adjacency matrix (graph)
  connectivity::UGraph adj{connectivity::adjacency_matrix(dd, molecule)};
  
  // Compute distance matrix and predecessor matrix
  Matrix dist, predecessors;
  std::tie(dist, predecessors) = connectivity::distance_matrix<Matrix>(adj);
  
  // Compute bonds
  std::vector<connectivity::Bond>
  bonds{
      connectivity::bonds(dist, molecule)};
  
  // Compute angles
  std::vector<connectivity::Angle>
  angles{
      connectivity::angles(dist, predecessors, molecule)};
  
  // Compute dihedrals
  std::vector<connectivity::Dihedral>
  dihedrals{
      connectivity::dihedrals(dist, predecessors, molecule)};
  
  // Return Wilson's B matrix
  return wilson_matrix<Vector3, Vector, Matrix>(
      molecule::to_cartesian<Vector3, Vector>(molecule),
      bonds, angles, dihedrals);
}

template<typename Matrix>
Matrix projector(const Matrix &B) {
  // TODO: Pass iB instead of computing it
  return B * linalg::pseudo_inverse(B);
}

} // namespace wilson

} // namespace irc

#endif //IRC_WILSON_H
