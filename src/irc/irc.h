#ifndef IRC_IRC_H
#define IRC_IRC_H

#include "../atoms/molecule.h"
#include "../connectivity/connectivity.h"
#include "../connectivity/wilson.h"
#include "../linear_algebra/linalg.h"
#include "../transformation/transformation.h"

#include <utility>

namespace irc {

template <typename Vector3, typename Vector, typename Matrix>
class IRC {
 public:
  IRC(const molecule::Molecule<Vector3>& molecule);
  
  Matrix projected_initial_hessian_inv(double alpha = 1000) const;
  
  Vector grad_cartesian_to_projected_irc(const Vector& grad_c) const;
  
  Vector irc_to_cartesian(const Vector& dq_irc,
                          const Vector& x_c_old,
                          size_t max_iters = 25);
  
 private:
  /// Distance matrix
  ///
  /// Contains information about connectivity and shortest paths within
  /// the initial molecule
  Matrix distance_m;
  
  /// Predecessor matrix
  ///
  /// Contains information to build the shortest path between atoms
  Matrix predecessors_m;
  
  /// List of bonds
  std::vector<connectivity::Bond<Vector3>> bonds;
  
  /// List of angles
  std::vector<connectivity::Angle<Vector3>> angles;
  
  /// List of dihedral angles
  std::vector<connectivity::Dihedral<Vector3>> dihedrals;
  
  /// Number of internal coordinates
  size_t n_irc;
  
  /// Number of cartesian coordinates
  size_t n_c;
  
  /// Wilson B matrix
  Matrix B;
  
  /// G matrix
  Matrix G;
  
  /// Generalized inverse of the G matrix
  Matrix iG;
  
  /// Projector
  Matrix P;
};

template <typename Vector3, typename Vector, typename Matrix>
IRC<Vector3, Vector, Matrix>::IRC(const molecule::Molecule<Vector3>& molecule){
  // Number of cartesian coordinates
  n_c = 3 * molecule.size();
  
  // Compute interatomic distances
  Matrix dd{ connectivity::distances<Vector3, Matrix>(molecule) };
  
  // Compute adjacency matrix (graph)
  connectivity::UGraph adj{ connectivity::adjacency_matrix(dd, molecule) };
  
  // Compute distance matrix and predecessor matrix
  std::tie(distance_m, predecessors_m) =
      connectivity::distance_matrix<Matrix>(adj) ;
  
  // Compute bonds
  bonds = connectivity::bonds(distance_m, molecule);
  
  // Compute angles
  angles = connectivity::angles(distance_m, predecessors_m, molecule);
  
  // Compute dihedrals
  dihedrals = connectivity::dihedrals(distance_m, predecessors_m, molecule);

  // Count the number of internal coordinates
  n_irc = bonds.size() + angles.size() + dihedrals.size();
  
  // Store initial Wilson's B matrix
  B = wilson::wilson_matrix<Vector3, Matrix>(molecule.size(),
                                             bonds, angles, dihedrals);
  
  // Compute G and iG (generalized inverse)
  std::tie(G, iG) = wilson::G_matirces(B);
  
  // Compute projector P
  P = wilson::projector(G, iG);
}

/// Initial estimate of the Hessian in internal redundant coordinates
///
/// \return
///
/// V. Bakken and T. Helgaker, J. Chem. Phys 117, 9160 (2002).
template <typename Vector3, typename Vector, typename Matrix>
Matrix IRC<Vector3, Vector, Matrix>::projected_initial_hessian_inv(
    double alpha) const {
  Matrix H0( linalg::zeros<Matrix>(n_irc, n_irc) );
  Matrix I( linalg::identity<Matrix>(n_irc) );
  
  size_t offset{0};
  
  for(size_t i{0}; i < bonds.size(); i++){
    H0(i,i) = 0.5;
  }
  
  offset = bonds.size();
  for(size_t i{0}; i < angles.size(); i++){
    H0(i + offset, i + offset) = 0.2;
  }
  
  offset = bonds.size() + angles.size();
  for(size_t i{0}; i < dihedrals.size(); i++){
    H0(i + offset, i + offset) = 0.1;
  }
  
  return linalg::inv<Matrix>( P * H0 * P + alpha * (I - P) );
}

/// Transform gradient in cartesian coordinates to gradient in internal
/// redundant coordinates and project the latter in the non-redundant
/// part of the internal coordinate space
///
/// \param grad_c Gradient in cartesian coordinates
/// \return Projected gradient in internal redundant coordinates
///
/// The gradient in redundant internal coordinates is given by
/// \f\[
///   \mathbf{g}_q = \mathbf{G}^-\mathbf{B} \mathbf{g}_x,
/// \f\]
/// where \f$\mathbf{g}_q\f$ and \f$\mathbf{g}_x\f$ are respectively the
/// gradient in internal redundant coordinates and the gradient in cartesian
/// coordinates. \f$\mathbf{B}\f$ is the Wilson B matrix, \f$\mathbf{G}\f$ is
/// the matrix defined by \f$\mathbf{G} = \mathbf{B}\mathbf{B}^T\f$ and
/// \f$\mathbf{G}^-\f$ is the pseudo-inverse of \f$\mathbf{G}\f$.
template <typename Vector3, typename Vector, typename Matrix>
Vector IRC<Vector3, Vector, Matrix>::grad_cartesian_to_projected_irc(
    const Vector& grad_c) const{
  return P *
      transformation::gradient_cartesian_to_irc<Vector,Matrix>(grad_c, B, iG);
}

template <typename Vector3, typename Vector, typename Matrix>
Vector IRC<Vector3, Vector, Matrix>::irc_to_cartesian(
    const Vector& dq_irc,
    const Vector& x_c_old,
    size_t max_iters){
  
}

}

#endif //IRC_IRC_H
