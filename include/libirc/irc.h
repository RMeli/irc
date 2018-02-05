#ifndef IRC_IRC_H
#define IRC_IRC_H

#include "connectivity.h"
#include "linalg.h"
#include "molecule.h"
#include "transformation.h"
#include "wilson.h"

#include <utility>

namespace irc {

template<typename Vector3, typename Vector, typename Matrix>
class IRC {
public:
  IRC(const molecule::Molecule<Vector3> &molecule = {});

  Matrix projected_initial_hessian_inv() const;

  Matrix projected_hessian_inv(const Matrix &H) const;

  Vector grad_cartesian_to_projected_irc(const Vector &grad_c) const;

  Vector cartesian_to_irc(const Vector &x_c) const;

  Vector irc_to_cartesian(const Vector &q_irc_old,
                          const Vector &dq_irc,
                          const Vector &x_c_old,
                          size_t max_iters = 25,
                          double tolerance = 1e-6);

private:
  /// List of bonds
  std::vector<connectivity::Bond> bonds;

  /// List of angles
  std::vector<connectivity::Angle> angles;

  /// List of dihedral angles
  std::vector<connectivity::Dihedral> dihedrals;

  /// Number of internal coordinates
  size_t n_irc;

  /// Number of cartesian coordinates
  size_t n_c;

  /// Wilson B matrix
  Matrix B;

  /// Projector
  Matrix P;
};

template<typename Vector3, typename Vector, typename Matrix>
IRC<Vector3, Vector, Matrix>::IRC(const molecule::Molecule<Vector3> &molecule) {
  // Number of cartesian coordinates
  n_c = 3 * molecule.size();

  // Compute interatomic distances
  Matrix dd{connectivity::distances<Vector3, Matrix>(molecule)};

  // Compute adjacency matrix (graph)
  connectivity::UGraph adj{connectivity::adjacency_matrix(dd, molecule)};

  // Compute distance matrix and predecessor matrix
  Matrix distance_m, predecessors_m;
  std::tie(distance_m, predecessors_m) =
      connectivity::distance_matrix<Matrix>(adj);

  // Compute bonds
  bonds = connectivity::bonds(distance_m, molecule);

  // Compute angles
  angles = connectivity::angles(distance_m, predecessors_m, molecule);

  // Compute dihedrals
  dihedrals = connectivity::dihedrals(distance_m, predecessors_m, molecule);

  // Count the number of internal coordinates
  n_irc = bonds.size() + angles.size() + dihedrals.size();

  // Store initial Wilson's B matrix
  B = wilson::wilson_matrix<Vector3, Vector, Matrix>(
      molecule::to_cartesian<Vector3, Vector>(molecule),
      bonds,
      angles,
      dihedrals);

  // Compute projector P
  P = wilson::projector(B);
}

/// Initial estimate of the Hessian in internal redundant coordinates
///
/// \return
///
/// V. Bakken and T. Helgaker, J. Chem. Phys 117, 9160 (2002).
template<typename Vector3, typename Vector, typename Matrix>
Matrix IRC<Vector3, Vector, Matrix>::projected_initial_hessian_inv() const {
  Matrix H0(linalg::zeros<Matrix>(n_irc, n_irc));

  size_t offset{0};

  for (size_t i{0}; i < bonds.size(); i++) {
    H0(i, i) = 0.5;
  }

  offset = bonds.size();
  for (size_t i{0}; i < angles.size(); i++) {
    H0(i + offset, i + offset) = 0.2;
  }

  offset = bonds.size() + angles.size();
  for (size_t i{0}; i < dihedrals.size(); i++) {
    H0(i + offset, i + offset) = 0.1;
  }

  return P * linalg::inv<Matrix>(H0) * P;
}

template<typename Vector3, typename Vector, typename Matrix>
Matrix
IRC<Vector3, Vector, Matrix>::projected_hessian_inv(const Matrix &Hinv) const {

  if (linalg::size(Hinv) != n_irc * n_irc) {
    throw std::length_error("ERROR: Wrong Hessian size.");
  }

  return P * Hinv * P;
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
template<typename Vector3, typename Vector, typename Matrix>
Vector IRC<Vector3, Vector, Matrix>::grad_cartesian_to_projected_irc(
    const Vector &grad_c) const {
  if (linalg::size(grad_c) != n_c) {
    throw std::length_error("ERROR: Wrong cartesian gradient size.");
  }

  return P *
         transformation::gradient_cartesian_to_irc<Vector, Matrix>(grad_c, B);
}

template<typename Vector3, typename Vector, typename Matrix>
Vector IRC<Vector3, Vector, Matrix>::cartesian_to_irc(const Vector &x_c) const {
  if (linalg::size(x_c) != n_c) {
    throw std::length_error("ERROR: Wrong cartesian coordinates size.");
  }

  return connectivity::cartesian_to_irc<Vector3, Vector>(
      x_c, bonds, angles, dihedrals);
}

template<typename Vector3, typename Vector, typename Matrix>
Vector IRC<Vector3, Vector, Matrix>::irc_to_cartesian(const Vector &q_irc_old,
                                                      const Vector &dq_irc,
                                                      const Vector &x_c_old,
                                                      size_t max_iters,
                                                      double tolerance) {

  if (linalg::size(q_irc_old) != n_irc) {
    throw std::length_error("ERROR: Wrong old IRC coordinates size.");
  }

  if (linalg::size(dq_irc) != n_irc) {
    throw std::length_error("ERROR: Wrong IRC displacement size.");
  }

  if (linalg::size(x_c_old) != n_c) {
    throw std::length_error("ERROR: Wrong old cartesian coordinates size.");
  }

  Vector x_c_new{
      transformation::irc_to_cartesian<Vector3, Vector, Matrix>(q_irc_old,
                                                                dq_irc,
                                                                x_c_old,
                                                                bonds,
                                                                angles,
                                                                dihedrals,
                                                                max_iters,
                                                                tolerance)};

  // TODO: This computation can be avoided; B is computed in irc_to_cartesian
  // Update Wilson's B matrix
  B = wilson::wilson_matrix<Vector3, Vector, Matrix>(
      x_c_new, bonds, angles, dihedrals);

  // Update projector P
  P = wilson::projector(B);

  // Return new cartesian coordinates
  return x_c_new;
}

} // namespace irc

#endif // IRC_IRC_H
