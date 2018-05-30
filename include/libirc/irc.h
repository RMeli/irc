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
  IRC(const molecule::Molecule<Vector3>& molecule = {},
      const std::vector<connectivity::Bond>& mybonds = {},
      const std::vector<connectivity::Angle>& myangles = {},
      const std::vector<connectivity::Dihedral>& mydihedrals = {});

  /// Compute initial projected inverted Hessian estimate
  ///
  /// \return Projected inverted initial Hessian
  Matrix projected_initial_hessian_inv(double k_bond = 0.5,
                                       double k_angle = 0.2,
                                       double k_dihedral = 0.1) const;

  /// Compute initial projected Hessian estimate
  ///
  /// \return Projected initial Hessian
  Matrix projected_initial_hessian(double k_bond = 0.5,
                                   double k_angle = 0.2,
                                   double k_dihedral = 0.1) const;

  /// Project inverted Hessian
  ///
  /// \param Hinv Inverted Hessian
  /// \return Projected inverted Hessian
  Matrix projected_hessian_inv(const Matrix& Hinv) const;
  
  /// Project Hessian
  ///
  /// \param H Hessian
  /// \return Projected Hessian
  Matrix projected_hessian(const Matrix& H) const;

  /// Transform gradient from cartesian coordinates to projected redundant
  /// internal coordinates.
  ///
  /// \param grad_c Gradient in cartesian coordinates
  /// \return Gradient in redundant internal coordinates
  Vector grad_cartesian_to_projected_irc(const Vector& grad_c) const;

  /// Transform cartesian coordinates to redundant internal coordinates
  ///
  /// \param x_c Cartesian coordinates
  /// \return Redundant internal coordinates
  Vector cartesian_to_irc(const Vector& x_c) const;

  /// Tranform redundant internal coordinates to cartesian coordinates
  ///
  /// \param q_irc_old Old redundant internal coordinates
  /// \param dq_irc Change in redundant internal coordinates
  /// \param x_c_old Old cartesian coordinates
  /// \param max_iters Maximum number of iterations
  /// \param tolerance Convergence tolerance
  /// \return New cartesian coordinates
  Vector irc_to_cartesian(const Vector& q_irc_old,
                          const Vector& dq_irc,
                          const Vector& x_c_old,
                          std::size_t max_iters = 25,
                          double tolerance = 1e-6);

private:
  /// List of bonds
  std::vector<connectivity::Bond> bonds;

  /// List of angles
  std::vector<connectivity::Angle> angles;

  /// List of dihedral angles
  std::vector<connectivity::Dihedral> dihedrals;

  /// Number of internal coordinates
  std::size_t n_irc;

  /// Number of cartesian coordinates
  std::size_t n_c;

  /// Wilson B matrix
  Matrix B;

  /// Projector
  Matrix P;
};

template<typename Vector3, typename Vector, typename Matrix>
IRC<Vector3, Vector, Matrix>::IRC(
    const molecule::Molecule<Vector3>& molecule,
    const std::vector<connectivity::Bond>& mybonds,
    const std::vector<connectivity::Angle>& myangles,
    const std::vector<connectivity::Dihedral>& mydihedrals) {

  // Number of cartesian coordinates
  n_c = 3 * molecule.size();

  // Compute interatomic distances
  const Matrix dd{connectivity::distances<Vector3, Matrix>(molecule)};

  // Compute adjacency matrix (graph)
  const connectivity::UGraph adj{connectivity::adjacency_matrix(dd, molecule)};

  // Compute distance matrix and predecessor matrix
  Matrix distance_m{connectivity::distance_matrix<Matrix>(adj)};

  // Compute bonds
  bonds = connectivity::bonds(distance_m, molecule);

  // Add user-defined bonds
  if (!mybonds.empty()) { // For CodeCov, can be removed after tests
    bonds.insert(bonds.cend(), mybonds.cbegin(), mybonds.cend());
  }

  // Compute angles
  angles = connectivity::angles(distance_m, molecule);

  // Add user-defined angles
  if (!myangles.empty()) { // For CodeCov, can be removed after tests
    angles.insert(angles.cend(), myangles.cbegin(), myangles.cend());
  }

  // Compute dihedrals
  dihedrals = connectivity::dihedrals(distance_m, molecule);

  // Add user-defined dihedrals
  if (!mydihedrals.empty()) { // For CodeCov, can be removed after tests
    dihedrals.insert(
        dihedrals.cend(), mydihedrals.cbegin(), mydihedrals.cend());
  }

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

/// Initial estimate of the inverse Hessian in internal redundant coordinates
///
/// \return
///
/// V. Bakken and T. Helgaker, J. Chem. Phys 117, 9160 (2002).
template<typename Vector3, typename Vector, typename Matrix>
Matrix IRC<Vector3, Vector, Matrix>::projected_initial_hessian_inv(
    double k_bond,
    double k_angle,
    double k_dihedral) const {
  Matrix iH0(linalg::zeros<Matrix>(n_irc, n_irc));

  std::size_t offset{0};

  for (std::size_t i{0}; i < bonds.size(); i++) {
    iH0(i, i) = 1. / k_bond;
  }

  offset = bonds.size();
  for (std::size_t i{0}; i < angles.size(); i++) {
    iH0(i + offset, i + offset) = 1. / k_angle;
  }

  offset = bonds.size() + angles.size();
  for (std::size_t i{0}; i < dihedrals.size(); i++) {
    iH0(i + offset, i + offset) = 1. / k_dihedral;
  }

  return P * iH0 * P;
}

/// Initial estimate of the Hessian in internal redundant coordinates
///
/// \return
///
/// V. Bakken and T. Helgaker, J. Chem. Phys 117, 9160 (2002).
template<typename Vector3, typename Vector, typename Matrix>
Matrix IRC<Vector3, Vector, Matrix>::projected_initial_hessian(
    double k_bond,
    double k_angle,
    double k_dihedral) const {
  Matrix H0(linalg::zeros<Matrix>(n_irc, n_irc));

  std::size_t offset{0};

  for (std::size_t i{0}; i < bonds.size(); i++) {
    H0(i, i) = k_bond;
  }

  offset = bonds.size();
  for (std::size_t i{0}; i < angles.size(); i++) {
    H0(i + offset, i + offset) = k_angle;
  }

  offset = bonds.size() + angles.size();
  for (std::size_t i{0}; i < dihedrals.size(); i++) {
    H0(i + offset, i + offset) = k_dihedral;
  }

  return P * H0 * P;
}

template<typename Vector3, typename Vector, typename Matrix>
Matrix
IRC<Vector3, Vector, Matrix>::projected_hessian_inv(const Matrix& Hinv) const {

  if (linalg::size(Hinv) != n_irc * n_irc) {
    throw std::length_error("ERROR: Wrong Hessian size.");
  }

  return P * Hinv * P;
}

template<typename Vector3, typename Vector, typename Matrix>
Matrix
IRC<Vector3, Vector, Matrix>::projected_hessian(const Matrix& H) const {
  
  if (linalg::size(H) != n_irc * n_irc) {
    throw std::length_error("ERROR: Wrong Hessian size.");
  }
  
  return P * H * P;
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
    const Vector& grad_c) const {
  if (linalg::size(grad_c) != n_c) {
    throw std::length_error("ERROR: Wrong cartesian gradient size.");
  }

  return P *
         transformation::gradient_cartesian_to_irc<Vector, Matrix>(grad_c, B);
}

template<typename Vector3, typename Vector, typename Matrix>
Vector IRC<Vector3, Vector, Matrix>::cartesian_to_irc(const Vector& x_c) const {
  if (linalg::size(x_c) != n_c) {
    throw std::length_error("ERROR: Wrong cartesian coordinates size.");
  }

  return connectivity::cartesian_to_irc<Vector3, Vector>(
      x_c, bonds, angles, dihedrals);
}

template<typename Vector3, typename Vector, typename Matrix>
Vector IRC<Vector3, Vector, Matrix>::irc_to_cartesian(const Vector& q_irc_old,
                                                      const Vector& dq_irc,
                                                      const Vector& x_c_old,
                                                      std::size_t max_iters,
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

  const auto itc_result =
      transformation::irc_to_cartesian<Vector3, Vector, Matrix>(q_irc_old,
                                                                dq_irc,
                                                                x_c_old,
                                                                bonds,
                                                                angles,
                                                                dihedrals,
                                                                max_iters,
                                                                tolerance);

  // TODO: This computation can be avoided; B is computed in irc_to_cartesian
  // Update Wilson's B matrix
  B = wilson::wilson_matrix<Vector3, Vector, Matrix>(
      itc_result.x_c, bonds, angles, dihedrals);

  // Update projector P
  P = wilson::projector(B);

  // Return new cartesian coordinates
  return itc_result.x_c;
}

} // namespace irc

#endif // IRC_IRC_H
