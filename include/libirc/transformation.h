#ifndef IRC_TRANSFORMATION_H
#define IRC_TRANSFORMATION_H

#include "connectivity.h"
#include "linalg.h"
#include "mathtools.h"
#include "wilson.h"

#include <iostream>

namespace irc {

namespace transformation {

/// Compute the root mean square value of \param v
///
/// \tparam Vector Vector
/// \param v Vector
/// \return Root mean square value of \param v
template<typename Vector>
double rms(const Vector &v) {

  size_t size{linalg::size<Vector>(v)};

  double sum{0};

  // TODO: Use linalg::norm instead
  for (size_t i{0}; i < size; i++) {
    sum += v(i) * v(i);
  }

  return std::sqrt(sum / size);
}

// TODO: Avoid transpose?
/// Transform gradient from cartesian to internal redundant coordinates,
/// without projection
///
/// \tparam Vector
/// \tparam Matrix
/// \param grad_c Gradient in cartesian coordinates
/// \param B Wilson \f$\mathbf{B}\f$ matrix
/// \return Gradient in internal redundant coordinates
template<typename Vector, typename Matrix>
Vector gradient_cartesian_to_irc(const Vector &grad_c, const Matrix &B) {
  return linalg::pseudo_inverse(linalg::transpose(B)) * grad_c;
}

// TODO: Avoid transpose?
/// Transform gradient from internal redundant coordinates (non-projected)
/// to cartesian coordinates
///
/// \tparam Vector
/// \tparam Matrix
/// \param grad_irc Gradient in internal redundant coordinates
/// \param B Wilson \f$\mathbf{B}\f$ matrix
/// \return Gradient in cartesian coordinates
template<typename Vector, typename Matrix>
Vector gradient_irc_to_cartesian(const Vector &grad_irc, const Matrix &B) {
  return linalg::transpose(B) * grad_irc;
}

/// Transform internal redundant displacements to cartesian coordinates
///
/// \tparam Vector3
/// \tparam Vector
/// \tparam Matrix
/// \param q_irc_old Old internal reaction coordinates
/// \param dq_irc Change in internar reaction coordinates
/// \param x_c_old Old cartesian coordinates
/// \param bonds List of bonds
/// \param angles List of angles
/// \param dihedrals List of dihedral angles
/// \param B Wilson \f$\mathbf{B}\f$ matrix
/// \param iG Pseudoinverse of the \f$\mathbf{G}\f$ matrix
/// \param max_iters Maximum number of iterations
/// \param tolerance Tolerance on change in cartesian coordinates
/// \return New cartesian coordinates
///
/// Since Cartesian coordinates are rectilinear and the internal coordinates are
/// curvilinear, the transformation must be done iteratively.
template<typename Vector3, typename Vector, typename Matrix>
Vector irc_to_cartesian(const Vector &q_irc_old,
                        const Vector &dq_irc,
                        const Vector &x_c_old,
                        const std::vector<connectivity::Bond> &bonds,
                        const std::vector<connectivity::Angle> &angles,
                        const std::vector<connectivity::Dihedral> &dihedrals,
                        size_t max_iters = 25,
                        double tolerance = 1e-6) {
  // Number of internal redundant coordinates
  size_t n_irc{bonds.size() + angles.size() + dihedrals.size()};

  // Convergence flag
  bool converged{false};

  // Cartesian coordinates
  Vector x_c{x_c_old};
  
  // Store change in internal redundant coordinates
  Vector dq{dq_irc};

  // New internal coordinates
  Vector q_new{q_irc_old};

  // Change in cartesian coordinates
  Vector dx{linalg::zeros<Vector>(linalg::size(x_c_old))};

  // Compute Wilson's B matrix
  Matrix B{wilson::wilson_matrix<Vector3, Vector, Matrix>(
      x_c, bonds, angles, dihedrals)};

  // Compute the transpose of B
  Matrix iB{linalg::pseudo_inverse(B)};

  double RMS{0};

  // Start iterative search
  for (size_t i{0}; i < max_iters; i++) {

    // Compute displacement in cartesian coordinates
    dx = iB * dq;

    // Check for convergence
    RMS = rms<Vector>(dx);
    if (RMS < tolerance) {
      converged = true;
      break;
    }

    // Update cartesian coordinates
    x_c += dx;

    // Update Wilson B matrix
    //B = wilson::wilson_matrix<Vector3, Vector, Matrix>(
    //    x_c, bonds, angles, dihedrals);

    // Update transpose of the Wilson B matrix
    //iB = linalg::pseudo_inverse(B);

    // Compute new internal coordinates
    q_new = connectivity::cartesian_to_irc<Vector3, Vector>(x_c, bonds, angles, dihedrals);

    // Check change in dihedral angles (in radians)
    size_t offset{bonds.size() + angles.size()};
    for (size_t i{offset}; i < n_irc; i++) {
      // Restrain dihedral angle on the interval [-pi,pi]
      q_new(i) = tools::math::pirange_rad(q_new(i));
    }

    // New difference in internal coordinates
    dq = dq_irc - (q_new - q_irc_old);
  }

  // TODO: Store first iteration to avoid computation
  // If iteration does not converge, use first estimate
  if (!converged) {
    // Re-compute original B matrix
    B = wilson::wilson_matrix<Vector3, Vector, Matrix>(
        x_c_old, bonds, angles, dihedrals);

    // Compute first estimate
    x_c = x_c_old + linalg::pseudo_inverse(B) * dq_irc;

    // TODO: Something better?
    std::cerr << "WARNING: IRC_TO_CARTESIAN not converged." << std::endl;
  }

  return x_c;
}

} // namespace transformation

} // namespace irc

#endif // IRC_TRANSFORMATION_H