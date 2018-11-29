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
double rms(const Vector& v) {
  return linalg::norm(v) / std::sqrt(linalg::size(v));
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
Vector gradient_cartesian_to_irc(const Vector& grad_c, const Matrix& B) {
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
Vector gradient_irc_to_cartesian(const Vector& grad_irc, const Matrix& B) {
  return linalg::transpose(B) * grad_irc;
}

template<typename Vector>
struct IrcToCartesianResult {
  Vector x_c;
  bool converged;
  std::size_t n_iterations;
};

/// Transform internal redundant displacements to cartesian coordinates
///
/// \tparam Vector3
/// \tparam Vector
/// \tparam Matrix
/// \param q_irc_old Old internal redundant coordinates
/// \param dq_irc Change in internal redundant coordinates
/// \param x_c_old Old cartesian coordinates
/// \param bonds List of bonds
/// \param angles List of angles
/// \param dihedrals List of dihedral angles
/// \param max_iters Maximum number of iterations
/// \param tolerance Tolerance on change in cartesian coordinates
/// \return New cartesian coordinates
///
/// Since Cartesian coordinates are rectilinear and the internal coordinates are
/// curvilinear, the transformation must be done iteratively.
template<typename Vector3, typename Vector, typename Matrix>
IrcToCartesianResult<Vector>
irc_to_cartesian(const Vector& q_irc_old,
                 const Vector& dq_irc,
                 const Vector& x_c_old,
                 const std::vector<connectivity::Bond>& bonds,
                 const std::vector<connectivity::Angle>& angles,
                 const std::vector<connectivity::Dihedral>& dihedrals,
                 const std::vector<connectivity::LinearAngle<Vector3>>& linear_angles,
                 const std::vector<connectivity::OutOfPlaneBend>& out_of_plane_bends,
                 std::size_t max_iters = 25,
                 double tolerance = 1e-6) {

  bool converged{false};

  // Cartesian coordinates
  Vector x_c{x_c_old};

  // Store change in internal redundant coordinates
  Vector dq{dq_irc};

  // New internal coordinates
  Vector q_new{q_irc_old};

  // Change in cartesian coordinates
  Vector dx{linalg::zeros<Vector>(linalg::size(x_c_old))};

  // Wilson's B matrix
  const Matrix B{wilson::wilson_matrix<Vector3, Vector, Matrix>(
      x_c, bonds, angles, dihedrals, linear_angles, out_of_plane_bends)};

  // Transpose of B
  const Matrix iB{linalg::pseudo_inverse(B)};

  double RMS{0};

  // Offset for dihedral angles in q_irc
  const std::size_t offset{bonds.size() + angles.size()};

  std::size_t n_iterations{0};
  for (; n_iterations < max_iters; n_iterations++) {
    // Displacement in cartesian coordinates
    dx = iB * dq;

    // Check for convergence
    RMS = rms<Vector>(dx);
    if (RMS < tolerance) {
      converged = true;
      break;
    }

    // Update cartesian coordinates
    x_c += dx;

    // Compute new internal coordinates
    q_new = connectivity::cartesian_to_irc<Vector3, Vector>(
        x_c, bonds, angles, dihedrals, linear_angles, out_of_plane_bends);

    // Restrain dihedral angle on the interval (-pi,pi]
    for (std::size_t i{offset}; i < offset + dihedrals.size(); i++) {
      q_new(i) = tools::math::pirange_rad(q_new(i));
    }

    // New difference in internal coordinates
    dq = dq_irc - (q_new - q_irc_old);
  }

  // TODO: Store first iteration to avoid computation
  // If iteration does not converge, use first estimate
  if (!converged) {
    // Compute first estimate
    x_c = x_c_old + iB * dq_irc;
  }

  return {x_c, converged, n_iterations};
}

} // namespace transformation

} // namespace irc

#endif // IRC_TRANSFORMATION_H
