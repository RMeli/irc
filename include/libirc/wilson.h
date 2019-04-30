#ifndef IRC_WILSON_H
#define IRC_WILSON_H

#include "connectivity.h"

#include "constants.h"
#include "mathtools.h"
#include "molecule.h"

#include <cmath>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

namespace irc {

namespace wilson {

/// Compute bond gradients
///
/// A pair of vectors act to increase the distance between \p p1 and \p p2
/// when added to their respective cartesian coordinates.
/// The displacement vectors are:
/// \f[
///    \left(\frac{p_1 - p_2}{d}, -\frac{p_1 - p_2}{d}\right)
/// \f]
/// where \f$d = \lVert p_2 - p_1\rVert \f$.
///
/// \tparam Vector3
/// \param p1 Point 1
/// \param p2 Point 2
/// \return A pair of cartesian displacements
template<typename Vector3>
std::pair<Vector3, Vector3> bond_gradient(const Vector3& p1, ///< Point 1
                                          const Vector3& p2) {
  const double d{connectivity::distance(p1, p2)};

  const Vector3 u{(p1 - p2) / d};

  return {u, -u};
}

/*! Compute angle gradients
 *
 * Three vectors act to increase the angle between \p p1, \p p2 and \p p3
 * when added to their respective cartesian coordinates.
 *
 * Bakken and Helgaker, J. Chem. Phys., 117, 9160 (2002).
 *
 * \tparam Vector3
 * \param p1 Point 1
 * \param p2 Point 2
 * \param p3 Point 3
 * \return Angle gradients
 */
template<typename Vector3>
std::tuple<Vector3, Vector3, Vector3> angle_gradient(const Vector3& p1,
                                                     const Vector3& p2,
                                                     const Vector3& p3,
                                                     double tolerance = 1e-6) {
  using tools::math::collinear;

  const double angle{connectivity::angle(p1, p2, p3)};

  Vector3 u{p1 - p2};
  Vector3 v{p3 - p2};

  const double bond21{linalg::norm(u)};
  const double bond23{linalg::norm(v)};

  u = u / bond21;
  v = v / bond23;

  // Deal with linear angles
  Vector3 w;
  Vector3 pmp{1, -1, 1}, mpp{-1, 1, 1};
  if (std::abs(angle - tools::constants::pi) > tolerance) {
    w = linalg::cross(u, v);
  } else if (collinear(u, pmp, tolerance) && collinear(v, pmp, tolerance)) {
    w = linalg::cross(u, pmp);
  } else if (collinear(u, mpp, tolerance) && collinear(v, mpp, tolerance)) {
    w = linalg::cross(u, mpp);
  } else {
    throw std::runtime_error("Problem with linear angle.");
  }

  w = w / linalg::norm(w);

  const Vector3 v1{linalg::cross(u, w) / bond21};
  const Vector3 v3{linalg::cross(w, v) / bond23};
  const Vector3 v2{-v1 - v3};

  return std::make_tuple(v1, v2, v3);
}

/*! Compute dihedral angle gradients
 *
 * Four vectors act to increase the dihedral angle between \p p1, \p p2, \p p3,
 * and \p p4 when added to their respective cartesian coordinates.
 * The dihedral is the rotation about \f$(p_3 - p_2)\f$ that maps \p p1 on to
 * \p p4 when projected on to a plane with a normal vector \f$(p_3 - p_2)\f$.
 * The displacement vectors are:
 * \f{eqnarray*}{
 *    v_1 &=& - \frac{b_{12} \times b_{23}}{b_{12} \sin^2 \phi_2} \\
 *    v_2 &=& \frac{b_{23} - b_{12} \cos \phi_2 }{b_{23} b_{12} \sin \phi_2}
 *            \frac{b_{12} \times b_{23}}{\sin \phi_2}
 *            + \frac{\cos \phi_3 }{b_{23} \sin \phi_3}
 *              \frac{b_{43} \times b_{32}}{\sin \phi_3} \\
 *    v_3 &=& \frac{b_{23} - b_{43} \cos \phi_3 }{b_{32} b_{43} \sin \phi_3}
 *            \frac{b_{43} \times b_{32}}{\sin \phi_3}
 *            + \frac{\cos \phi_2 }{b_{32} \sin \phi_2}
 *              \frac{b_{12} \times b_{23}}{\sin \phi_2} \\
 *    v_4 &=& - \frac{b_{43} \times b_{32}}{b_{43} \sin^2 \phi_3}
 * \f}
 * where \f$d_{ij} = \lVert p_i - p_j\rVert\f$,
 * \f$b_{ij} = \frac{p_i - p_j}{d_{ij}}\f$,
 * \f$\phi_2 is the angle between \f$(p_1, p_2, p_3)\f$ and
 * \f$\phi_3 is the angle between \f$(p_2, p_3, p_4)\f$.
 *
 * \tparam Vector3
 * \param p1 Point 1
 * \param p2 Point 2
 * \param p3 Point 3
 * \param p3 Point 4
 * \return Dihedral angle gradients
 */
template<typename Vector3>
std::tuple<Vector3, Vector3, Vector3, Vector3>
dihedral_gradient(const Vector3& p1,
                  const Vector3& p2,
                  const Vector3& p3,
                  const Vector3& p4) {

  const double angle123{connectivity::angle(p1, p2, p3)};
  const double sin_angle123{std::sin(angle123)};
  const double cos_angle123{std::cos(angle123)};

  const double angle234{connectivity::angle(p2, p3, p4)};
  const double sin_angle234{std::sin(angle234)};
  const double cos_angle234{std::cos(angle234)};

  Vector3 b12{p2 - p1};
  Vector3 b23{p3 - p2};
  Vector3 b34{p4 - p3};

  const double bond12{linalg::norm(b12)};
  const double bond23{linalg::norm(b23)};
  const double bond34{linalg::norm(b34)};

  b12 = b12 / bond12;
  b23 = b23 / bond23;
  b34 = b34 / bond34;

  const Vector3 b32{-b23};
  const Vector3 b43{-b34};

  const Vector3 v1{-linalg::cross(b12, b23) /
                   (bond12 * sin_angle123 * sin_angle123)};

  double vc1{0.}, vc2{0.};
  Vector3 vv1, vv2{0., 0., 0.};

  vc1 = (bond23 - bond12 * cos_angle123) / (bond12 * bond23 * sin_angle123);
  vc2 = cos_angle234 / (bond23 * sin_angle234);
  vv1 = linalg::cross(b12, b23) / sin_angle123;
  vv2 = linalg::cross(b43, b32) / sin_angle234;

  const Vector3 v2{vc1 * vv1 + vc2 * vv2};

  vc1 = (bond23 - bond34 * cos_angle234) / (bond23 * bond34 * sin_angle234);
  vc2 = cos_angle123 / (bond23 * sin_angle123);
  vv1 = linalg::cross(b43, b32) / sin_angle234;
  vv2 = linalg::cross(b12, b23) / sin_angle123;

  const Vector3 v3{vc1 * vv1 + vc2 * vv2};

  const Vector3 v4{-linalg::cross(b43, b32) /
                   (bond34 * sin_angle234 * sin_angle234)};

  return std::make_tuple(v1, v2, v3, v4);
}


/*! Compute out of plane angle gradients
 */
template<typename Vector3>
std::tuple<Vector3, Vector3, Vector3, Vector3>
out_of_plane_gradient(const Vector3& vc,
                      const Vector3& v1,
                      const Vector3& v2,
                      const Vector3& v3) {
  const Vector3 b1{v1 - vc};
  const Vector3 b2{v2 - vc};
  const Vector3 b3{v3 - vc};

  const Vector3 e1{linalg::normalize(b1)};
  const Vector3 e2{linalg::normalize(b2)};
  const Vector3 e3{linalg::normalize(b3)};

  const double a1{connectivity::angle(v2, vc, v3)};
  const double a2{connectivity::angle(v3, vc, v1)};
  const double a3{connectivity::angle(v1, vc, v2)};

  const double sin_a1{std::sin(a1)};
  const double cos_a1{std::cos(a1)};
  const double cos_a2{std::cos(a2)};
  const double cos_a3{std::cos(a3)};

  const double ir1{1/linalg::norm(b1)};
  const double ir2{1/linalg::norm(b2)};
  const double ir3{1/linalg::norm(b3)};

  const Vector3 t1{linalg::cross(e2, e3)/std::sin(a1)};

  // Out of plane angle
  const double angle = std::asin(linalg::dot(t1, e1));

  const double cos_angle{std::cos(angle)};
  const double tan_angle{std::tan(angle)};

  const Vector3 s1 = ir1 * (t1/cos_angle - tan_angle * e1);
  const double denominator = cos_angle * sin_a1 * sin_a1;
  const Vector3 s2 = ir2 * t1 * (cos_a1 * cos_a2 - cos_a3)/denominator;
  const Vector3 s3 = ir3 * t1 * (cos_a1 * cos_a3 - cos_a2)/denominator;
  const Vector3 sc = - s1 - s2 - s3;

  return std::make_tuple(sc, s1, s2, s3);
}

/*!
 * \brief The linear angle gradient bending in the plane with the \p
 * othogonal_direction
 *
 * The gradient contribution \f$ (g_1, g_2, g_3) \f$ is computed at positions
 * \p p1, \p p2 and \p p3. These are formed by first computing the standard
 * angle gradients \f$ (g_1, g_{2a}, g_{0a}) \f$ at the positions
 * \f$ (p1, p2, p2+d_{orth}) \f$ and \f$ (g_{0b}, g_{2b}, g_3) \f$ at the
 * positions \f$ (p2+d_{orth}, p2, p3) \f$. The linear angle gradient is then
 * formed as \f$ (g_1, g_{2a}+g_{2b}, g_3) \f$
 *
 * \tparam Vector3
 * \param p1 Point 1
 * \param p2 Point 2
 * \param p3 Point 3
 * \param orthogonal_direction orthogonal director to the \p p1 to \p p3 vector
 * \param tolerance
 * \return
 */
template<typename Vector3>
std::tuple<Vector3, Vector3, Vector3>
linear_angle_gradient(const Vector3& p1,
                      const Vector3& p2,
                      const Vector3& p3,
                      const Vector3& orthogonal_direction,
                      double tolerance = 1e-6) {

  Vector3 v1, v2, v2add, v3, vOrth;

  const Vector3 p0 = p2 + orthogonal_direction;

  std::tie(v1, v2, vOrth) = angle_gradient(p1, p2, p0, tolerance);
  std::tie(vOrth, v2add, v3) = angle_gradient(p0, p2, p3, tolerance);

  return std::make_tuple(v1, v2 + v2add, v3);
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
/// defines the transformation from Cartesian displacements
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
    const std::vector<connectivity::Dihedral>& dihedrals = {},
    const std::vector<connectivity::LinearAngle<Vector3>>& linear_angles = {},
              const std::vector<connectivity::OutOfPlaneBend>& out_of_plane_bends = {}) {
  const std::size_t n_atoms{linalg::size<Vector>(x_cartesian) / 3};

  const std::size_t n_irc{bonds.size() + angles.size() + dihedrals.size() +
                          linear_angles.size() + out_of_plane_bends.size()};

  // Wilson's B matrix
  Matrix B{linalg::zeros<Matrix>(n_irc, 3 * n_atoms)};

  // Utility vectors for atomic positions
  Vector3 p1, p2, p3, p4;

  // Utility vector for gradients storage
  Vector3 g1, g2, g3, g4;

  // B-matrix rows offset
  std::size_t offset{0};

  // Populate B matrix's rows corresponding to bonds
  for (std::size_t i{0}; i < bonds.size(); i++) {
    auto bond = bonds[i];

    for (std::size_t m{0}; m < 3; m++) {
      p1(m) = x_cartesian(3 * bond.i + m);
      p2(m) = x_cartesian(3 * bond.j + m);
    }

    std::tie(g1, g2) = bond_gradient(p1, p2);

    for (std::size_t idx{0}; idx < 3; idx++) {
      B(i, 3 * bond.i + idx) = g1(idx);
      B(i, 3 * bond.j + idx) = g2(idx);
    }
  }

  // Populate B matrix's rows corresponding to angles
  offset = bonds.size();
  for (std::size_t i{0}; i < angles.size(); i++) {
    auto angle = angles[i];

    for (std::size_t m{0}; m < 3; m++) {
      p1(m) = x_cartesian(3 * angle.i + m);
      p2(m) = x_cartesian(3 * angle.j + m);
      p3(m) = x_cartesian(3 * angle.k + m);
    }

    std::tie(g1, g2, g3) = angle_gradient(p1, p2, p3);

    for (std::size_t idx{0}; idx < 3; idx++) {
      B(i + offset, 3 * angle.i + idx) = g1(idx);
      B(i + offset, 3 * angle.j + idx) = g2(idx);
      B(i + offset, 3 * angle.k + idx) = g3(idx);
    }
  }

  // Populate B matrix's rows corresponding to dihedrals
  offset = bonds.size() + angles.size();
  for (std::size_t i{0}; i < dihedrals.size(); i++) {
    auto dihedral = dihedrals[i];

    for (std::size_t m{0}; m < 3; m++) {
      p1(m) = x_cartesian(3 * dihedral.i + m);
      p2(m) = x_cartesian(3 * dihedral.j + m);
      p3(m) = x_cartesian(3 * dihedral.k + m);
      p4(m) = x_cartesian(3 * dihedral.l + m);
    }

    std::tie(g1, g2, g3, g4) = dihedral_gradient(p1, p2, p3, p4);

    for (std::size_t idx{0}; idx < 3; idx++) {
      B(i + offset, 3 * dihedral.i + idx) = g1(idx);
      B(i + offset, 3 * dihedral.j + idx) = g2(idx);
      B(i + offset, 3 * dihedral.k + idx) = g3(idx);
      B(i + offset, 3 * dihedral.l + idx) = g4(idx);
    }
  }

  // Populate B matrix's rows corresponding to linear angles
  offset = bonds.size() + angles.size() + dihedrals.size();
  for (std::size_t i{0}; i < linear_angles.size(); i++) {
    auto linear_angle = linear_angles[i];

    for (std::size_t m{0}; m < 3; m++) {
      p1(m) = x_cartesian(3 * linear_angle.i + m);
      p2(m) = x_cartesian(3 * linear_angle.j + m);
      p3(m) = x_cartesian(3 * linear_angle.k + m);
    }

    std::tie(g1, g2, g3) =
        linear_angle_gradient(p1, p2, p3, linear_angle.orthogonal_direction);

    for (std::size_t idx{0}; idx < 3; idx++) {
      B(i + offset, 3 * linear_angle.i + idx) = g1(idx);
      B(i + offset, 3 * linear_angle.j + idx) = g2(idx);
      B(i + offset, 3 * linear_angle.k + idx) = g3(idx);
    }
  }

  // Populate B matrix's rows corresponding to out of plane bends
  offset = bonds.size() + angles.size() + dihedrals.size() + linear_angles.size();
  for (std::size_t i{0}; i < out_of_plane_bends.size(); i++) {
    auto bend = out_of_plane_bends[i];

    for (std::size_t m{0}; m < 3; m++) {
      p4(m) = x_cartesian(3 * bend.c + m);
      p1(m) = x_cartesian(3 * bend.i + m);
      p2(m) = x_cartesian(3 * bend.j + m);
      p3(m) = x_cartesian(3 * bend.k + m);
    }

    std::tie(g4, g1, g2, g3) = out_of_plane_gradient(p4, p1, p2, p3);

    for (std::size_t idx{0}; idx < 3; idx++) {
      B(i + offset, 3 * bend.c + idx) = g4(idx);
      B(i + offset, 3 * bend.i + idx) = g1(idx);
      B(i + offset, 3 * bend.j + idx) = g2(idx);
      B(i + offset, 3 * bend.k + idx) = g3(idx);
    }
  }
  return B;
}

template<typename Vector3, typename Vector, typename Matrix>
Matrix wilson_matrix_numerical(
    const Vector& x_c,
    const std::vector<connectivity::Bond>& bonds,
    const std::vector<connectivity::Angle>& angles = {},
    const std::vector<connectivity::Dihedral>& dihedrals = {},
    const std::vector<connectivity::LinearAngle<Vector3>>& linear_angles = {},
    const std::vector<connectivity::OutOfPlaneBend>& out_of_plane_bends = {},
    double dx = 1.e-6) {

  const std::size_t n_c{linalg::size(x_c)};

  const std::size_t n_irc{bonds.size() + angles.size() + dihedrals.size() +
                          linear_angles.size() + out_of_plane_bends.size()};

  // Wilson B matrix
  Matrix B{linalg::zeros<Matrix>(n_irc, n_c)};

  // Displaced cartesian coordinates
  Vector x_c_pm{x_c};

  // Displaced IRC
  Vector q_irc_plus{linalg::zeros<Vector>(n_irc)};
  Vector q_irc_minus{linalg::zeros<Vector>(n_irc)};

  for (std::size_t j{0}; j < n_c; j++) {
    // Positive displacement for cartesian coordinate j
    x_c_pm(j) += dx;

    // IRC corresponding to positive displacement of x_c(j)
    q_irc_plus = connectivity::cartesian_to_irc<Vector3, Vector>(
        x_c_pm, bonds, angles, dihedrals, linear_angles, out_of_plane_bends);

    // Negative displacement for cartesian coordinate j
    x_c_pm(j) -= 2 * dx;

    // IRC corresponding to negative displacement of x_c(j)
    q_irc_minus = connectivity::cartesian_to_irc<Vector3, Vector>(
        x_c_pm, bonds, angles, dihedrals, linear_angles, out_of_plane_bends);

    for (std::size_t i{0}; i < n_irc; i++) {
      // Compute derivative (centered finite difference)
      B(i, j) = (q_irc_plus(i) - q_irc_minus(i)) / (2 * dx);
    }

    // Reset original cartesian coordinates
    x_c_pm(j) = x_c(j);
  }

  return B;
}

/// Compute projector from the \param B matrix
/// \tparam Matrix
/// \param B Wilson's B matrix
/// \return Projector
template<typename Matrix>
Matrix projector(const Matrix& B) {
  return B * linalg::pseudo_inverse(B);
}

/// Compute projector from the \param B matrix with constraint \param C
/// \tparam Matrix
/// \param B Wilson's B matrix
/// \param C Constraint matrix
/// \return Projector
template<typename Matrix>
Matrix projector(const Matrix& B, const Matrix& C) {
  // Standard projector
  Matrix P{B * linalg::pseudo_inverse(B)};

  // Projector with constraints
  return P - P * C * linalg::inv<Matrix>(C * P * C) * C * P;
}

} // namespace wilson

} // namespace irc

#endif // IRC_WILSON_H
