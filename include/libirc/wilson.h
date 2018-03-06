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
Matrix
wilson_matrix(const Vector& x_cartesian,
              const std::vector<connectivity::Bond>& bonds,
              const std::vector<connectivity::Angle>& angles = {},
              const std::vector<connectivity::Dihedral>& dihedrals = {}) {
  // Get number of atoms
  const size_t n_atoms{linalg::size<Vector>(x_cartesian) / 3};

  // Get the total number of internal redundant coordinates
  const size_t n_irc{bonds.size() + angles.size() + dihedrals.size()};

  // Allocate Wilson's B matrix
  Matrix B{linalg::zeros<Matrix>(n_irc, 3 * n_atoms)};

  // Utility vectors for atomic positions
  Vector3 p1, p2, p3, p4;

  // Utility vector for gradients storage
  Vector3 g1, g2, g3, g4;

  // B-matrix rows offset
  size_t offset{0};

  // Populate B matrix's rows corresponding to bonds
  connectivity::Bond bond;
  for (size_t i{0}; i < bonds.size(); i++) {
    bond = bonds[i];

    for (size_t m{0}; m < 3; m++) {
      p1(m) = x_cartesian(3 * bond.i + m);
      p2(m) = x_cartesian(3 * bond.j + m);
    }

    std::tie(g1, g2) = bond_gradient(p1, p2);

    for (size_t idx{0}; idx < 3; idx++) {
      B(i, 3 * bond.i + idx) = g1(idx);
      B(i, 3 * bond.j + idx) = g2(idx);
    }
  }

  // Populate B matrix's rows corresponding to angles
  offset = bonds.size();
  connectivity::Angle angle;
  for (size_t i{0}; i < angles.size(); i++) {
    angle = angles[i];

    for (size_t m{0}; m < 3; m++) {
      p1(m) = x_cartesian(3 * angle.i + m);
      p2(m) = x_cartesian(3 * angle.j + m);
      p3(m) = x_cartesian(3 * angle.k + m);
    }

    std::tie(g1, g2, g3) = angle_gradient(p1, p2, p3);

    for (size_t idx{0}; idx < 3; idx++) {
      B(i + offset, 3 * angle.i + idx) = g1(idx);
      B(i + offset, 3 * angle.j + idx) = g2(idx);
      B(i + offset, 3 * angle.k + idx) = g3(idx);
    }
  }

  // Populate B matrix's rows corresponding to dihedrals
  offset = bonds.size() + angles.size();
  connectivity::Dihedral dihedral;
  for (size_t i{0}; i < dihedrals.size(); i++) {
    dihedral = dihedrals[i];

    for (size_t m{0}; m < 3; m++) {
      p1(m) = x_cartesian(3 * dihedral.i + m);
      p2(m) = x_cartesian(3 * dihedral.j + m);
      p3(m) = x_cartesian(3 * dihedral.k + m);
      p4(m) = x_cartesian(3 * dihedral.l + m);
    }

    std::tie(g1, g2, g3, g4) = dihedral_gradient(p1, p2, p3, p4);

    for (size_t idx{0}; idx < 3; idx++) {
      B(i + offset, 3 * dihedral.i + idx) = g1(idx);
      B(i + offset, 3 * dihedral.j + idx) = g2(idx);
      B(i + offset, 3 * dihedral.k + idx) = g3(idx);
      B(i + offset, 3 * dihedral.l + idx) = g4(idx);
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
    double dx = 1.e-6) {

  // Number of cartesian coordinates
  const size_t n_c{linalg::size(x_c)};

  // Number of IRC
  const size_t n_irc{bonds.size() + angles.size() + dihedrals.size()};

  // Allocate Wilson B matrix
  Matrix B{linalg::zeros<Matrix>(n_irc, n_c)};

  // Allocate displaced cartesian coordinates
  Vector x_c_pm{x_c};

  // Allocate displaced IRC
  Vector q_irc_plus{linalg::zeros<Vector>(n_irc)};
  Vector q_irc_minus{linalg::zeros<Vector>(n_irc)};

  for (size_t j{0}; j < n_c; j++) {
    // Compute positive displacement for cartesian coordinate j
    x_c_pm(j) += dx;

    // Compute IRC corresponding to positive displacement of x_c(j)
    q_irc_plus = connectivity::cartesian_to_irc<Vector3, Vector>(
        x_c_pm, bonds, angles, dihedrals);

    // Compute negative displacement for cartesian coordinate j
    x_c_pm(j) -= 2 * dx;

    // Compute IRC corresponding to negative displacement of x_c(j)
    q_irc_minus = connectivity::cartesian_to_irc<Vector3, Vector>(
        x_c_pm, bonds, angles, dihedrals);

    for (size_t i{0}; i < n_irc; i++) {
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
  // TODO: Pass iB instead of computing it
  return B * linalg::pseudo_inverse(B);
}

} // namespace wilson

} // namespace irc

#endif // IRC_WILSON_H
