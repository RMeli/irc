#ifndef IRC_MATH_H
#define IRC_MATH_H

#include "libirc/constants.h"
#include "libirc/conversion.h"
#include "libirc/linalg.h"

namespace irc {

namespace tools {

namespace math {

/// Returns \param angle in the range \f$(-\pi,\pi]\f$
///
/// \param angle
/// \return
inline double pirange_rad(double angle) noexcept {
  if (angle > constants::pi) {
    return pirange_rad(angle - 2. * constants::pi);
  } else if (angle <= -constants::pi) {
    return pirange_rad(angle + 2. * constants::pi);
  } else {
    return angle;
  }
}

/// Returns \param angle in the range \f$(-180,180]\f$
///
/// \param angle
/// \return
inline double pirange_deg(double angle) noexcept {
  return pirange_rad(angle * conversion::deg_to_rad) * conversion::rad_to_deg;
}

/*! Check if two vectors @param v1 and @param v2 are collinear
 *
 * @tparam Vector3
 * @param v1 Vector 1
 * @param v2 Vector 2
 * @param tolerance
 * @return True if @param v1 and @param v2 are collinear, false othwerwise
 */
template<typename Vector3>
inline bool collinear(Vector3 v1, Vector3 v2, double tolerance = 1e-6) {
  const double l1{linalg::norm(v1)};
  const double l2{linalg::norm(v2)};

  const double angle{std::acos(linalg::dot(v1 / l1, v2 / l2))};

  bool c{false};

  if (std::abs(angle) < tolerance) {
    c = true;
  } else if (std::abs(angle - tools::constants::pi) < tolerance) {
    c = true;
  }

  return c;
}

} // namespace math

} // namespace tools

} // namespace irc

#endif // IRC_MATH_H
