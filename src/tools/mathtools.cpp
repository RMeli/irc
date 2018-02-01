#include "libirc/mathtools.h"
#include "libirc/conversion.h"

namespace irc {

namespace tools {

namespace math {

double pirange_rad(double angle) noexcept {
  if (angle > constants::pi) {
    return pirange_rad(angle - 2. * constants::pi);
  } else if (angle <= -constants::pi) {
    return pirange_rad(angle + 2. * constants::pi);
  } else {
    return angle;
  }
}

double pirange_deg(double angle) noexcept {
  return pirange_rad(angle * conversion::deg_to_rad) * conversion::rad_to_deg;
}

} // namespace math

} // namespace tools

} // namespace irc