#include "libirc/mathtools.h"

namespace irc {

namespace tools {

namespace math {

double pirange_rad(double angle) {
  if (angle > constants::pi) {
    return pirange_rad(angle - 2. * constants::pi);
  } else if (angle <= -constants::pi) {
    return pirange_rad(angle + 2. * constants::pi);
  } else {
    return angle;
  }
}

double pirange_deg(double angle){
  return rad_to_deg( pirange_rad( deg_to_rad(angle) ) );
}

// TODO: Move to conversion
double deg_to_rad(double angle_deg){
  return angle_deg / 180. * tools::constants::pi;
}

// TODO: Move to conversion
double rad_to_deg(double angle_rad){
  return angle_rad * 180 / tools::constants::pi;
}

} // namespace math

} // namespace tools

} // namespace irc