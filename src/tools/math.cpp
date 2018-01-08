#include "math.h"

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
  if (angle > 180.) {
    return pirange_deg(angle - 360.);
  } else if (angle <= -180.) {
    return pirange_deg(angle + 360);
  } else {
    return angle;
  }
}

} // namespace math

} // namespace tools

} // namespace irc