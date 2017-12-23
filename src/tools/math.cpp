#include "math.h"

namespace irc {

namespace tools {

namespace math {

double pirange(double angle) {
  
  if (angle > constants::pi) {
    return pirange(angle - 2. * constants::pi);
  } else if (angle <= -constants::pi) {
    return pirange(angle + 2. * constants::pi);
  } else {
    return angle;
  }
}

} // namespace math

} // namespace tools

} // namespace irc