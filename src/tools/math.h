#ifndef IRC_MATH_H
#define IRC_MATH_H

#include "constants.h"

namespace irc {

namespace tools {

namespace math {

/// Returns \param angle in the range \f$(-\pi,\pi]\f$
///
/// \param angle
/// \return
double pirange(double angle);

} // namespace math

} // namespace tools

} // namespace irc

#endif //IRC_MATH_H
