#ifndef IRC_CONVERSION_H
#define IRC_CONVERSION_H

#include "constants.h"

namespace irc {

namespace tools {

namespace conversion {

/// Angstrom to bohr conversion factor
constexpr double angstrom_to_bohr{1.889725989};

/// Bohr to angstrom conversion factor
constexpr double bohr_to_angstrom{1. / angstrom_to_bohr};

/// Degrees to radians conversion factor
constexpr double deg_to_rad{ constants::pi / 180. };

/// Radians to degree conversion factor
constexpr double rad_to_deg{ 180. / constants::pi };

} // namespace irc

} // namespace tools

} // namespace conversion

#endif //IRC_CONVERSION_H
