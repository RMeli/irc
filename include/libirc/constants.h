#ifndef IRC_CONSTANTS_H
#define IRC_CONSTANTS_H

#include <cmath>

namespace irc {

namespace tools {

namespace constants {

/// \f$\pi\f$
constexpr double pi{M_PI};

/// Covalent bond multiplier
constexpr double covalent_bond_multiplier{1.3};

/// Van der Waals bond multiplier
constexpr double vdw_bond_multiplier{0.9};

} // namespace constants

} // namespace tools

} // namespace irc

#endif //IRC_CONSTANTS_H
