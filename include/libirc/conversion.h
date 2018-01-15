#ifndef IRC_CONVERSION_H
#define IRC_CONVERSION_H

namespace irc {

namespace tools {

namespace conversion {

/// Angstrom to bohr conversion factor
constexpr double angstrom_to_bohr{1.889725989};

/// Bohr to angstrom conversion factor
constexpr double bohr_to_angstrom{1. / angstrom_to_bohr};

} // namespace irc

} // namespace tools

} // namespace conversion

#endif //IRC_CONVERSION_H
