#ifndef IRC_ATOM_H
#define IRC_ATOM_H

#include "periodic_table.h"

#include <iomanip>
#include <ostream>
#include <string>

namespace irc {

namespace atom {

/// Atomic number
///
/// The \class AtomicNumber represents a valid atomic number, for which
/// all the quantities accessible with dedicated functions are available
/// in \headerfile periodic_table.h
struct AtomicNumber {
  std::size_t atomic_number;

  AtomicNumber(std::size_t an);

  AtomicNumber(const std::string& symbol);
};

inline AtomicNumber::AtomicNumber(std::size_t an) {
  if (!periodic_table::valid_atomic_number(an)) {
    throw std::logic_error("Invalid atomic number.");
  }

  atomic_number = an;
}

inline AtomicNumber::AtomicNumber(const std::string& symbol)
  : AtomicNumber(periodic_table::atomic_number(symbol)) {}

/// Get an atomic symbol from \class AtomicNumber
///
/// \param an \class AtomicNumber
/// \return Atomic symbol corresponding to \class AtomicNumber
inline std::string symbol(const AtomicNumber& an) noexcept {
  return periodic_table::symbols[an.atomic_number];
}

// TODO: constexpr for C++14
/// Get atomic symbol from \class AtomicNumber
///
/// \param an \class AtomicNumber
/// \return Mass corresponding to \class AtomicNumber
inline double mass(const AtomicNumber& an) noexcept {
  return periodic_table::masses[an.atomic_number];
}

// TODO: constexpr for C++14
/// Get covalent radius from \class AtomicNumber
///
/// \param an \class AtomicNumber
/// \return Covalent radius corresponding to \class AtomicNumber
inline double covalent_radius(const AtomicNumber& an) noexcept {
  return periodic_table::covalent_radii[an.atomic_number];
}

// TODO: constexpr for C++14
/// Get Van der Waals radius from \class AtomicNumber
///
/// \param an \class AtomicNumber
/// \return Van der Waals radius corresponding to \class AtomicNumber
inline double vdw_radius(const AtomicNumber& an) noexcept {
  return periodic_table::vdw_radii[an.atomic_number];
}

// TODO: constexpr for C++14
/// Check if an atom is either N, O, F, P, S or Cl
///
/// \param an \class AtomicNumber
/// \return
inline bool is_NOFPSCl(const AtomicNumber& an) noexcept {
  const std::size_t n{an.atomic_number};

  return (n == 7 or n == 8 or n == 9 or n == 15 or n == 16 or n == 17);
}

/// Check if an atom is a hydrogen atom
///
/// \param an \class AtomicNumber
/// \return
inline constexpr bool is_H(const AtomicNumber& an) noexcept {
  return an.atomic_number == 1;
}

/// Class representing an atom
///
/// \tparam Vector3 Mathematical vector (possibly of three dimensions)
///
/// An atom is defined only by its atomic number and its position in space.
/// Its properties can be looked up in \headerfile periodic_table.h
template<typename Vector3>
struct Atom {
  /// Atomic number
  AtomicNumber atomic_number;

  /// Position (in 3D space)
  Vector3 position;

  /// Constructor from \class AtomicNumber
  ///
  /// \param an \class AtomicNumber
  /// \param pos Position (in 3D space)
  Atom(const AtomicNumber& an, const Vector3& pos = {0., 0., 0.});

  /// Constructor from atomic symbol
  ///
  /// \param symbol Atomic symbol
  /// \param pos Position (in 3D space)
  Atom(const std::string& symbol, const Vector3& pos = {0., 0., 0.});
};

/*! Construct and atom from atomic number @param an
 *
 * @tparam Vector3
 * @param an Atomic number
 * @param pos Atomic position
 */
template<typename Vector3>
Atom<Vector3>::Atom(const AtomicNumber& an, const Vector3& pos)
  : atomic_number(an), position(pos) {}

/*! Construct and atom from atomic symbol @param symbol
 *
 * @tparam Vector3
 * @param symbol Atomic symbol
 * @param pos Atomic position
 */
template<typename Vector3>
Atom<Vector3>::Atom(const std::string& symbol, const Vector3& pos)
  : atomic_number(symbol), position(pos) {}

} // namespace atom

} // namespace irc

#endif // IRC_ATOM_H
