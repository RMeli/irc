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
  size_t atomic_number;

  AtomicNumber(size_t an);

  AtomicNumber(const std::string& symbol);
};

AtomicNumber::AtomicNumber(size_t an) {
  if (!periodic_table::valid_atomic_number(an)) {
    throw std::logic_error("Invalid atomic number.");
  }

  atomic_number = an;
}

AtomicNumber::AtomicNumber(const std::string& symbol)
  : AtomicNumber(periodic_table::atomic_number(symbol)) {}

/// Get an atomic symbol from \class AtomicNumber
///
/// \param an \class AtomicNumber
/// \return Atomic symbol corresponding to \class AtomicNumber
std::string symbol(const AtomicNumber& an) noexcept {
  return periodic_table::symbols[an.atomic_number];
}

/// Get atomic symbol from \class AtomicNumber
///
/// \param an \class AtomicNumber
/// \return Mass corresponding to \class AtomicNumber
double mass(const AtomicNumber& an) noexcept {
  return periodic_table::masses[an.atomic_number];
}

/// Get covalent radius from \class AtomicNumber
///
/// \param an \class AtomicNumber
/// \return Covalent radius corresponding to \class AtomicNumber
double covalent_radius(const AtomicNumber& an) noexcept {
  return periodic_table::covalent_radii[an.atomic_number];
}

/// Get Van der Waals radius from \class AtomicNumber
///
/// \param an \class AtomicNumber
/// \return Van der Waals radius corresponding to \class AtomicNumber
double vdw_radius(const AtomicNumber& an) noexcept {
  return periodic_table::vdw_radii[an.atomic_number];
}

// TODO: constexpr for C++14
/// Check if an atom is either N, O, F, P, S or Cl
///
/// \param an \class AtomicNumber
/// \return
bool is_NOFPSCl(const AtomicNumber& an) noexcept {
  const size_t n{an.atomic_number};

  return (n == 7 or n == 8 or n == 9 or n == 15 or n == 16 or n == 17);
}

/// Check if an atom is a hydrogen atom
///
/// \param an \class AtomicNumber
/// \return
constexpr bool is_H(const AtomicNumber& an) noexcept {
  return an.atomic_number == 1;
}

/// Print \class AtomicNumber
///
/// \param out Output stream
/// \param an \class AtomicNumber
/// \return Output stream
std::ostream& operator<<(std::ostream& out, const AtomicNumber& an);

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

template<typename Vector3>
Atom<Vector3>::Atom(const AtomicNumber& an, const Vector3& pos)
  : atomic_number(an), position(pos) {}

template<typename Vector3>
Atom<Vector3>::Atom(const std::string& symbol, const Vector3& pos)
  : atomic_number(symbol), position(pos) {}

/// Output operator for an atom
///
/// \tparam T Mathematical vector (possibly of three dimensions)
/// \param out Output stream
/// \param a \class Atom<T>
/// \return Output stream
template<typename Vector3>
std::ostream& operator<<(std::ostream& out, const Atom<Vector3>& a) {
  // Print top line
  out << std::left << std::setw(15) << std::setfill('-') << '+';
  out << ' ';
  out << std::left << std::setw(2) << std::setfill(' ');
  out << symbol(a.atomic_number) << ' ';
  out << std::right << std::setw(15) << std::setfill('-') << '+' << std::endl;

  // Print atomic number
  out << std::setfill(' ');
  out << std::left << std::setw(20) << "| atomic number:";
  out << std::right << std::setw(12) << a.atomic_number;
  out << " |" << std::endl;

  // Print atomic mass
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "| atomic mass:";
  out << std::right << std::setw(12);
  out << mass(a.atomic_number) << " |" << std::endl;

  // Print covalent radius
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "| covalent radius:";
  out << std::right << std::setw(12);
  out << covalent_radius(a.atomic_number) << " |" << std::endl;

  // Print position
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "| position:";
  out << std::right << std::setw(14);
  out << " |" << std::endl;

  // x coordinate
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "|    x =";
  out << std::right << std::setw(12) << a.position(0);
  out << " |" << std::endl;

  // y coordinate
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "|    y =";
  out << std::right << std::setw(12) << a.position(1);
  out << " |" << std::endl;

  // z coordinate
  out << std::setfill(' ');
  out << std::fixed << std::setprecision(5) << std::scientific;
  out << std::left << std::setw(20) << "|    x =";
  out << std::right << std::setw(12) << a.position(2);
  out << " |" << std::endl;

  // Print bottom line
  out << std::left << std::setw(17) << std::setfill('-') << '+';
  out << std::right << std::setw(17) << std::setfill('-') << '+' << std::endl;

  return out;
}

} // namespace atom

} // namespace irc

#endif // IRC_ATOM_H
