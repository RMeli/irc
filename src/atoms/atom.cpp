#include "libirc/atom.h"

#include <stdexcept>

namespace irc {

namespace atom {

AtomicNumber::AtomicNumber(size_t an) {
  if (!periodic_table::valid_atomic_number(an)) {
    throw std::logic_error("Invalid atomic number.");
  }

  atomic_number = an;
}

AtomicNumber::AtomicNumber(const std::string &symbol)
  : AtomicNumber(periodic_table::atomic_number(symbol)) {}

std::string symbol(const AtomicNumber &an) noexcept {
  return periodic_table::pt_symbols[an.atomic_number];
}

double mass(const AtomicNumber &an) noexcept {
  return periodic_table::pt_masses[an.atomic_number];
}

double covalent_radius(const AtomicNumber &an) noexcept {
  return periodic_table::pt_covalent_radii[an.atomic_number];
}

double vdw_radius(const AtomicNumber &an) noexcept {
  return periodic_table::pt_vdv_radii[an.atomic_number];
}

// TODO: constexpr for C++14
bool is_NOFPSCl(const AtomicNumber &an) noexcept {
  size_t n{an.atomic_number};

  return (n == 7 or n == 8 or n == 9 or n == 15 or n == 16 or n == 17);
}

bool is_H(const AtomicNumber &an) noexcept {
  return an.atomic_number == 1;
}

std::ostream &operator<<(std::ostream &out, const AtomicNumber &an) {
  out << an.atomic_number;

  return out;
}

} // namespace atom

} // namespace irc
