#ifndef IRC_MOLECULE_H
#define IRC_MOLECULE_H

#include "atom.h"
#include "linalg.h"

#include <ostream>
#include <utility>
#include <vector>

namespace irc {

namespace molecule {

/// Molecule as collection of atoms
template<typename Vector3>
using Molecule = std::vector<atom::Atom<Vector3>>;

/// Compute the total mass of a molecule
///
/// \tparam T 3D vector
/// \param molecule Molecule (collection of atoms)
/// \return Mass of the molecule
template<typename Vector3>
double mass(const Molecule<Vector3> &molecule) noexcept {
  double m{0};
  for (const auto &atom : molecule) {
    m += atom::mass(atom.atomic_number);
  }

  return m;
}

// TODO: Avoid modification to the molecule and define operator*
/// Multiply all atomic positions within a molecule by a given \param multiplier
///
/// \tparam T 3D vector
/// \param molecule Molecule
/// \param multiplier Multiplier for atomic positions
template<typename T, typename Vector3>
void multiply_positions(Molecule<Vector3> &molecule, T multiplier) {
  for (auto &atom : molecule) {
    atom.position = atom.position * multiplier;
  }
}

/// Get cartesian coordinates of all atoms in \param molecule
///
/// \tparam Vector3
/// \tparam Vector
/// \param molecule Molecule
/// \return Cartesian coordinates of all atoms in \param molecule
///
/// The cartesian coordinates are stored in a linear vector so that the first
/// three entries are the (x,y,z) coordinates of the first atom and so on.
template<typename Vector3, typename Vector>
Vector to_cartesian(const Molecule<Vector3> &molecule) {
  size_t n_atoms{molecule.size()};

  Vector x_cartesian{linalg::zeros<Vector>(3 * n_atoms)};

  for (size_t i{0}; i < n_atoms; i++) {
    for (size_t idx{0}; idx < 3; idx++) {
      x_cartesian(3 * i + idx) = molecule[i].position(idx);
    }
  }

  return x_cartesian;
}

/// Print a molecule
///
/// \tparam T 3D vector
/// \param out Output stream
/// \param molecule Molecule
/// \return Output stream
template<typename Vector3>
std::ostream &operator<<(std::ostream &out, const Molecule<Vector3> &molecule) {
  for (const auto &atom : molecule) {
    out << atom;
  }

  return out;
}

} // namespace molecule

} // namespace irc

#endif // IRC_MOLECULE_H
