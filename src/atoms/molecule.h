#ifndef IRC_MOLECULE_H
#define IRC_MOLECULE_H

#include "atom.h"

#include <ostream>
#include <utility>
#include <vector>

namespace molecule{

/// Molecule as collection of atoms
template<typename Vector3>
using Molecule = std::vector<atom::Atom<Vector3>>;

template<typename Vector3>
Molecule<Vector3> make_molecule(
    const std::vector<std::tuple<size_t, double, double, double>>& atoms){
  size_t n_atoms{ atoms.size() };
  
  Molecule<Vector3> molecule;
  
  size_t an;
  double x, y, z;
  for(size_t i{0}; i < n_atoms; i++){
    std::tie(an, x, y, z) = atoms[i];
    
    // TODO: Preallocate memory
    molecule.push_back({an, {x,y,z}});
  }
  
  return molecule;
}

/// Compute the total mass of a molecule
///
/// \tparam T 3D vector
/// \param molecule Molecule (collection of atoms)
/// \return Mass of the molecule
template<typename Vector3>
double mass(const Molecule<Vector3>& molecule){
  double m{0};
  for(const auto& atom : molecule){
    m += atom::mass(atom.atomic_number);
  }
  
  return m;
}

// TODO: Define operator*?
/// Multiply all atomic positions within a molecule by a given \param multiplier
///
/// \tparam T 3D vector
/// \param molecule Molecule
/// \param multiplier Multiplier for atomic positions
template<typename T, typename Vector3>
void multiply_positions(Molecule<Vector3>& molecule, T multiplier){
  for(auto& atom : molecule){
    atom.position = atom.position * multiplier;
  }
}

/// Print a molecule
///
/// \tparam T 3D vector
/// \param out Output stream
/// \param molecule Molecule
/// \return Output stream
template<typename Vector3>
std::ostream& operator<<(std::ostream& out, const Molecule<Vector3>& molecule){
  for(const auto& atom : molecule){
    out << atom;
  }
  
  return out;
}

}

#endif //IRC_MOLECULE_H
