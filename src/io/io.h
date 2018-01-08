#ifndef IRC_IO_H
#define IRC_IO_H

#include "atoms/molecule.h"

#include <istream>
#include <stdexcept>

namespace irc {

namespace io {

/// Load molecule in XYZ format from input stream
///
/// \tparam Vector3 3D vector
/// \param in Input stream
/// \return Molecule
template<typename Vector3>
molecule::Molecule <Vector3> load_xyz(std::istream &in) {
  size_t n_atoms{0};
  std::string dummy{""};
  
  // Read header
  in >> n_atoms;
  std::getline(in, dummy);
  
  std::string atom{""};
  double x{0.}, y{0.}, z{0.};
  
  molecule::Molecule<Vector3> molecule;
  while (in >> atom >> x >> y >> z) {
    molecule.push_back(atom::Atom<Vector3>{atom, {x, y, z}});
  }
  
  return molecule;
}

/// Load molecule in XYZ format from file
///
/// \tparam Vector3 3D vector
/// \param in Input stream
/// \return Molecule
template<typename Vector3>
molecule::Molecule <Vector3> load_xyz(std::string fname) {
  std::ifstream in(fname);
  
  if (!in.is_open()) {
    throw std::runtime_error("Impossible to open file " + fname);
  }
  
  return load_xyz<Vector3>(in);
}

} // namespace io

} // namespace irc

#endif //IRC_IO_H
