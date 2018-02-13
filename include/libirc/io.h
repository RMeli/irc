#ifndef IRC_IO_H
#define IRC_IO_H

#include "connectivity.h"
#include "molecule.h"

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <vector>

namespace irc {

namespace io {

/// Load molecule in XYZ format from input stream
///
/// \tparam Vector3 3D vector
/// \param in Input stream
/// \return Molecule
template<typename Vector3>
molecule::Molecule<Vector3> load_xyz(std::istream &in) {
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
molecule::Molecule<Vector3> load_xyz(std::string fname) {
  std::ifstream in{fname};

  if (!in.is_open()) {
    throw std::runtime_error("Impossible to open file " + fname);
  }

  return load_xyz<Vector3>(in);
}

/// Print bonds
template<typename Vector3, typename Vector>
void print_bonds(const Vector &x_c,
                 const std::vector<connectivity::Bond> &bonds,
                 std::ostream &out = std::cout) {

  // Total number of bonds
  const size_t n_bonds{bonds.size()};

  out << n_bonds << " bonds (\u212B):" << std::endl;

  // Atomic positions
  Vector3 p1{0., 0., 0.}, p2{0., 0., 0.};

  // Indices
  size_t idx_i{0}, idx_j{0};

  out.precision(3);
  out.fill(' ');

  for (size_t i{0}; i < n_bonds; i++) {
    // Get indices
    idx_i = bonds[i].i;
    idx_j = bonds[i].j;

    // Print indices
    out.width(1);
    out << std::left << '(';
    out.width(4);
    out << std::right << idx_i << ',';
    out.width(4);
    out << std::right << idx_j << ')';

    // Get positions
    p1 = {x_c(3 * idx_i), x_c(3 * idx_i + 1), x_c(3 * idx_i + 2)};
    p2 = {x_c(3 * idx_j), x_c(3 * idx_j + 1), x_c(3 * idx_j + 2)};

    // Print distance

    out << std::setw(14) << std::fixed << " " << std::right
        << connectivity::distance(p1, p2) * tools::conversion::bohr_to_angstrom
        << std::endl;
  }
}

template<typename Vector3, typename Vector>
void print_angles(const Vector &x_c,
                  const std::vector<connectivity::Angle> &angles,
                  std::ostream &out = std::cout) {

  // Total number of angles
  const size_t n_angles{angles.size()};

  out << n_angles << " angles (\u00B0):" << std::endl;

  // Atomic positions
  Vector3 p1{0., 0., 0.}, p2{0., 0., 0.}, p3{0., 0., 0.};

  // Indices
  size_t idx_i{0}, idx_j{0}, idx_k{0};

  out.precision(2);
  out.fill(' ');

  for (size_t i{0}; i < n_angles; i++) {
    // Get indices
    idx_i = angles[i].i;
    idx_j = angles[i].j;
    idx_k = angles[i].k;

    // Print indices
    out.width(1);
    out << std::left << '(';
    out.width(4);
    out << std::right << idx_i << ',';
    out.width(4);
    out << std::right << idx_j << ',';
    out.width(4);
    out << std::right << idx_k << ')';

    // Get positions
    p1 = {x_c(3 * idx_i), x_c(3 * idx_i + 1), x_c(3 * idx_i + 2)};
    p2 = {x_c(3 * idx_j), x_c(3 * idx_j + 1), x_c(3 * idx_j + 2)};
    p3 = {x_c(3 * idx_k), x_c(3 * idx_k + 1), x_c(3 * idx_k + 2)};

    // Print distance
    out << std::setw(14) << std::fixed << std::right
        << connectivity::angle(p1, p2, p3) * tools::conversion::rad_to_deg
        << std::endl;
  }
}

template<typename Vector3, typename Vector>
void print_dihedrals(const Vector &x_c,
                     const std::vector<connectivity::Dihedral> &dihedrals,
                     std::ostream &out = std::cout) {

  // Total number of angles
  const size_t n_dihedrals{dihedrals.size()};

  out << n_dihedrals << " dihedrals (\u00B0):" << std::endl;

  // Atomic positions
  Vector3 p1{0., 0., 0.}, p2{0., 0., 0.}, p3{0., 0., 0.}, p4{0., 0., 0.};

  // Indices
  size_t idx_i{0}, idx_j{0}, idx_k{0}, idx_l{0};

  out.precision(2);
  out.fill(' ');

  for (size_t i{0}; i < n_dihedrals; i++) {
    // Get indices
    idx_i = dihedrals[i].i;
    idx_j = dihedrals[i].j;
    idx_k = dihedrals[i].k;
    idx_l = dihedrals[i].l;

    // Print indices
    out.width(1);
    out << std::left << '(';
    out.width(4);
    out << std::right << idx_i << ',';
    out.width(4);
    out << std::right << idx_j << ',';
    out.width(4);
    out << std::right << idx_k << ',';
    out.width(4);
    out << std::right << idx_l << ')';

    // Get positions
    p1 = {x_c(3 * idx_i), x_c(3 * idx_i + 1), x_c(3 * idx_i + 2)};
    p2 = {x_c(3 * idx_j), x_c(3 * idx_j + 1), x_c(3 * idx_j + 2)};
    p3 = {x_c(3 * idx_k), x_c(3 * idx_k + 1), x_c(3 * idx_k + 2)};
    p4 = {x_c(3 * idx_l), x_c(3 * idx_l + 1), x_c(3 * idx_l + 2)};

    // Print distance

    out << std::setw(9) << std::fixed << std::right
        << connectivity::dihedral(p1, p2, p3, p4) *
               tools::conversion::rad_to_deg
        << std::endl;
  }
}

} // namespace io

} // namespace irc

#endif // IRC_IO_H
