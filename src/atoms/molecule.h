#ifndef IRC_MOLECULE_H
#define IRC_MOLECULE_H

#include "atom.h"

#include <ostream>
#include <vector>

namespace molecule{

using Molecule = std::vector<atom::Atom>;

double mass(const Molecule& molecule);

void multiply_positions(Molecule& molecule, double multiplier);

std::ostream& operator<<(std::ostream& out, const Molecule& molecule);

}


#endif //IRC_MOLECULE_H
