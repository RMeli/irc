#ifndef IRC_MOLECULE_H
#define IRC_MOLECULE_H

#include "atom.h"

#include <vector>

namespace molecule{

using Molecule = std::vector<atom::Atom>;

double mass(const Molecule& molecule);

}


#endif //IRC_MOLECULE_H
