#ifndef IRC_CONNECTIVITY_H
#define IRC_CONNECTIVITY_H

#include "../atoms/atom.h"
#include "../atoms/molecule.h"

#include <string>
#include <vector>

#include <armadillo>

namespace connectivity {

arma::mat connectivity_matrix(const molecule::Molecule &molecule);

double bond(const atom::Atom &a1, const atom::Atom &a2);

}


#endif //IRC_CONNECTIVITY_H
