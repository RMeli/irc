#ifndef IRC_CONNECTIVITY_H
#define IRC_CONNECTIVITY_H

#include "../atoms/atom.h"
#include "../atoms/molecule.h"

#include <string>
#include <vector>

#include <armadillo>

namespace connectivity {

arma::mat connectivity_matrix(const molecule::Molecule &molecule);

double distance(const arma::vec &v1, const arma::vec& v2);

double angle(const arma::vec &v1, const arma::vec& v2, const arma::vec& v3);

double bond(const atom::Atom &a1, const atom::Atom& a2);

std::vector<double> bonds(const arma::mat& connectivity);

std::vector<double> angles(const molecule::Molecule& molecule,
                           const arma::mat& connectivity);

}


#endif //IRC_CONNECTIVITY_H
