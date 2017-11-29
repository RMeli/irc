#ifndef IRC_CONNECTIVITY_H
#define IRC_CONNECTIVITY_H

#include <string>
#include <vector>

#include <armadillo>

namespace connectivity{

  arma::vec get_bonds(std::vector<std::string> atoms_symbols,
                      const arma::mat& atoms_positions);
}


#endif //IRC_CONNECTIVITY_H
