#include "connectivity.h"

#include "../atoms/periodic_table.h"

#include <stdexcept>

namespace connectivity{

arma::vec get_bonds(std::vector<std::string> atoms_symbols,
                    const arma::mat& atoms_positions){
  
  if(atoms_symbols.size() != atoms_positions.n_rows){
    throw std::length_error("");
  }
  
  if(atoms_positions.n_cols != 3){
    throw std::length_error("");
  }
  
  size_t n_atoms{atoms_positions.n_rows};
  
  arma::vec bonds;
  
  size_t i_idx{0}, j_idx{0};
  double r{0};
  double i_cr{0}, j_cr{0};
  for(size_t i{0}; i < n_atoms; i++){
    i_idx = atoms::atomic_number(atoms_symbols[i]);
    i_cr = atoms::pt_covalent_radii[i_idx];
    
    for(size_t j{i+1}; j < n_atoms; j++){
      j_idx = atoms::atomic_number(atoms_symbols[j]);
      j_cr = atoms::pt_covalent_radii[j_idx];
      
      r = arma::norm( atoms_positions.row(i) - atoms_positions.row(j) );
      
      if( r < 1.3 * (i_cr + j_cr)){
        bonds = arma::join_cols(bonds, arma::vec{r});
      }
    }
  }
  
  return bonds;
}

}
