#include "connectivity.h"

namespace connectivity{

constexpr double bond_multiplier{1.3};

arma::mat connectivity_matrix(const molecule::Molecule& molecule){
  size_t n_atoms{ molecule.size() };
  
  arma::mat connectivity{n_atoms, n_atoms, arma::fill::zeros};
  
  for(size_t j{0}; j < n_atoms; j++){
    for(size_t i{j+1}; i < n_atoms; i++){
      connectivity(i,j) = bond(molecule[i], molecule[j]);
      connectivity(j,i) = connectivity(i,j);
    }
  }
  
  return connectivity;
}

double bond(const atom::Atom& a1, const atom::Atom& a2){
  double r{ arma::norm(a1.position - a2.position)};
  double sum_covalent_radii{ atom::covalent_radius(a1.atomic_number) +
                             atom::covalent_radius(a2.atomic_number) };
  
  double bond{0.};
  
  if( r < sum_covalent_radii * bond_multiplier){
    bond = r;
  }
  
  return bond;
}

}
