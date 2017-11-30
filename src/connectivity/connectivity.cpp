#include "connectivity.h"

#include "../tools/comparison.h"
#include "../tools/constants.h"

#include <cmath>

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

double distance(const arma::vec &v1, const arma::vec& v2){
  return arma::norm(v1 - v2);
}


double angle(const arma::vec& v1, const arma::vec& v2, const arma::vec& v3){
  arma::vec r1{v1 - v2};
  arma::vec r2{v3 - v2};
  
  double angle = std::acos( arma::norm_dot(r1, r2) );
  
  return angle * 180.0 /  tools::constants::pi;
}


double bond(const atom::Atom& a1, const atom::Atom& a2){
  double r{ distance(a1.position, a2.position) };
  double sum_covalent_radii{ atom::covalent_radius(a1.atomic_number) +
                             atom::covalent_radius(a2.atomic_number) };
  
  double bond{0.};
  
  if( r < sum_covalent_radii * bond_multiplier){
    bond = r;
  }
  
  return bond;
}

std::vector<double> bonds(const arma::mat& connectivity){
  size_t n_atoms{ connectivity.n_rows };
  
  std::vector<double> b;
  
  double d{0};
  for(size_t j{0}; j < n_atoms; j++){
    for(size_t i{j+1}; i < n_atoms; i++){
      d = connectivity(i,j);
      
      if( !tools::comparison::nearly_zero(d) ){
        b.push_back(d);
      }
    }
  }
  
  return b;
}

std::vector<double> angles(const molecule::Molecule& molecule,
                           const arma::mat& connectivity){
  size_t n_atoms{ connectivity.n_rows };
  
  std::vector<double> ang;
  
  arma::vec p1{3, arma::fill::zeros};
  arma::vec p2{3, arma::fill::zeros};
  arma::vec p3{3, arma::fill::zeros};
  
  double d1{0}, d2{0};
  for(size_t j{0}; j < n_atoms; j++){
    for(size_t i{j+1}; i < n_atoms; i++){
      d1 = connectivity(i,j);
      
      if( !tools::comparison::nearly_zero(d1) ){
        p1 = molecule[i].position;
        p2 = molecule[j].position;
        
        for(size_t k{i+1}; k < n_atoms; k++){
          d2 = connectivity(j,k);
          if( !tools::comparison::nearly_zero(d2) ){
            p3 = molecule[k].position;
            
            ang.push_back( angle(p1,p2,p3) );
          }
        }
      }
    }
  }
  
  return ang;
}

}
