#ifndef IRC_CONNECTIVITY_H
#define IRC_CONNECTIVITY_H

#include "../atoms/atom.h"
#include "../atoms/molecule.h"
#include "../linear_algebra/linalg.h"
#include "../tools/constants.h"

#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace connectivity {

constexpr double bond_multiplier{1.3};

template <typename Vector3>
double distance(const Vector3& v1, const Vector3& v2){
  return linalg::norm(v1 - v2);
}

template <typename Vector3>
double angle(const Vector3& v1, const Vector3& v2, const Vector3& v3){
  Vector3 r1{v1 - v2};
  Vector3 r2{v3 - v2};
  
  double N{ linalg::norm(r1) * linalg::norm(r2) };
  
  double angle = std::acos( linalg::dot(r1, r2) / N );
  
  return angle * 180.0 /  tools::constants::pi;
}

template <typename Vector3>
double bond(const atom::Atom<Vector3> &a1, const atom::Atom<Vector3>& a2){
  double r{ distance(a1.position, a2.position) };
  
  double sum_covalent_radii{ atom::covalent_radius(a1.atomic_number) +
                             atom::covalent_radius(a2.atomic_number) };
  
  double bond{0.};
  
  if( r < sum_covalent_radii * bond_multiplier){
    bond = r;
  }
  
  return bond;
}

template <typename Vector3, typename Matrix>
Matrix connectivity_matrix(const molecule::Molecule<Vector3>& molecule){
  size_t n_atoms{ molecule.size() };
  
  Matrix connectivity{ linalg::zeros<Matrix>(n_atoms, n_atoms) };
  
  for(size_t i{0}; i < n_atoms; i++){
    for(size_t j{i+1}; j < n_atoms; j++){
      connectivity(i,j) = bond(molecule[i], molecule[j]);
      connectivity(j,i) = connectivity(i,j);
    }
  }
  
  return std::move(connectivity);
}

template <typename Matrix>
std::vector<double> bonds(const Matrix& connectivity, double epsilon = 1e-12){
  size_t n_atoms{ linalg::n_rows(connectivity) };
  
  std::vector<double> b;
  
  double d{0};
  for(size_t j{0}; j < n_atoms; j++){
    for(size_t i{j+1}; i < n_atoms; i++){
      d = connectivity(j,i);
      
      if( std::abs(d) > epsilon ){
        b.push_back(d);
      }
    }
  }
  
  return std::move(b);
}

template <typename Vector3, typename Matrix>
std::vector<double> angles(const molecule::Molecule<Vector3>& molecule,
                           const Matrix& connectivity){
  size_t n_atoms{ linalg::n_rows(connectivity) };
  
  std::vector<double> ang;
  
  Vector3 p1{0., 0., 0.};
  Vector3 p2{0., 0., 0.};
  Vector3 p3{0., 0., 0.};
  
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
  
  return std::move(ang);
}

}


#endif //IRC_CONNECTIVITY_H
