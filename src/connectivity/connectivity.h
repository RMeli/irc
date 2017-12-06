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

constexpr double covalent_bond_multiplier{1.3};

template <typename Vector3>
struct Bond{
  size_t i;
  size_t j;
  Vector3 p1;
  Vector3 p2;
  double bond;
};

template <typename Vector3>
struct Angle{
  size_t i;
  size_t j;
  size_t k;
  Vector3 p1;
  Vector3 p2;
  Vector3 p3;
  double angle;
};

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

template <typename Vector3, typename Matrix>
Matrix connectivity_matrix(const molecule::Molecule<Vector3>& molecule){
  size_t n_atoms{ molecule.size() };
  
  Matrix connectivity{ linalg::zeros<Matrix>(n_atoms, n_atoms) };
  
  double r{0.};
  double sum_covalent_radii{0.};
  for(size_t i{0}; i < n_atoms; i++){
    for(size_t j{i+1}; j < n_atoms; j++){
      
      r = distance(molecule[i].position, molecule[j].position);
      
      sum_covalent_radii = atom::covalent_radius(molecule[i].atomic_number) +
                           atom::covalent_radius(molecule[j].atomic_number);
      
      if( r < sum_covalent_radii * covalent_bond_multiplier){
        connectivity(i,j) = r;
        connectivity(j,i) = r;
      }
    }
  }
  
  return std::move(connectivity);
}

template <typename Vector3, typename Matrix>
std::vector<Bond<Vector3>> bonds(const molecule::Molecule<Vector3>& molecule,
                                 const Matrix& connectivity,
                                 double epsilon = 1e-12){
  
  size_t n_atoms{ linalg::n_rows(connectivity) };
  
  std::vector<Bond<Vector3>> b;
  
  double d{0};
  for(size_t j{0}; j < n_atoms; j++){
    for(size_t i{j+1}; i < n_atoms; i++){
      d = connectivity(i,j);
      
      if( std::abs(d) > epsilon ){
        b.push_back(Bond<Vector3>{i, j,
                                  molecule[i].position,
                                  molecule[j].position, d});
      }
    }
  }
  
  return std::move(b);
}

template <typename Vector3, typename Matrix>
std::vector<Angle<Vector3>> angles(const molecule::Molecule<Vector3>& molecule,
                                   const Matrix& connectivity,
                                   double epsilon = 1e-12){
  size_t n_atoms{ linalg::n_rows(connectivity) };
  
  std::vector<Angle<Vector3>> ang;
  
  Vector3 p1{0., 0., 0.};
  Vector3 p2{0., 0., 0.};
  Vector3 p3{0., 0., 0.};
  
  double d1{0}, d2{0};
  for(size_t j{0}; j < n_atoms; j++){
    for(size_t i{j+1}; i < n_atoms; i++){
      d1 = connectivity(i,j);
      
      if( std::abs(d1) > epsilon ){
        p1 = molecule[i].position;
        p2 = molecule[j].position;
        
        for(size_t k{i+1}; k < n_atoms; k++){
          d2 = connectivity(j,k);
          if( std::abs(d2) > epsilon ){
            p3 = molecule[k].position;
            
            ang.push_back(Angle<Vector3>{i, j, k,
                                molecule[i].position,
                                molecule[j].position,
                                molecule[k].position,
                                angle(p1,p2,p3)});
          }
        }
      }
    }
  }
  
  return std::move(ang);
}

}

#endif //IRC_CONNECTIVITY_H
