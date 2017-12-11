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

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

namespace connectivity {

constexpr double covalent_bond_multiplier{1.3};

using EdgeProperty = boost::property<boost::edge_weight_t, double>;

using UGraph =
  boost::adjacency_list<
      boost::vecS, //
      boost::vecS, //
      boost::undirectedS, // Graph type
      boost::no_property, // Vertex property
      EdgeProperty // Edge property
  >;

using Vertex = boost::graph_traits<UGraph>::vertex_descriptor;
using Edge = boost::graph_traits<UGraph>::edge_descriptor;

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

/// Compute the distance matrix for \param molecule
///
/// \tparam Vector3 3D vector
/// \tparam Matrix Matrix
/// \param molecule Molecule (collection of atoms)
/// \return Distance matrix
///
/// The distance matrix \f$\mathbf{D}\f$ for \param molecule is given by
/// \f[
///   D_{ij} = |\mathbf{r}_i - \mathbf{r}_j|,
/// \f]
/// where the matrix element \f$D_{ij}\f$ is the distance between atom at
/// position \f$\mathbf{r}_i\f$ and the atom at position \f$\mathbf{r}_j\f$.
template <typename Vector3, typename Matrix>
Matrix distance_matrix(const molecule::Molecule<Vector3>& molecule){
  size_t n_atoms{ molecule.size() };
  
  Matrix distance_m{ linalg::zeros<Matrix>(n_atoms, n_atoms) };
  
  double r{0.};
  for(size_t i{0}; i < n_atoms; i++){
    for(size_t j{0}; j < n_atoms; j++){
      
      r = distance(molecule[i].position, molecule[j].position);
      
      distance_m(i,j) = r;
      distance_m(j,i) = r;
    }
  }
  
  return distance_m;
}

/// Compute adjacency matrix for \param molecule
///
/// \tparam Vector3 3D vector
/// \tparam Matrix Matrix
/// \param distance_m Distance matrix for \param molecule
/// \param molecule Molecule
/// \return Adjacency matrix
///
/// The adjacency matrix is represented here by a boost::adjacency_list object
/// as implemented in the Boost Graph Library (BGL).
/// The number of vertices corresponds to the number of atoms, while the
/// number of edges is determined by covalent bonding.
/// Hydrogen bonding is taken into account separately.
template <typename Vector3, typename Matrix>
UGraph adjacency_matrix(const Matrix& distance_m,
                        const molecule::Molecule<Vector3>& molecule){
  
  // Extract number of atoms
  size_t n_atoms{ molecule.size() };
  
  // Define a undirected graph with n_atoms vertices
  UGraph ug(n_atoms);
  
  double d{0.};
  double sum_covalent_radii{0.};
  for(size_t j{0}; j < n_atoms; j++){
    for(size_t i{j+1}; i < n_atoms; i++){
  
      // Extract distance between atom i and atom j
      d = distance_m(i,j);
      
      // Compute sum of covalent radii for atoms i and j
      sum_covalent_radii = atom::covalent_radius(molecule[i].atomic_number) +
                           atom::covalent_radius(molecule[j].atomic_number);
      
      // Determine if atoms i and j are bonded
      if( d < covalent_bond_multiplier * sum_covalent_radii ){
        // Add edge to boost::adjacency_list between vertices i and j
        // Store the distance d between atoms i and j as the edge weight
        boost::add_edge(i, j, d, ug);
      }
    }
  }
  
  return ug;
}

/*
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
*/

template <typename Vector3>
std::vector<Bond<Vector3>> bonds(const UGraph& ug,
                                 const molecule::Molecule<Vector3>& molecule){
  
  size_t n_atoms{ molecule.size() };
  
  std::vector<Bond<Vector3>> b;
  
  Edge e_ij;
  bool exists{false};
  for(size_t i{0}; i < n_atoms; i++){
    for(size_t j{0}; j < i; j++){
      
      std::tie(e_ij, exists) = edge(i, j, ug);
      
      if( exists ){
        b.push_back(Bond<Vector3>{i, j,
                                  molecule[i].position,
                                  molecule[j].position,
                                  boost::get(boost::edge_weight, ug, e_ij)});
      }
    }
  }
  
  return b;
}

/*
// TODO: SOLVE BUG ON ANGLE ORDER!
template <typename Vector3, typename Matrix>
std::vector<Angle<Vector3>> angles(const molecule::Molecule<Vector3>& molecule,
                                   const Matrix& connectivity,
                                   double epsilon = 1e-12){
  
  size_t n_atoms{ molecule.size() };
  
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
          d2 = connectivity(k,j);
          
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
  
  return ang;
}
*/

}

#endif //IRC_CONNECTIVITY_H
