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
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>


namespace connectivity {

constexpr double covalent_bond_multiplier{1.3};

using EdgeProperty = boost::property<boost::edge_weight_t, size_t>;

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

using DistanceProperty = boost::exterior_vertex_property<UGraph, size_t>;

using DistanceMatrix = DistanceProperty::matrix_type;

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
Matrix distances(const molecule::Molecule<Vector3>& molecule){
  const size_t n_atoms{ molecule.size() };
  
  Matrix distances_m{ linalg::zeros<Matrix>(n_atoms, n_atoms) };
  
  double r{0.};
  for(size_t i{0}; i < n_atoms; i++){
    for(size_t j{0}; j < n_atoms; j++){
      
      r = distance(molecule[i].position, molecule[j].position);
      
      distances_m(i,j) = r;
      distances_m(j,i) = r;
    }
  }
  
  return distances_m;
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
  const size_t n_atoms{ molecule.size() };
  
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
        boost::add_edge(i, j, 1, ug);
      }
    }
  }
  
  return ug;
}

template<typename Matrix>
Matrix distance_matrix(const UGraph& ug){

  const size_t n_vertices{ boost::num_vertices(ug) };

  // Allocate distance matrix
  DistanceMatrix d{n_vertices};
  
  // Find shortest distance between every pair of vertices in the graph
  boost::johnson_all_pairs_shortest_paths(ug, d);
  
  // Allocate Matrix
  Matrix dist{ linalg::zeros<Matrix>(n_vertices, n_vertices) };
  
  // Compy DistanceMatrix in standard Matrix
  for(size_t j{0}; j < n_vertices; j++){
    for(size_t i{0}; i < n_vertices; i++){
      dist(i,j) = d[i][j];
    }
  }
  
  // Return distance matrix
  return dist;
}

template <typename Vector3, typename Matrix>
std::vector<Bond<Vector3>> bonds(const Matrix& distance_m,
                                 const molecule::Molecule<Vector3>& molecule){
  
  // Extract number of atoms
  const size_t n_atoms{ molecule.size() };
  
  // Declare bond list
  std::vector<Bond<Vector3>> b;
  
  double d{0};
  for(size_t j{0}; j < n_atoms; j++){
    for(size_t i{0}; i < j; i++){
      
      if( distance_m(i,j) == 1 ){
        d = distance(molecule[i].position, molecule[j].position);
        
        b.push_back(Bond<Vector3>{i, j,
                                  molecule[i].position,
                                  molecule[j].position,
                                  d});
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
  
  const size_t n_atoms{ molecule.size() };
  
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
