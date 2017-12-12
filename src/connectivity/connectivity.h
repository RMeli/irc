#ifndef IRC_CONNECTIVITY_H
#define IRC_CONNECTIVITY_H

#include "../atoms/atom.h"
#include "../atoms/molecule.h"
#include "../linear_algebra/linalg.h"
#include "../tools/constants.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/array.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/graph_traits.hpp>


namespace connectivity {

constexpr double covalent_bond_multiplier{1.3};

using EdgeProperty = boost::property<boost::edge_weight_t, int>;

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

using DistanceProperty = boost::exterior_vertex_property<UGraph, int>;
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
        // The weights are set to 1 for all edges.
        boost::add_edge(i, j, 1, ug);
      }
    }
  }
  
  return ug;
}

/// Find the distance and predecessors matrices of the graph \param ug
/// \tparam Matrix
/// \param ug Graph
/// \return Distance and predecessors matrices
///
/// The element \f$(i,j)\f$ of the distance matrix is an integer indicating
/// how many bonds are between atom \f$i\f$ and atom \f$j\f$, since the
/// weight of each edge is set to 1 in \function adjacency_matrix. This allow
/// to easily determine if two atoms are connected via one bond, two bonds
/// (they form an angle) or three bonds (they form a dihedral).
/// The element \f$(i,j)\f% of the predecessors matrix is an integer indicating
/// the index of the second to last vertex in the shortest path from i to j.
/// This information allow to reconstruct the shortest path from i to j.
template<typename Matrix>
std::pair<Matrix,Matrix> distance_matrix(const UGraph& ug){
  
  using namespace boost;

  // Store number of vertices (number of atoms)
  const size_t n_vertices{ boost::num_vertices(ug) };
  
  // Allocate distance matrix
  Matrix dist{ linalg::zeros<Matrix>(n_vertices, n_vertices) };
  
  // Allocate predecessors matrix
  Matrix predecessors{ linalg::zeros<Matrix>(n_vertices, n_vertices) };
  
  // Allocate distance map for single-source problem
  std::vector<int> d_map(n_vertices, 0);
  
  // Allocate predecessors map for single-source problem
  std::vector<int> p_map(n_vertices, 0);
  
  // Loop over vetrices
  for(size_t i{0}; i < n_vertices; i++){
    // Solve single-source problem for every vertex
    dijkstra_shortest_paths(ug, i,
                            distance_map(&d_map[0]).predecessor_map(&p_map[0]));
    
    // Store distance and predecessors maps
    for(size_t j{0}; j < n_vertices; j++){
      // Fill distance matrix
      dist(i,j) = d_map[j];
      
      // Fill predecessors matrix
      if(i != j){
        predecessors(i,j) = p_map[j];
      }
      else{
        predecessors(i,j) = -1;
      }
    }
  }
  
  // Return distance matrix
  return std::make_pair(dist, predecessors);
}

/// Returns the covalent covalent bonds in \param molecule
/// \tparam Vector3 3D vevtor
/// \tparam Matrix Matrix
/// \param distance_m Distance matrix
/// \param molecule Molecule
/// \return List of covalent bonds
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
        // Compute distance between atom i and atom j
        d = distance(molecule[i].position, molecule[j].position);
        
        // Store bond informations between atom i and atom j
        b.push_back(Bond<Vector3>{i, j,
                                  molecule[i].position,
                                  molecule[j].position,
                                  d});
      }
    }
  }
  
  return b;
}

/// Returns the angles between atoms bonded via a covalent bond in \param
/// molecule
///
/// \tparam Vector3 3D vector
/// \tparam Matrix Matrix
/// \param distance_m Distance matrix
/// \param predecessors_m Matrix of predecessors
/// \param molecule Molecule
/// \return List of angles
template <typename Vector3, typename Matrix>
std::vector<Angle<Vector3>> angles(const Matrix& distance_m,
                                   const Matrix& predecessors_m,
                                   const molecule::Molecule<Vector3>& molecule){
  
  // Extract number of atoms
  const size_t n_atoms{ molecule.size() };
  
  // Declare list of angles
  std::vector<Angle<Vector3>> ang;
  
  size_t k{0};
  double a{0.};
  for(size_t j{0}; j < n_atoms; j++){
    for(size_t i{0}; i < j; i++){
      
      if( distance_m(i,j) == 2){
        k = predecessors_m(i,j);
        
        // Compute angle (i,k,j)
        a = angle(molecule[i].position,
                  molecule[k].position,
                  molecule[j].position);
        
        // Store angle
        ang.push_back( Angle<Vector3>{i, k, j,
                                      molecule[i].position,
                                      molecule[k].position,
                                      molecule[j].position, a} );
      }
    }
  }
  
  return ang;
}

}

#endif //IRC_CONNECTIVITY_H
