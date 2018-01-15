#ifndef IRC_CONNECTIVITY_H
#define IRC_CONNECTIVITY_H

#include "atom.h"
#include "molecule.h"
#include "linalg.h"
#include "constants.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/graph_traits.hpp>

namespace irc {

namespace connectivity {

constexpr double covalent_bond_multiplier{1.3};
constexpr double vdw_bond_multiplier{0.9};

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

/// Couple of atoms forming a bond
/// Triplet of atoms forming an angle
///
/// Atoms are represented by their index in a list of coordinates.
struct Bond {
  size_t i;
  size_t j;
};

/// Triplet of atoms forming an angle
///
/// Atoms are represented by their index in a list of coordinates.
struct Angle {
  size_t i;
  size_t j;
  size_t k;
};

/// Quadruplet of atoms forming an angle
///
/// Atoms are represented by their index in a list of coordinates.
struct Dihedral {
  size_t i;
  size_t j;
  size_t k;
  size_t l;
};

/// Compute the distance between two points
///
/// \tparam Vector3
/// \param v1 Point 1
/// \param v2 Poin 2
/// \return Distance between point  1 and point 2
template<typename Vector3>
inline double distance(const Vector3 &v1, const Vector3 &v2) {
  return linalg::norm(v1 - v2);
}

/// Compute bond length
///
/// \tparam Vector3
/// \tparam Vector
/// \param b Bond
/// \param x_cartesian Atomic cartesian coordinates
/// \return Bond length
///
/// Given a (linear) vector of cartesian atomic coordinates \param x_cartesian
/// and a bond \param b, the corresponding bond length is computed.
template<typename Vector3, typename Vector>
inline double bond(const Bond& b, const Vector& x_cartesian){
  // Temporary positions
  Vector3 b1, b2;
  
  // TODO
  b1 = {x_cartesian(3 * b.i + 0), x_cartesian(3 * b.i + 1), x_cartesian(3 * b.i + 2)};
  b2 = {x_cartesian(3 * b.j + 0), x_cartesian(3 * b.j + 1), x_cartesian(3 * b.j + 2)};
  
  return distance(b1, b2);
}

/// Compute bond length
///
/// \tparam Vector3
/// \param b Bond
/// \param molecule Molecule
/// \return Bond length
template<typename Vector3>
inline double bond(const Bond& b, const molecule::Molecule<Vector3>& molecule){
  Vector3 b1{ molecule[b.i].position };
  Vector3 b2{ molecule[b.j].position };
  
  return distance(b1, b2);
}

/// Compute angle formed by three points
///
/// \tparam Vector3
/// \param v1 Point 1
/// \param v2 Point 2
/// \param v3 Point 3
/// \return Angle between points 1, 2 and 3
template<typename Vector3>
inline double angle(const Vector3 &v1, const Vector3 &v2, const Vector3 &v3) {
  Vector3 r1{v1 - v2};
  Vector3 r2{v3 - v2};
  
  double N{linalg::norm(r1) * linalg::norm(r2)};
  
  double angle{std::acos(linalg::dot(r1, r2) / N)};
  
  return angle * 180.0 / tools::constants::pi;
}

/// Compute angle
///
/// \tparam Vector3
/// \tparam Vector
/// \param a Angle
/// \param x_cartesian Atomic cartesian coordinates
/// \return Angle
///
/// Given a (linear) vector of cartesian atomic coordinates \param x_cartesian
/// and a bond \param b, the corresponding bond length is computed.
template<typename Vector3, typename Vector>
inline double angle(const Angle& a, const Vector& x_cartesian){
  // Temporary positions
  Vector3 a1, a2, a3;
  
  // TODO
  a1 = {x_cartesian(3 * a.i + 0), x_cartesian(3 * a.i + 1), x_cartesian(3 * a.i + 2)};
  a2 = {x_cartesian(3 * a.j + 0), x_cartesian(3 * a.j + 1), x_cartesian(3 * a.j + 2)};
  a3 = {x_cartesian(3 * a.k + 0), x_cartesian(3 * a.k + 1), x_cartesian(3 * a.k + 2)};
  
  return angle(a1, a2, a3);
}

/// Compute angle
///
/// \tparam Vector3
/// \param a Angle
/// \param molecule Molecule
/// \return Angle
template<typename Vector3>
inline double angle(const Angle& a,
                    const molecule::Molecule<Vector3>& molecule){
  Vector3 a1{ molecule[a.i].position };
  Vector3 a2{ molecule[a.j].position };
  Vector3 a3{ molecule[a.k].position };

  return angle(a1, a2, a3);
}

/// Compute dihedral angle formed by four points
///
/// \tparam Vector3
/// \param v1 Point 1
/// \param v2 Point 2
/// \param v3 Point 3
/// \param v4 Point 4
/// \return Dihedral angle
template<typename Vector3>
inline double dihedral(const Vector3 &v1,
                       const Vector3 &v2,
                       const Vector3 &v3,
                       const Vector3 &v4) {
  Vector3 b1{v1 - v2};
  Vector3 b2{v2 - v3};
  Vector3 b3{v3 - v4};
  
  Vector3 n1{linalg::cross(b1, b2)};
  Vector3 n2{linalg::cross(b2, b3)};
  
  n1 /= linalg::norm(n1);
  n2 /= linalg::norm(n2);
  
  Vector3 m{linalg::cross(n1, b2) / linalg::norm(b2)};
  
  double x{linalg::dot(n1, n2)};
  double y{linalg::dot(m, n2)};
  
  // Compute dihedral angle in radians (in the intervale [-pi,pi])
  double angle{std::atan2(y, x)};
  
  // Convert angle from radians to degrees
  angle *= 180.0 / tools::constants::pi;
  
  return angle;
}

/// Compute dihedral angle \param d, given cartesian coordinates
///
/// \tparam Vector3
/// \tparam Vector
/// \param d Dihedral
/// \param x_cartesian Atomic cartesian coordinates
/// \return Dihedral angle
///
/// Given a (linear) vector of cartesian atomic coordinates \param x_cartesian
/// and a bond \param b, the corresponding bond length is computed.
template<typename Vector3, typename Vector>
inline double dihedral(const Dihedral& d, const Vector& x_cartesian){
  // Temporary positions
  Vector3 d1, d2, d3, d4;
  
  // TODO
  d1 = {x_cartesian(3 * d.i + 0), x_cartesian(3 * d.i + 1), x_cartesian(3 * d.i + 2)};
  d2 = {x_cartesian(3 * d.j + 0), x_cartesian(3 * d.j + 1), x_cartesian(3 * d.j + 2)};
  d3 = {x_cartesian(3 * d.k + 0), x_cartesian(3 * d.k + 1), x_cartesian(3 * d.k + 2)};
  d4 = {x_cartesian(3 * d.l + 0), x_cartesian(3 * d.l + 1), x_cartesian(3 * d.l + 2)};
  
  return dihedral(d1, d2, d3, d4);
}

/// Compute dihedral angle \param d, given a molecule
///
/// \tparam Vector3
/// \param d Dihedral
/// \param molecule Molecule
/// \return Dihedral angle
template<typename Vector3>
inline double dihedral(const Dihedral& d,
                       const molecule::Molecule<Vector3>& molecule){
  Vector3 d1{ molecule[d.i].position };
  Vector3 d2{ molecule[d.j].position };
  Vector3 d3{ molecule[d.k].position };
  Vector3 d4{ molecule[d.l].position };
  
  return dihedral(d1, d2, d3, d4);
}

/// Compute all distances between atoms in \param molecule
///
/// \tparam Vector3 3D vector
/// \tparam Matrix Matrix
/// \param molecule Molecule (collection of atoms)
/// \return Distances
///
/// The distances \f$\mathbf{D}\f$ for atoms in \param molecule is given by
/// \f[
///   D_{ij} = |\mathbf{r}_i - \mathbf{r}_j|,
/// \f]
/// where the matrix element \f$D_{ij}\f$ is the distance between atom at
/// position \f$\mathbf{r}_i\f$ and the atom at position \f$\mathbf{r}_j\f$.
template<typename Vector3, typename Matrix>
Matrix distances(const molecule::Molecule<Vector3> &molecule) {
  const size_t n_atoms{molecule.size()};
  
  Matrix distances_m{linalg::zeros<Matrix>(n_atoms, n_atoms)};
  
  double r{0.};
  for (size_t i{0}; i < n_atoms; i++) {
    for (size_t j{0}; j < n_atoms; j++) {
      
      r = distance(molecule[i].position, molecule[j].position);
      
      distances_m(i, j) = r;
      distances_m(j, i) = r;
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
/// number of edges is determined by bonding.
template<typename Vector3, typename Matrix>
UGraph adjacency_matrix(const Matrix &distance_m,
                        const molecule::Molecule<Vector3> &molecule) {
  
  // Extract number of atoms
  const size_t n_atoms{molecule.size()};
  
  // Define a undirected graph with n_atoms vertices
  UGraph ug(n_atoms);
  
  double d{0.};
  double sum_covalent_radii{0.};
  double sum_vdw_radii{0.};
  for (size_t j{0}; j < n_atoms; j++) {
    for (size_t i{j + 1}; i < n_atoms; i++) {
      
      // Extract distance between atom i and atom j
      d = distance_m(i, j);
      
      // Compute sum of covalent radii for atoms i and j
      sum_covalent_radii = atom::covalent_radius(molecule[i].atomic_number) +
                           atom::covalent_radius(molecule[j].atomic_number);
      
      // Determine if atoms i and j are bonded
      if (d < covalent_bond_multiplier * sum_covalent_radii) {
        // Add edge to boost::adjacency_list between vertices i and j
        // The weights are set to 1 for all edges.
        boost::add_edge(i, j, 1, ug);
        
        // TODO: Better ways of doing this...
        // Search for H-bonds: XH...Y
        if ((atom::is_NOFPSCl(molecule[i].atomic_number) and
             atom::is_H(molecule[j].atomic_number))
            or
            (atom::is_NOFPSCl(molecule[j].atomic_number) and
             atom::is_H(molecule[i].atomic_number))) { // Possible H-bond
          
          size_t idx{0}; // X atom index
          size_t h_idx{0}; // Hydrogen bond index
          
          double a{0}; // Angle between X, H and Y in XH...Y
          
          // Assign correct indices to X and H
          if (atom::is_H(molecule[j].atomic_number)) {
            idx = i;
            h_idx = j;
          } else {
            idx = j;
            h_idx = i;
          }
          
          // Loop over all other atoms, excluding i and j, to find Y
          for (size_t k{0}; k < n_atoms; k++) {
            if (atom::is_NOFPSCl(molecule[k].atomic_number) and
                k != idx and k != h_idx) {
              
              // Load distance
              d = distance_m(h_idx, k);
              
              // Compute sum of Van der Waals radii
              sum_vdw_radii =
                  atom::vdw_radius(molecule[h_idx].atomic_number) +
                  atom::vdw_radius(molecule[k].atomic_number);
              
              // Compute sum of covalent radii
              sum_covalent_radii =
                  atom::covalent_radius(molecule[h_idx].atomic_number) +
                  atom::covalent_radius(molecule[k].atomic_number);
              
              a = angle(molecule[idx].position,
                        molecule[h_idx].position,
                        molecule[k].position);
              
              // Check H-bond properties
              if (d > sum_covalent_radii and
                  d < sum_vdw_radii * vdw_bond_multiplier and
                  a > 90) {
                // Add hydrogen bond
                boost::add_edge(h_idx, k, 1, ug);
              }
            }
          } // End H-bond search
        }
      }
    }
  }
  
  // Allocate storage for fragment indices
  std::vector<size_t> fragments( boost::num_vertices(ug) );
  
  // Fill component std::vector and return number of different fragments
  // If num_fragments == 1 the graph is connected
  size_t num_fragments{ boost::connected_components (ug, &fragments[0])};
  
  if(num_fragments > 1){
    // Print fragments
    for(size_t idx : fragments) {
      std::cout << idx << ' ';
    }
    std::cout << std::endl;
    
    // TODO: Support fragments
    throw std::logic_error("Fragment recognition not implemented.");
  }
  
  return ug;
}

/// Find the distance and predecessors matrices of the graph \param ug
/// \tparam Matrix
/// \param ug Graph
/// \return Distance and predecessors matrices
///
/// The element \f$(i,j)\f$ of the distance matrix is an integer indicating
/// how many bonds (along the shortest path) are between atom \f$i\f$ and atom
/// \f$j\f$, since the weight of each edge is set to 1 in \function
/// adjacency_matrix. This allow to easily determine if two atoms are connected
/// via one bond, two bonds (they form an angle) or three bonds (they form a
/// dihedral). The element \f$(i,j)\f% of the predecessors matrix is an integer
/// indicating the index of the second to last vertex in the shortest path from
/// i to j. This information allow to reconstruct the shortest path from i to j.
template<typename Matrix>
std::pair<Matrix, Matrix> distance_matrix(const UGraph &ug) {
  
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
  for (size_t i{0}; i < n_vertices; i++) {
    // Solve single-source problem for every vertex
    dijkstra_shortest_paths(ug, i,
                            distance_map(&d_map[0]).predecessor_map(&p_map[0]));
    
    // Store distance and predecessors maps
    for (size_t j{0}; j < n_vertices; j++) {
      // Fill distance matrix
      dist(i, j) = d_map[j];
      
      // Fill predecessors matrix
      if (i != j) {
        predecessors(i, j) = p_map[j];
      } else {
        predecessors(i, j) = -1;
      }
    }
  }
  
  // Return distance matrix
  return std::make_pair(dist, predecessors);
}

/// Returns the bonds in \param molecule
///
/// \tparam Vector3
/// \tparam Matrix
/// \param distance_m Distance matrix
/// \param molecule Molecule
/// \return List of bonds
///
/// The bonds can be covalent bonds, hydrogen bonds or inter-fragment bonds.
template<typename Vector3, typename Matrix>
std::vector<Bond> bonds(const Matrix &distance_m,
                        const molecule::Molecule<Vector3> &molecule) {
  
  // Extract number of atoms
  const size_t n_atoms{molecule.size()};
  
  // Declare bond list
  std::vector<Bond> b;

  for (size_t j{0}; j < n_atoms; j++) {
    for (size_t i{0}; i < j; i++) {
      
      if (distance_m(i, j) == 1) {
        // Store bond informations between atom i and atom j
        b.push_back(Bond{i, j});
      }
    }
  }
  
  // Return list of bonds
  return b;
}

/// Returns the angles between bonded atoms in \param molecule
///
/// \tparam Vector3
/// \tparam Matrix
/// \param distance_m Distance matrix
/// \param predecessors_m Matrix of predecessors
/// \param molecule Molecule
/// \return List of angles
template<typename Vector3, typename Matrix>
std::vector<Angle> angles(const Matrix &distance_m,
                          const Matrix &predecessors_m,
                          const molecule::Molecule<Vector3> &molecule){
  
  // Extract number of atoms
  const size_t n_atoms{molecule.size()};
  
  // Declare list of angles
  std::vector<Angle> ang;
  
  size_t k{0};
  for (size_t j{0}; j < n_atoms; j++) {
    for (size_t i{0}; i < j; i++) {
      
      if (distance_m(i, j) == 2) {
        k = predecessors_m(i, j);
        
        // Store angle
        ang.push_back(Angle{i, k, j});
      }
    }
  }
  
  // Return list of angles
  return ang;
}

/// Returns the dihedral angles between bonded atoms in \param molecule
///
/// \tparam Vector3
/// \tparam Matrix
/// \param distance_m Distance matrix
/// \param predecessors_m Matrix of predecessors
/// \param molecule Molecule
/// \return List of dihedral angles
template<typename Vector3, typename Matrix>
std::vector<Dihedral> dihedrals(const Matrix &distance_m,
                                const Matrix &predecessors_m,
                                const molecule::Molecule<Vector3> &molecule) {
  
  // Extract number of atoms
  const size_t n_atoms{molecule.size()};
  
  // Declare list of dihedrals
  std::vector<Dihedral> dih;
  
  size_t k{0}, l{0};
  for (size_t j{0}; j < n_atoms; j++) {
    for (size_t i{0}; i < j; i++) {
      
      if (distance_m(i, j) == 3) {
        k = predecessors_m(i, j);
        l = predecessors_m(i, k);
        
        // Store dihedral angle
        dih.push_back(Dihedral{i, l, k, j});
      }
    }
  }
  
  // Return list of dihedral angles
  return dih;
}

} // namespace connectivity

} // namespace irc

#endif //IRC_CONNECTIVITY_H
