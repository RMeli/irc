#ifndef IRC_CONNECTIVITY_H
#define IRC_CONNECTIVITY_H

#include "atom.h"
#include "constants.h"
#include "linalg.h"
#include "molecule.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/math/special_functions/round.hpp>

namespace irc {

/// Connectivity
namespace connectivity {

using EdgeProperty = boost::property<boost::edge_weight_t, int>;

using UGraph = boost::adjacency_list<boost::vecS,        //
                                     boost::vecS,        //
                                     boost::undirectedS, // Graph type
                                     boost::no_property, // Vertex property
                                     EdgeProperty        // Edge property
                                     >;

using Vertex = boost::graph_traits<UGraph>::vertex_descriptor;
using Edge = boost::graph_traits<UGraph>::edge_descriptor;

using DistanceProperty = boost::exterior_vertex_property<UGraph, int>;
using DistanceMatrix = DistanceProperty::matrix_type;

enum class Constraint { constrained, unconstrained };

/// Couple of atoms forming a bond
/// Triplet of atoms forming an angle
///
/// Atoms are represented by their index in a list of coordinates.
class Bond {
public:
  Bond(const std::size_t& i_,
       const std::size_t& j_,
       const Constraint& constraint_ = Constraint::unconstrained)
    : i(i_), j(j_), constraint(constraint_) {
    if (i == j) {
      throw std::logic_error("Bond error.");
    }

    // Ordering (needed for hash function and comparison operators)
    if (j < i) {
      std::swap(i, j);
    }
  }

  std::size_t i;
  std::size_t j;

  Constraint constraint{Constraint::unconstrained};

  bool operator==(const Bond& b) const { return i == b.i && j == b.j; }

  bool operator!=(const Bond& b) const { return !(*this == b); }
};

/// Triplet of atoms forming an angle
///
/// Atoms are represented by their index in a list of coordinates.
class Angle {
public:
  Angle(const std::size_t& i_,
        const std::size_t& j_,
        const std::size_t& k_,
        const Constraint& constraint_ = Constraint::unconstrained)
    : i(i_), j(j_), k(k_), constraint(constraint_) {
    if (i == j || i == k || j == k) {
      throw std::logic_error("Angle error.");
    }

    // Ordering (needed for hash function and comparison operators)
    if (k < i) {
      std::swap(i, k);
    }
  }

  std::size_t i;
  std::size_t j;
  std::size_t k;
  // enum class linear{NONE, XY, YZ}; // Use switch

  Constraint constraint{Constraint::unconstrained};

  bool operator==(const Angle& a) const {
    return i == a.i && j == a.j && k == a.k;
  }

  bool operator!=(const Angle& a) const { return !(*this == a); }
};

enum class LinearAngleTag { First, Second };

/// \brief String form of \p tag
std::string inline to_string(const LinearAngleTag tag) {
  switch (tag) {
    case LinearAngleTag::First:
      return "First";
    case LinearAngleTag::Second:
      return "Second";
  }
}

/// \brief Triplet of atoms forming a linear angle and an orthogonal direction
///
/// Atoms are represented by their index in a list of coordinates.
/// \p orthogonal_direction_ vector should be orthogonal to the \p i_ to \p k_
/// atoms directions.
/// A \p tag_ is also required, because the linear angles should be defined in
/// pairs such that bending can occur in any direction.
/// The \p orthogonal_direction_ of these two must themselves be orthogonal.
template<typename Vector3>
class LinearAngle {
public:
  LinearAngle(const std::size_t i_,
              const std::size_t j_,
              const std::size_t k_,
              const Vector3 orthogonal_direction_,
              const LinearAngleTag tag_,
              const Constraint constraint_ = Constraint::unconstrained)
    : i(i_), j(j_), k(k_), tag(tag_),
      orthogonal_direction(orthogonal_direction_), constraint(constraint_) {
    if (i == j || i == k || j == k) {
      throw std::logic_error("Angle error.");
    }

    // Ordering (needed for hash function and comparison operators)
    if (k < i) {
      std::swap(i, k);
    }
  }

  std::size_t i;
  std::size_t j;
  std::size_t k;
  LinearAngleTag tag;
  Vector3 orthogonal_direction;

  Constraint constraint{Constraint::unconstrained};

  bool operator==(const LinearAngle<Vector3>& a) const {
    return i == a.i && j == a.j && k == a.k && tag == a.tag;
  }

  bool operator!=(const LinearAngle<Vector3>& a) const { return !(*this == a); }
};

/// Quadruplet of atoms forming an angle
///
/// Atoms are represented by their index in a list of coordinates.
class Dihedral {
public:
  Dihedral(const std::size_t& i_,
           const std::size_t& j_,
           const std::size_t& k_,
           const std::size_t& l_,
           const Constraint& constraint_ = Constraint::unconstrained)
    : i(i_), j(j_), k(k_), l(l_), constraint(constraint_) {
    if (i == j || i == k || i == l || j == k || j == l || k == l) {
      throw std::logic_error("Dihedral error.");
    }

    // Ordering (needed for hash function and comparison operators)
    if (l < i) {
      std::swap(i, l);
      std::swap(j, k);
    }
  }

  std::size_t i;
  std::size_t j;
  std::size_t k;
  std::size_t l;

  Constraint constraint{Constraint::unconstrained};

  bool operator==(const Dihedral& d) const {
    return i == d.i && j == d.j && k == d.k && l == d.l;
  }

  bool operator!=(const Dihedral& d) const { return !(*this == d); }
};

/// Quadruplet of atoms forming an out of plane bend
///
/// Atoms are represented by their index in a list of coordinates.
class OutOfPlaneBend {
public:
  OutOfPlaneBend(const std::size_t& c_,
                 const std::size_t& i_,
                 const std::size_t& j_,
                 const std::size_t& k_,
                 const Constraint& constraint_ = Constraint::unconstrained)
      : c(c_), i(i_), j(j_), k(k_), constraint(constraint_) {
    if (c == i || c == j || c == k || i == j || i == k || j == k) {
      throw std::logic_error("OutOfPlaneBend error.");
    }

    using std::swap;
    // Sort such that i <= j <= k
    if (i > k) {
      swap(i, k);
    }

    if (i > j) {
      swap(i, j);
    }

    if (j > k) {
      swap(j, k);
    }
  }

  std::size_t c; ///<! Central atom
  std::size_t i;
  std::size_t j;
  std::size_t k;

  Constraint constraint{Constraint::unconstrained};

  bool operator==(const OutOfPlaneBend& d) const {
    return c == d.c && i == d.i && j == d.j && k == d.k;
  }

  bool operator!=(const OutOfPlaneBend& d) const { return !(*this == d); }
};

/// Compute the distance between two points
///
/// \tparam Vector3
/// \param v1 Point 1
/// \param v2 Point 2
/// \return Distance between point  1 and point 2
///
/// The distance \f$d\f$ between two points \f$\vec{v}_1\f$ and \f$\vec{v}_2\f$
/// is defined as
/// \f[
///   d = |\vec{v}_1 - \vec{v}_2|.
/// \f]
template<typename Vector3>
inline double distance(const Vector3& v1, const Vector3& v2) {
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
inline double bond(const Bond& b, const Vector& x_cartesian) {
  // Temporary positions
  const Vector3 b1{x_cartesian(3 * b.i + 0),
                   x_cartesian(3 * b.i + 1),
                   x_cartesian(3 * b.i + 2)};

  const Vector3 b2{x_cartesian(3 * b.j + 0),
                   x_cartesian(3 * b.j + 1),
                   x_cartesian(3 * b.j + 2)};

  return distance(b1, b2);
}

/// Compute bond length
///
/// \tparam Vector3
/// \param b Bond
/// \param molecule Molecule
/// \return Bond length
template<typename Vector3>
inline double bond(const Bond& b, const molecule::Molecule<Vector3>& molecule) {
  const Vector3 b1{molecule[b.i].position};
  const Vector3 b2{molecule[b.j].position};

  return distance(b1, b2);
}

/// Compute angle formed by three points
///
/// \tparam Vector3
/// \param v1 Point 1
/// \param v2 Point 2
/// \param v3 Point 3
/// \return Angle between points 1, 2 and 3
///
/// The angle \f$a\f$ between three points \f$\vec{v}_1\f$, \f$\vec{v}_2\f$
/// and \f$\vec{v}_3\f$ is defined as
/// \f[
///   a = \cos^{-1}\left( \frac{\vec{r}_{21}\cdot\vec{r}_{23}}
///       {|\vec{r}_{21}||\vec{r}_{23}|} \right).
/// \f]
/// where \f$\vec{r}_{21}=\vec{v}_1-\vec{v}_2\f$ and
/// \f$\vec{r}_{23}= \vec{v}_3-\vec{v}_2\f$
template<typename Vector3>
inline double angle(const Vector3& v1, const Vector3& v2, const Vector3& v3) {
  const Vector3 r1{v1 - v2};
  const Vector3 r2{v3 - v2};

  const double N{linalg::norm(r1) * linalg::norm(r2)};

  const double dot_r1r2 = linalg::dot(r1, r2);
  if (dot_r1r2 / N <= -1.0) {
    return tools::constants::pi;
  } else if (dot_r1r2 / N >= 1.0) {
    return 0.0;
  } else {
    return std::acos(linalg::dot(r1, r2) / N);
  }
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
inline double angle(const Angle& a, const Vector& x_cartesian) {
  // Temporary positions
  const Vector3 a1{x_cartesian(3 * a.i + 0),
                   x_cartesian(3 * a.i + 1),
                   x_cartesian(3 * a.i + 2)};

  const Vector3 a2{x_cartesian(3 * a.j + 0),
                   x_cartesian(3 * a.j + 1),
                   x_cartesian(3 * a.j + 2)};

  const Vector3 a3{x_cartesian(3 * a.k + 0),
                   x_cartesian(3 * a.k + 1),
                   x_cartesian(3 * a.k + 2)};

  return angle(a1, a2, a3);
}

template<typename Vector3, typename Vector>
inline double angle(const LinearAngle<Vector3>& a, const Vector& x_cartesian) {

  // Temporary positions
  const Vector3 a1{x_cartesian(3 * a.i + 0),
                   x_cartesian(3 * a.i + 1),
                   x_cartesian(3 * a.i + 2)};

  const Vector3 a2{x_cartesian(3 * a.j + 0),
                   x_cartesian(3 * a.j + 1),
                   x_cartesian(3 * a.j + 2)};

  const Vector3 a3{x_cartesian(3 * a.k + 0),
                   x_cartesian(3 * a.k + 1),
                   x_cartesian(3 * a.k + 2)};

  return angle(a1, a2, a.orthogonal_direction) +
         angle(a.orthogonal_direction, a2, a3);
}

/// Compute angle
///
/// \tparam Vector3
/// \param a Angle
/// \param molecule Molecule
/// \return Angle
template<typename Vector3>
inline double angle(const Angle& a,
                    const molecule::Molecule<Vector3>& molecule) {
  const Vector3 a1{molecule[a.i].position};
  const Vector3 a2{molecule[a.j].position};
  const Vector3 a3{molecule[a.k].position};

  return angle(a1, a2, a3);
}
template<typename Vector3>
inline double angle(const LinearAngle<Vector3>& a,
                    const molecule::Molecule<Vector3>& molecule) {
  const Vector3 a1{molecule[a.i].position};
  const Vector3 a2{molecule[a.j].position};
  const Vector3 a3{molecule[a.k].position};

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
inline double dihedral(const Vector3& v1,
                       const Vector3& v2,
                       const Vector3& v3,
                       const Vector3& v4) {
  const Vector3 b1{v1 - v2};
  const Vector3 b2{v2 - v3};
  const Vector3 b3{v3 - v4};

  Vector3 n1{linalg::cross(b1, b2)};
  Vector3 n2{linalg::cross(b2, b3)};

  n1 /= linalg::norm(n1);
  n2 /= linalg::norm(n2);

  const Vector3 m{linalg::cross(n1, b2) / linalg::norm(b2)};

  const double x{linalg::dot(n1, n2)};
  const double y{linalg::dot(m, n2)};

  // Compute dihedral angle in radians (in the interval [-pi,pi])
  const double angle{std::atan2(y, x)};

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
inline double dihedral(const Dihedral& d, const Vector& x_cartesian) {
  // Temporary positions
  const Vector3 d1{x_cartesian(3 * d.i + 0),
                   x_cartesian(3 * d.i + 1),
                   x_cartesian(3 * d.i + 2)};

  const Vector3 d2{x_cartesian(3 * d.j + 0),
                   x_cartesian(3 * d.j + 1),
                   x_cartesian(3 * d.j + 2)};

  const Vector3 d3{x_cartesian(3 * d.k + 0),
                   x_cartesian(3 * d.k + 1),
                   x_cartesian(3 * d.k + 2)};

  const Vector3 d4{x_cartesian(3 * d.l + 0),
                   x_cartesian(3 * d.l + 1),
                   x_cartesian(3 * d.l + 2)};

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
                       const molecule::Molecule<Vector3>& molecule) {
  const Vector3 d1{molecule[d.i].position};
  const Vector3 d2{molecule[d.j].position};
  const Vector3 d3{molecule[d.k].position};
  const Vector3 d4{molecule[d.l].position};

  return dihedral(d1, d2, d3, d4);
}


/// Compute dihedral angle formed by four points
template<typename Vector3>
double out_of_plane_angle(const Vector3& vc,
                          const Vector3& v1,
                          const Vector3& v2,
                          const Vector3& v3) {
  const Vector3 b1{v1 - vc};
  const Vector3 b2{v2 - vc};
  const Vector3 b3{v3 - vc};

  const Vector3 e1{linalg::normalize(b1)};
  const Vector3 e2{linalg::normalize(b2)};
  const Vector3 e3{linalg::normalize(b3)};

  const double a1{angle(v2, vc, v3)};
  const Vector3 t1{linalg::cross(e2, e3)/std::sin(a1)};

  return std::asin(linalg::dot(t1, e1));
}

/// Compute dihedral angle \param d, given cartesian coordinates
///
/// Given a (linear) vector of cartesian atomic coordinates \param x_cartesian
/// and a bond \param b, the corresponding bond length is computed.
template<typename Vector3, typename Vector>
double out_of_plane_angle(const OutOfPlaneBend& d,
                          const Vector& x_cartesian) {
  // Temporary positions
  const Vector3 dc{x_cartesian(3 * d.c + 0),
                   x_cartesian(3 * d.c + 1),
                   x_cartesian(3 * d.c + 2)};

  const Vector3 d1{x_cartesian(3 * d.i + 0),
                   x_cartesian(3 * d.i + 1),
                   x_cartesian(3 * d.i + 2)};

  const Vector3 d2{x_cartesian(3 * d.j + 0),
                   x_cartesian(3 * d.j + 1),
                   x_cartesian(3 * d.j + 2)};

  const Vector3 d3{x_cartesian(3 * d.k + 0),
                   x_cartesian(3 * d.k + 1),
                   x_cartesian(3 * d.k + 2)};

  return out_of_plane_angle(dc, d1, d2, d3);
}

/// Compute dihedral angle \param d, given a molecule
template<typename Vector3>
double out_of_plane_angle(const OutOfPlaneBend& d,
                          const molecule::Molecule<Vector3>& molecule) {

  const Vector3 dc{molecule[d.c].position};
  const Vector3 d1{molecule[d.i].position};
  const Vector3 d2{molecule[d.j].position};
  const Vector3 d3{molecule[d.k].position};

  return out_of_plane_angle(dc, d1, d2, d3);
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
Matrix distances(const molecule::Molecule<Vector3>& molecule) {
  const std::size_t n_atoms{molecule.size()};

  Matrix distances_m{linalg::zeros<Matrix>(n_atoms, n_atoms)};

  double r{0.};
  for (std::size_t j{0}; j < n_atoms; j++) {
    for (std::size_t i{0}; i < j; i++) {

      r = distance(molecule[i].position, molecule[j].position);

      distances_m(i, j) = r;
      distances_m(j, i) = r;
    }
  }

  return distances_m;
}

// TODO: Improve algorithm
// TODO: Test
/*!
 *
 * @tparam Matrix
 * @param i First fragment index
 * @param j Second fragment index
 * @param fragments Fragment indices
 * @param distances Distance matrix
 * @return Tuple containing the indices of the atoms with minimum interfragment
 * distance (between fragment @param i and fragment @j) and such distance.
 */
template<typename Matrix>
std::tuple<std::size_t, std::size_t, double>
min_interfragment_distance(std::size_t i,
                           std::size_t j,
                           const std::vector<std::size_t>& fragments,
                           const Matrix& distances) {
  const std::size_t n_atoms{fragments.size()};

  // Interfragment distance
  double d{0};

  // Minimal interfragment distance
  double min_distance{std::numeric_limits<double>::max()};

  std::size_t k_min{0}, l_min{0};
  for (std::size_t k{0}; k < n_atoms; k++) {
    for (std::size_t l{0}; l < n_atoms; l++) {
      if (k != l and fragments[k] == i and fragments[l] == j) {

        d = distances(l, k);

        if (d < min_distance) {
          min_distance = d;
          k_min = k;
          l_min = l;
        }
      }
    }
  }

  return std::make_tuple(k_min, l_min, min_distance);
}

// TODO: Test
/*! Identify fragments (connected components)
 *
 * The function returns the number of fragments and a vector containing
 * the index of the fragment each atom belongs to.
 *
 * @param ug Unsigned graph
 * @return Number of fragments and fragments indices
 */
inline std::pair<std::size_t, std::vector<std::size_t>>
identify_fragments(const UGraph& ug) {
  // Fragment indices
  std::vector<std::size_t> fragments(boost::num_vertices(ug));

  // Fill component std::vector and return number of different fragments
  // If num_fragments == 1 the graph is connected
  const std::size_t num_fragments{
      boost::connected_components(ug, &fragments[0])};

  return {num_fragments, fragments};
}

/*! Search for regular bonds (covalent bonds)
 *
 * @tparam Vector3
 * @tparam Matrix
 * @param ug Adjacency matrix
 * @param distances Distance matrix
 * @param molecule Molecule
 */
template<typename Vector3, typename Matrix>
void add_regular_bonds(UGraph& ug,
                       const Matrix& distances,
                       const molecule::Molecule<Vector3>& molecule) {
  const std::size_t n_atoms{molecule.size()};

  double d{0.};

  double sum_covalent_radii{0.};

  for (std::size_t j{0}; j < n_atoms; j++) {
    for (std::size_t i{j + 1}; i < n_atoms; i++) {
      d = distances(i, j);

      sum_covalent_radii = atom::covalent_radius(molecule[i].atomic_number) +
                           atom::covalent_radius(molecule[j].atomic_number);

      // Determine if atoms i and j are bonded
      if (d < tools::constants::covalent_bond_multiplier * sum_covalent_radii) {
        // Add edge to boost::adjacency_list between vertices i and j
        // The weights are set to 1 for all edges.
        boost::add_edge(i, j, 1, ug);
      }
    }
  }
}

// TODO: Improve algorithm
// TODO: Test
/*! Recursive search of interfragment bonds
 *
 * At each iteration of the recursive search the interfragment bonds (and
 * auxiliary interfragment bonds) are added between closest fragments.
 *
 * @tparam Matrix
 * @param ug Adjacency matrix
 * @param distances Distance matrix
 */
template<typename Matrix>
void add_interfragment_bonds(UGraph& ug, const Matrix& distances) {

  const size_t n_atoms{boost::num_vertices(ug)};

  std::size_t num_fragments;
  std::vector<std::size_t> fragments;

  // Get number of fragments and fragment indices
  std::tie(num_fragments, fragments) = identify_fragments(ug);

  while (num_fragments > 1) {

    struct InterfragmentDistance {
      double d;
      size_t i;
      size_t j;
    };

    // Minimum interfragment distances
    // Matrix
    // min_dist_fragments{linalg::zeros<Matrix>(num_fragments,num_fragments)};
    std::vector<std::vector<InterfragmentDistance>> min_dist_fragments(
        num_fragments, std::vector<InterfragmentDistance>(num_fragments));

    // Determine minimal interfragment distances
    std::size_t i_min{0}, j_min{0};
    double min_d{0};
    for (std::size_t j{0}; j < num_fragments; j++) {
      for (std::size_t i{0}; i < j; i++) {
        std::tie(i_min, j_min, min_d) =
            min_interfragment_distance<Matrix>(i, j, fragments, distances);

        min_dist_fragments[i][j] = {min_d, i_min, j_min};
        min_dist_fragments[j][i] = {min_d, i_min, j_min};
      }
    }

    // Add interfragment distances between closest fragments
    double d{0.};
    for (std::size_t j{0}; j < num_fragments; j++) {
      size_t i_min_fragment{0};
      double d_min{std::numeric_limits<double>::max()};
      for (std::size_t i{0}; i < num_fragments; i++) {
        d = min_dist_fragments[i][j].d;
        if (d < d_min && i != j) {
          i_min_fragment = i;
          i_min = min_dist_fragments[i][j].i;
          j_min = min_dist_fragments[i][j].j;
          d_min = d;
        }
      }

      // Add shortest interfragment bond
      boost::add_edge(i_min, j_min, 1, ug);

      // Add auxiliary interfragment distances
      double d{0};
      for (std::size_t k{0}; k < n_atoms; k++) {
        for (std::size_t l{0}; l < n_atoms; l++) {
          if (k != l and fragments[k] == i_min_fragment and fragments[l] == j) {
            d = distances(l, k);

            // TODO: Check
            if (d < std::min(
                        min_d * tools::constants::interfragment_bond_multiplier,
                        2. * tools::conversion::angstrom_to_bohr)) {
              boost::add_edge(l, k, 1, ug);
            }
          }
        }
      }
    }

    // Recursive search of fragments
    // add_interfragment_bonds(ug, distances);
    // Get number of fragments and fragment indices
    std::tie(num_fragments, fragments) = identify_fragments(ug);
  }

  return;
}

// TODO: Better strategy to look for H-bonds (regular bonds are known)
/*! Search for hydrogen bonds
 *
 * @tparam Vector3
 * @tparam Matrix
 * @param ug Adjacency matrix
 * @param distances Distance matrix
 * @param molecule Molecule
 */
template<typename Vector3, typename Matrix>
void add_hydrogen_bonds(UGraph& ug,
                        const Matrix& distances,
                        const molecule::Molecule<Vector3>& molecule) {

  const std::size_t n_atoms{molecule.size()};

  double d{0.};

  double sum_covalent_radii{0.};
  double sum_vdw_radii{0.};
  for (std::size_t j{0}; j < n_atoms; j++) {
    for (std::size_t i{j + 1}; i < n_atoms; i++) {

      d = distances(i, j);

      sum_covalent_radii = atom::covalent_radius(molecule[i].atomic_number) +
                           atom::covalent_radius(molecule[j].atomic_number);

      // Determine if atoms i and j are bonded
      if (d < tools::constants::covalent_bond_multiplier * sum_covalent_radii) {

        // TODO: Better ways of doing this?
        // Search for H-bonds: XH...Y
        if ((atom::is_NOFPSCl(molecule[i].atomic_number) and
             atom::is_H(molecule[j].atomic_number)) or
            (atom::is_NOFPSCl(molecule[j].atomic_number) and
             atom::is_H(molecule[i].atomic_number))) { // Possible H-bond
          // On atom is H, while the other is either N, O, F, P, S or Cl

          std::size_t idx{0};   // X atom index
          std::size_t h_idx{0}; // Hydrogen bond index

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
          for (std::size_t k{0}; k < n_atoms; k++) {
            if (atom::is_NOFPSCl(molecule[k].atomic_number) and k != idx and
                k != h_idx) {

              d = distances(h_idx, k);

              sum_vdw_radii = atom::vdw_radius(molecule[h_idx].atomic_number) +
                              atom::vdw_radius(molecule[k].atomic_number);

              sum_covalent_radii =
                  atom::covalent_radius(molecule[h_idx].atomic_number) +
                  atom::covalent_radius(molecule[k].atomic_number);

              a = angle(molecule[idx].position,
                        molecule[h_idx].position,
                        molecule[k].position);

              // Check H-bond properties
              if (d > sum_covalent_radii and
                  d < sum_vdw_radii * tools::constants::vdw_bond_multiplier and
                  a > tools::constants::pi / 2.) {
                // Add hydrogen bond
                boost::add_edge(h_idx, k, 1, ug);
              }
            }
          }
        }
      }
    }
  } // End search for hydrogen bonds
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
UGraph adjacency_matrix(const Matrix& distances,
                        const molecule::Molecule<Vector3>& molecule) {
  const std::size_t n_atoms{molecule.size()};

  // Define a undirected graph with n_atoms vertices
  UGraph ug(n_atoms);

  add_regular_bonds(ug, distances, molecule);

  add_interfragment_bonds(ug, distances);

  add_hydrogen_bonds(ug, distances, molecule);

  // TODO: Extra redundant coordinates.

  return ug;
}

/// Find the distance and predecessors matrices of the graph \param ug
/// \tparam Matrix
/// \param ug Graph
/// \return Distance matrix
///
/// The element \f$(i,j)\f$ of the distance matrix is an integer indicating
/// how many bonds (along the shortest path) are between atom \f$i\f$ and atom
/// \f$j\f$, since the weight of each edge is set to 1 in \function
/// adjacency_matrix. This allow to easily determine if two atoms are connected
/// via one bond, two bonds (they form an angle) or three bonds (they form a
/// dihedral).
template<typename Matrix>
Matrix distance_matrix(const UGraph& ug) {

  using namespace boost;

  // Store number of vertices (number of atoms)
  const std::size_t n_vertices{boost::num_vertices(ug)};

  // Allocate distance matrix
  Matrix dist{linalg::zeros<Matrix>(n_vertices, n_vertices)};

  // Allocate distance map for single-source problem
  std::vector<int> d_map(n_vertices, 0);

  // Allocate predecessors map for single-source problem
  std::vector<int> p_map(n_vertices, 0);

  // Loop over vertices
  for (std::size_t i{0}; i < n_vertices; i++) {
    // Solve single-source problem for every vertex
    dijkstra_shortest_paths(
        ug, i, distance_map(&d_map[0]).predecessor_map(&p_map[0]));

    // Store distance and predecessors maps
    for (std::size_t j{0}; j < n_vertices; j++) {
      // Fill distance matrix
      dist(i, j) = d_map[j];
    }
  }

  // Return distance matrix
  return dist;
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
std::vector<Bond> bonds(const Matrix& distance_m,
                        const molecule::Molecule<Vector3>& molecule) {

  using boost::math::iround;

  const std::size_t n_atoms{molecule.size()};

  std::vector<Bond> b;

  for (std::size_t j{0}; j < n_atoms; j++) {
    for (std::size_t i{0}; i < j; i++) {

      if (iround(distance_m(i, j)) == 1) {
        b.emplace(b.end(), i, j);
      }
    }
  }

  // Return list of bonds
  return b;
}

/// Determine all possible angles between atoms i and j
///
/// \tparam Matrix
/// \param i
/// \param j
/// \param distance
/// \return
///
/// Dijkstra shortest paths algorithm returns only one shortest path. In some
/// cases however, there might be two different angles between the same two
/// end atoms.
template<typename Matrix>
std::vector<Angle>
angles(std::size_t i, std::size_t j, const Matrix& distance) {

  using boost::math::iround;

  std::vector<Angle> angles;

  const std::size_t n_atoms{linalg::n_rows(distance)};

  for (std::size_t k{0}; k < n_atoms; k++) {
    if (iround(distance(k, i)) == 1 and iround(distance(k, j)) == 1) {
      angles.emplace(angles.end(), i, k, j);
    }
  }

  return angles;
}

/// \brief Returns list from \p angles whose angle is not quasi linear
///
/// Quasi linear angles are those whose angle in \p molecule is greater than
/// tools::constants::quasi_linear_angle.
///
/// \tparam Vector3
/// \param angles Set of potential angles
/// \param molecule Molecule
/// \return List of angles
template<typename Vector3>
std::vector<Angle> valid_angles(const std::vector<Angle>& angles,
                                const molecule::Molecule<Vector3>& molecule) {
  std::vector<Angle> new_angles;

  for (const auto& aa : angles) {
    const double a = angle<Vector3>(aa, molecule);
    if (a <= tools::constants::quasi_linear_angle) {
      new_angles.push_back(aa);
    }
  }

  // Return list of angles
  return new_angles;
}

/// Determine paths between nodes seperated by 1 other
///
/// Finds all elements of \p distance_m that has value <= 2.
/// It then determines all connections from these two indices
///
/// \tparam Matrix
/// \param distance_m Distance matrix
/// \return List of angles
template<typename Matrix>
std::vector<Angle> all_angles(const Matrix& distance_m) {

  using boost::math::iround;

  assert(linalg::n_rows(distance_m) == linalg::n_cols(distance_m));
  const std::size_t n_rows = linalg::n_rows(distance_m);

  std::vector<Angle> angs;

  for (std::size_t j{0}; j < n_rows; j++) {
    for (std::size_t i{0}; i < j; i++) {

      if (iround(distance_m(i, j)) <= 2) {

        std::vector<Angle> A = angles(i, j, distance_m);

        for (const auto& aa : A) {
          angs.push_back(aa);
        }
      }
    }
  }

  // Return list of angles
  return angs;
}

/// Returns the angles between bonded atoms in \p molecule
///
/// \tparam Vector3
/// \tparam Matrix
/// \param distance_m Distance matrix
/// \param molecule Molecule
/// \return List of angles
template<typename Vector3, typename Matrix>
std::vector<Angle> angles(const Matrix& distance_m,
                          const molecule::Molecule<Vector3>& molecule) {
  // Return list of angles
  return valid_angles<Vector3>(all_angles(distance_m), molecule);
}

/// Determine all possible dihedrals between atoms i and j
///
/// \tparam Matrix
/// \param i First atom index
/// \param j Last atom index
/// \param distance Distance matrix
/// \return List of dihedral angles between atoms \param i and \param j
///
/// Dijkstra shortest paths algorithm returns only one shortest path. In some
/// cases however, there might be different dihedrals between the same two
/// end atoms.
template<typename Matrix>
std::vector<Dihedral>
dihedrals(std::size_t i, std::size_t j, const Matrix& distance) {

  using boost::math::iround;

  std::vector<Dihedral> dihedrals;

  const std::size_t n_atoms{linalg::n_rows(distance)};

  // Compute possible (i,k,l,j) dihedral angles
  for (std::size_t k{0}; k < n_atoms; k++) {
    if (iround(distance(k, i)) == 1 && iround(distance(k, j)) == 2) {
      for (std::size_t l{0}; l < n_atoms; l++) {
        if (iround(distance(l, i)) == 2 && iround(distance(l, j)) == 1 &&
            iround(distance(l, k)) == 1) {
          dihedrals.emplace(dihedrals.end(), i, k, l, j);
        }
      }
    }
  }

  return dihedrals;
}

/// Returns the dihedral angles between bonded atoms in \param molecule
///
/// \tparam Vector3
/// \tparam Matrix
/// \param distance_m Distance matrix
/// \param molecule Molecule
/// \return List of dihedral angles
template<typename Vector3, typename Matrix>
std::vector<Dihedral>
dihedrals(const Matrix& distance_m,
          const molecule::Molecule<Vector3>& molecule,
          const double linear_angle = tools::constants::quasi_linear_angle) {

  using boost::math::iround;

  const std::size_t n_atoms{molecule.size()};

  std::vector<Dihedral> dih;

  for (std::size_t j{0}; j < n_atoms; j++) {
    for (std::size_t i{0}; i < j; i++) {

      // A dihedral angle with terminal atoms i and j can still be present
      // when the shortest path between i and j is smaller than 3. This
      // happen when a pentagon is present (i.e. in caffeine)
      if (iround(distance_m(i, j)) <= 3) {

        const std::vector<Dihedral> D = dihedrals(i, j, distance_m);

        for (const auto& dd : D) {
          const double a1 = angle<Vector3>({dd.i, dd.j, dd.k}, molecule);
          if (a1 > linear_angle) {
            continue;
          }

          const double a2 = angle<Vector3>({dd.j, dd.k, dd.l}, molecule);
          if (a2 > linear_angle) {
            continue;
          }

          dih.push_back(dd);
        }
      }
    }
  }

  // TODO Check if enough coordinates found elsewhere
  
  // Return list of dihedral angles
  return dih;
}

template<typename Vector3, typename Matrix>
std::vector<OutOfPlaneBend>
out_of_plane_bends(const Matrix& distance_m,
                   const molecule::Molecule<Vector3>& molecule,
                   const double angle_threshold = 10.0 * tools::conversion::deg_to_rad) {

  using boost::math::iround;

  const std::size_t n_atoms{molecule.size()};

  std::vector<OutOfPlaneBend> bends;

  for (std::size_t c{0}; c < n_atoms; c++) {
    std::vector<size_t> bonded_to_c;

    for (std::size_t i{0}; i < n_atoms; i++) {
      if (iround(distance_m(i, c)) == 1) {
        bonded_to_c.push_back(i);
      }
    }
    const size_t n_bonded = bonded_to_c.size();
    if(n_bonded < 3) {
      continue;
    }

    // determine all set of 3 others with no linear angles
    for (std::size_t b_i{0}; b_i < n_bonded; b_i++) {
      const size_t i = bonded_to_c[b_i];
      for (std::size_t b_j{0}; b_j < b_i; b_j++) {
        const size_t j = bonded_to_c[b_j];

        // Check i-c-j angle not linear
        {
          const double a_icj = angle<Vector3>(Angle(i, c, j), molecule);
          if (a_icj > tools::constants::quasi_linear_angle) {
            continue;
          }
        }
        for (std::size_t b_k{0}; b_k < b_j; b_k++) {
          const size_t k = bonded_to_c[b_k];
          // Check i-c-k angle not linear
          {
            const double a_ick = angle<Vector3>(Angle(i, c, k), molecule);
            if (a_ick > tools::constants::quasi_linear_angle) {
              continue;
            }
          }

          // Check j-c-k angle not linear
          {
            const double a_jck = angle<Vector3>(Angle(j, c, k), molecule);
            if (a_jck > tools::constants::quasi_linear_angle) {
              continue;
            }
          }

          const OutOfPlaneBend bend = OutOfPlaneBend(c, i, j, k);
          const double angle = connectivity::out_of_plane_angle(bend, molecule);
          if(std::abs(angle) > angle_threshold) {
            continue;
          }
          // No linear angles in the out-of-plane bend
          bends.push_back(bend);
        }
      }
    }
  }

  return bends;
}

/// \brief Determines where x, y or z is most orthogonal to direction \p d.
template<typename Vector3>
inline Vector3 non_parallel_direction(const Vector3& d) {
  const std::vector<Vector3> directions = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  return *std::min_element(directions.begin(),
                           directions.end(),
                           [d](const Vector3& a, const Vector3& b) {
                             return linalg::dot(d, a) * linalg::dot(d, a) <
                                    linalg::dot(d, b) * linalg::dot(d, b);
                           });
}

/// \brief Determines where x, y or z is most orthogonal to angle \p a.
template<typename Vector3>
inline Vector3
non_parallel_direction(const Angle& a,
                       const molecule::Molecule<Vector3>& molecule) {
  const Vector3 d = molecule[a.k].position - molecule[a.i].position;
  return non_parallel_direction(d);
}

/// \brief Returns two vectors that are orthogonal to \p d and each other
///
/// The \p axis is required to form the first orthogonal vector. It must not
/// be parallel to \p d.
template<typename Vector3>
inline std::pair<Vector3, Vector3> orthogonal_axis(const Vector3& d,
                                                   const Vector3& axis) {

  const Vector3 first = linalg::normalize(linalg::cross(d, axis));
  const Vector3 second = linalg::normalize(linalg::cross(d, first));

  return {first, second};
}

/// \brief Returns two vectors that are orthogonal to the Angle \p a and each
/// other
///
/// The \p axis is required to form the first orthogonal vector. It must not
/// be parallel to \p d.
template<typename Vector3>
inline std::pair<Vector3, Vector3>
orthogonal_axis(const Angle& a,
                const molecule::Molecule<Vector3>& molecule,
                const Vector3& axis) {
  const Vector3 d = molecule[a.k].position - molecule[a.i].position;
  return orthogonal_axis(d, axis);
}

template<typename Vector3>
std::vector<LinearAngle<Vector3>>
valid_linear_angles(const std::vector<Angle>& angles,
                    const molecule::Molecule<Vector3>& molecule) {
  std::vector<LinearAngle<Vector3>> new_linear_angles;

  for (const auto& aa : angles) {
    const double a = angle<Vector3>(aa, molecule);
    if (a > tools::constants::quasi_linear_angle) {
      Vector3 direction = non_parallel_direction(aa, molecule);
      std::pair<Vector3, Vector3> axis =
          orthogonal_axis(aa, molecule, direction);
      new_linear_angles.push_back(LinearAngle<Vector3>{
          aa.i, aa.j, aa.k, axis.first, LinearAngleTag::First, aa.constraint});
      new_linear_angles.push_back(LinearAngle<Vector3>{aa.i,
                                                       aa.j,
                                                       aa.k,
                                                       axis.second,
                                                       LinearAngleTag::Second,
                                                       aa.constraint});
    }
  }

  // Return list of angles
  return new_linear_angles;
}

/// \brief Constructs all linear angles in the \p molecule.
///
/// An angle is consider linear if it is greater than
/// tools::constants::quasi_linear_angle. For each linear angle two instances
/// of a LinearAngle are formed to allow bending in any direction.
template<typename Vector3, typename Matrix>
std::vector<LinearAngle<Vector3>>
linear_angles(const Matrix& distance_m,
              const molecule::Molecule<Vector3>& molecule) {

  return valid_linear_angles(all_angles(distance_m), molecule);
}

// TODO: Move to transformation? (Circular dependency?)
/// Transform cartesian coordinates to internal redundant coordinates using
/// information contained in the lists of bonds, angles and dihedrals
///
/// \tparam Vector3
/// \tparam Vector
/// \param x_c Cartesian coordinates
/// \param bonds List of bonds
/// \param angles List of angles
/// \param dihedrals List of dihedral angles
/// \return
template<typename Vector3, typename Vector>
Vector cartesian_to_irc(
    const Vector& x_c,
    const std::vector<connectivity::Bond>& bonds,
    const std::vector<connectivity::Angle>& angles,
    const std::vector<connectivity::Dihedral>& dihedrals,
    const std::vector<connectivity::LinearAngle<Vector3>>& linear_angles,
    const std::vector<connectivity::OutOfPlaneBend>& out_of_plane_bends) {

  const auto n_bonds = bonds.size();
  const auto n_angles = angles.size();
  const auto n_dihedrals = dihedrals.size();
  const auto n_linear_angles = linear_angles.size();
  const auto n_bends = out_of_plane_bends.size();

  // Number of internal redundant coordinates
  const auto n_irc = n_bonds + n_angles + n_dihedrals + n_linear_angles + n_bends;

  // Allocate vector for internal redundant coordinates
  Vector q_irc{linalg::zeros<Vector>(n_irc)};

  // Offset
  std::size_t offset{0};

  // Compute bonds
  for (std::size_t i{0}; i < n_bonds; i++) {
    q_irc(i) = bond<Vector3, Vector>(bonds[i], x_c);
  }

  // Compute angles
  offset = n_bonds;
  for (std::size_t i{0}; i < n_angles; i++) {
    q_irc(i + offset) = angle<Vector3, Vector>(angles[i], x_c);
  }

  // Compute dihedrals
  offset = n_bonds + n_angles;
  for (std::size_t i{0}; i < n_dihedrals; i++) {
    q_irc(i + offset) = dihedral<Vector3, Vector>(dihedrals[i], x_c);
  }

  // Compute linear angles
  offset = n_bonds + n_angles + n_dihedrals;
  for (std::size_t i{0}; i < n_linear_angles; i++) {
    q_irc(i + offset) = angle<Vector3, Vector>(linear_angles[i], x_c);
  }

  // Compute out of plane bends
  offset = n_bonds + n_angles + n_dihedrals + n_linear_angles;
  for (std::size_t i{0}; i < n_bends; i++) {
    q_irc(i + offset) = out_of_plane_angle<Vector3, Vector>(out_of_plane_bends[i], x_c);
  }

  // Return internal redundant coordinates
  return q_irc;
}

} // namespace connectivity

} // namespace irc

#endif // IRC_CONNECTIVITY_H
