#ifndef IRC_IRC_H
#define IRC_IRC_H

#include "../atoms/molecule.h"
#include "../connectivity/connectivity.h"
#include "../connectivity/wilson.h"

#include <utility>

namespace irc {

template <typename Vector3, typename Vector, typename Matrix>
class IRC {
 public:
  template <typename Vector3>
  IRC(const Molecule<Vector3>& molecule);
  
  Matrix projected_initial_hessian();
  Vector grad_cartesian_to_projected_irc(const Vector& grad_c);
  Vector irc_to_cartesian(const Vector& q_irc);
  
 private:
  /// Distance matrix
  ///
  /// Contains information about connectivity and shortest paths within
  /// the initial molecule
  Matrix distance_m;
  
  /// Predecessor matrix
  ///
  /// Contains information to build the shortest path between atoms
  Matirx predecessors_m;
  
  /// Wilson B matrix
  Matrix B;
  
  /// G matrix
  Matrix G;
  
  /// Generalized inverse of the G matrix
  Matrix iG;
  
  /// Projector
  Matrix P;
};

template <typename Vector3, typename Vector, typename Matrix>
IRC::IRC<Vector3, Vector, Matrix>(const Molecule<Vector3>& molecule){
  // Compute interatomic distances
  Matrix dd{ connectivity::distances<Vector3, Matrix>(molecule) };
  
  // Compute adjacency matrix (graph)
  connectivity::UGraph adj{ connectivity::adjacency_matrix(dd, molecule) };
  
  // Compute distance matrix and predecessor matrix
  Matrix distance_m, predecessors_m;
  std::tie(distance_m, predecessors_m) =
      connectivity::distance_matrix<Matrix>(adj) ;
  
  // Compute bonds
  std::vector<connectivity::Bond<Vector3>> bonds{
      connectivity::bonds(dist, molecule)};
  
  // Compute angles
  std::vector<connectivity::Angle<Vector3>> angles{
      connectivity::angles(dist, predecessors, molecule)};
  
  // Compute dihedrals
  std::vector<connectivity::Dihedral<Vector3>> dihedrals{
      connectivity::dihedrals(dist, predecessors, molecule)};
  
  // Store initial Wilson's B matrix
  B = wilson::wilson_matrix<Vector3, Matrix>(molecule.size(),
                                             bonds, angles, dihedrals);
  
  // Compute G and iG (generalized inverse)
  std::tie(G, iG) = wilson::G_matirces(B);
  
  // Compute projector P
  P = wilson::projector(G, iG);
}

}


#endif //IRC_IRC_H
