#ifndef IRC_TRANSFORMATION_H
#define IRC_TRANSFORMATION_H

#include "../connectivity/connectivity.h"
#include "../linear_algebra/linalg.h"
#include "../connectivity/wilson.h"
#include "../tools/math.h"

#include<iostream>

namespace irc {

namespace transformation {

/// Get current internal redundant coordinates for list of bonds, angles and
/// dihedral angles
///
/// \tparam Vector3
/// \tparam Vector
/// \param bonds List of bonds
/// \param angles List of angles
/// \param dihedrals List of dihedrals
/// \return Current internal redundant coordinates
template<typename Vector3, typename Vector>
Vector irc_from_bad(
    const Vector& x_cartesian,
    const std::vector<connectivity::Bond>& bonds,
    const std::vector<connectivity::Angle>& angles,
    const std::vector<connectivity::Dihedral>& dihedrals) {

  // Get number of bonds, angles and dihedrals
  size_t n_bonds{bonds.size()};
  size_t n_angles{angles.size()};
  size_t n_dihedrals{dihedrals.size()};

  // Compute number of internal redundant coordinates
  size_t n_irc{n_bonds + n_angles + n_dihedrals};

  // Allocate vector for internal redundant coordinates
  Vector q_irc{linalg::zeros<Vector>(n_irc)};

  // Offset
  size_t offset{0};

  // Compute bonds
  for (size_t i{0}; i < n_bonds; i++) {
    q_irc(i) = connectivity::bond<Vector3, Vector>(bonds[i], x_cartesian);
  }

  // Compute angles
  offset = n_bonds;
  for (size_t i{0}; i < n_angles; i++) {
    q_irc(i + offset) =
        connectivity::angle<Vector3, Vector>(angles[i], x_cartesian);
  }

  // Compute dihedrals
  offset = n_bonds + n_angles;
  for (size_t i{0}; i < n_dihedrals; i++) {
    q_irc(i + offset) =
        connectivity::dihedral<Vector3, Vector>(dihedrals[i], x_cartesian);
  }

  // Return internal redundant coordinates
  return q_irc;
}

/// Compute the root mean square value of \param v
///
/// \tparam Vector Vector
/// \param v Vector
/// \return Root mean square value of \param v
template<typename Vector>
double rms(const Vector& v) {
  size_t size{linalg::size<Vector>(v)};
  
  double sum{0};
  
  // TODO: Use linalg::norm instead
  for (size_t i{0}; i < size; i++) {
    sum += v(i) * v(i);
  }
  
  return std::sqrt(sum / size);
}

/// Transform gradient from cartesian to internal redundant coordinates,
/// without projection
///
/// \tparam Vector
/// \tparam Matrix
/// \param grad_c Gradient in cartesian coordinates
/// \param B Wilson \f$\mathbf{B}\f$ matrix
/// \param iG Pseudo inverse of the \f$\mathbf{G}\f$ matrix
/// \return Gradient in internal redundant coordinates
template<typename Vector, typename Matrix>
Vector gradient_cartesian_to_irc(const Vector &grad_c,
                                 const Matrix &B, const Matrix &iG) {
  return iG * B * grad_c;
}

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
Vector cartesian_to_irc(const Vector &x_c,
                        const std::vector<connectivity::Bond> &bonds,
                        const std::vector<connectivity::Angle> &angles,
                        const std::vector<connectivity::Dihedral> &dihedrals) {
  
  // Compute the number of internal redundant coordinates
  size_t n_irc{bonds.size() + angles.size() + dihedrals.size()};
  
  // Allocate internal redundant coordinates
  Vector q_irc{linalg::zeros<Vector>(n_irc)};
  
  // Temporary indices
  size_t idx1, idx2, idx3, idx4;
  
  // Temporary positions
  Vector3 p1, p2, p3, p4;
  
  // Offset for internal coordinates vector
  size_t offset{0};
  
  // Compute bonds
  for (size_t i{0}; i < bonds.size(); i++) {
    idx1 = bonds[i].i;
    idx2 = bonds[i].j;
    
    p1 = {x_c(3 * idx1), x_c(3 * idx1 + 1), x_c(3 * idx1 + 2)};
    p2 = {x_c(3 * idx2), x_c(3 * idx2 + 1), x_c(3 * idx2 + 2)};
    
    q_irc(i + offset) = connectivity::distance(p1, p2);
  }

  // Compute angles
  offset = bonds.size();
  for (size_t i{0}; i < angles.size(); i++) {
    idx1 = angles[i].i;
    idx2 = angles[i].j;
    idx3 = angles[i].k;
    
    p1 = {x_c(3 * idx1), x_c(3 * idx1 + 1), x_c(3 * idx1 + 2)};
    p2 = {x_c(3 * idx2), x_c(3 * idx2 + 1), x_c(3 * idx2 + 2)};
    p3 = {x_c(3 * idx3), x_c(3 * idx3 + 1), x_c(3 * idx3 + 2)};
    
    q_irc(i + offset) = connectivity::angle(p1, p2, p3);
  }

  // Compute dihedrals
  offset = bonds.size() + angles.size();
  for (size_t i{0}; i < dihedrals.size(); i++) {
    idx1 = dihedrals[i].i;
    idx2 = dihedrals[i].j;
    idx3 = dihedrals[i].k;
    idx4 = dihedrals[i].l;
    
    p1 = {x_c(3 * idx1), x_c(3 * idx1 + 1), x_c(3 * idx1 + 2)};
    p2 = {x_c(3 * idx2), x_c(3 * idx2 + 1), x_c(3 * idx2 + 2)};
    p3 = {x_c(3 * idx3), x_c(3 * idx3 + 1), x_c(3 * idx3 + 2)};
    p4 = {x_c(3 * idx4), x_c(3 * idx4 + 1), x_c(3 * idx4 + 2)};
    
    q_irc(i + offset) = connectivity::dihedral(p1, p2, p3, p4);
  }
  
  return q_irc;
}

/// Transform internal redundant displacements to cartesian coordinates
///
/// \tparam Vector3
/// \tparam Vector
/// \tparam Matrix
/// \param q_irc_old Old internal reaction coordinates
/// \param dq_irc Change in internar reaction coordinates
/// \param x_c_old Old cartesian coordinates
/// \param bonds List of bonds
/// \param angles List of angles
/// \param dihedrals List of dihedral angles
/// \param B Wilson \f$\mathbf{B}\f$ matrix
/// \param iG Pseudoinverse of the \f$\mathbf{G}\f$ matrix
/// \param max_iters Maximum number of iterations
/// \param tolerance Tolerance on change in cartesian coordinates
/// \return New cartesian coordinates
///
/// Since Cartesian coordinates are rectilinear and the internal coordinates are
/// curvilinear, the transformation must be done iteratively.
template<typename Vector3, typename Vector, typename Matrix>
Vector irc_to_cartesian(const Vector &q_irc_old,
                        const Vector &dq_irc,
                        const Vector &x_c_old,
                        const std::vector<connectivity::Bond> &bonds,
                        const std::vector<connectivity::Angle> &angles,
                        const std::vector<connectivity::Dihedral> &dihedrals,
                        size_t max_iters = 25,
                        double tolerance = 1e-6) {
  // Number of internal redundant coordinates
  size_t n_irc{ bonds.size() + angles.size() + dihedrals.size() };
  
  // Convergence flag
  bool converged{false};
  
  // Cartesian coordinates
  Vector x_c{x_c_old};
  
  // Store change in internal redundant coordinates
  Vector dq{dq_irc};
  
  // Old internal coordinates
  Vector q_old{q_irc_old};
  
  // New internal coordinates
  Vector q_new{q_irc_old};
  
  // Change in cartesian coordinates
  Vector dx{ linalg::zeros<Vector>(linalg::size(x_c_old)) };

  // Compute Wilson's B matrix
  Matrix B{
      wilson::wilson_matrix<Vector3,Vector,Matrix>(x_c,
                                                   bonds,
                                                   angles,
                                                   dihedrals)
  };

  // Compute G matrices
  Matrix G, iG;
  std::tie(G, iG) = wilson::G_matrices(B);
  
  // Compute the transpose of B
  Matrix Bt{ linalg::transpose(B) };
  
  double RMS{0};
  
  // Start iterative search
  for (size_t i{0}; i < max_iters; i++) {
    
    // Compute displacement in cartesian coordinates
    dx = Bt * iG * dq;
    
    // Check for convergence
    RMS = rms<Vector>(dx);
    if ( RMS < tolerance) {
      converged = true;
      break;
    }
    
    // Update cartesian coordinates
    x_c += dx;

    // Update Wilson B matrix
    B = wilson::wilson_matrix<Vector3,Vector,Matrix>(x_c,
                                                     bonds,
                                                     angles,
                                                     dihedrals);
  
    // Update transpose of the Wilson B matrix
    Bt = linalg::transpose(B);

    // TODO: Only one matrix needed
    // Update G matrices
    std::tie(G, iG) = wilson::G_matrices(B);
    
    // Store old internal coordinates
    q_old = q_new;
    
    // TODO: Check change in angles and dihedrals
    // Compute new internal coordinates
    q_new = cartesian_to_irc<Vector3, Vector>(x_c, bonds, angles, dihedrals);

    // Check change in dihedral angles
    size_t offset{ bonds.size() + angles.size() };
    for(size_t i{offset}; i < n_irc; i++){
      // Restrain dihedral angle on the interval [-180,180]
      q_new(i) = tools::math::pirange_deg(q_new(i));
    }
    
    // New difference in internal coordinates
    dq = dq - (q_new - q_old);
  }
  
  // If iteration does not converge, use first estimate
  if (!converged) {
    // Re-compute original B matrix
    B = wilson::wilson_matrix<Vector3,Vector,Matrix>(x_c,
                                                     bonds,
                                                     angles,
                                                     dihedrals);

    // Re compute original G matrices
    std::tie(G, iG) = wilson::G_matrices(B);

    // Compute first estimate
    x_c = x_c_old + linalg::transpose(B) * iG * dq_irc;
    
    // TODO: Something better?
    std::cerr << "WARNING: IRC_TO_CARTESIAN not converged." << std::endl;
  }
  
  return x_c;
}

} // namespace transformation

} // namespace irc

#endif //IRC_TRANSFORMATION_H
