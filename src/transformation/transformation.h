#ifndef IRC_TRANSFORMATION_H
#define IRC_TRANSFORMATION_H

#include "../connectivity/connectivity.h"
#include "../linear_algebra/linalg.h"

#include<iostream>

namespace transformation{


/// Compute the root mean square value of \param v
///
/// \tparam Vector Vector
/// \param v Vector
/// \return Root mean square value of \param v
template <typename Vector>
double rms(const Vector& v){
  size_t size{linalg::size<Vector>(v)};
  
  double sum{0};
  
  // TODO: Use linalg::norm instead
  for(size_t i{0}; i < size; i++){
    sum += v(i)*v(i);
  }
  
  return std::sqrt( sum / size );
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
template <typename Vector, typename Matrix>
Vector gradient_cartesian_to_irc(const Vector& grad_c,
                                 const Matrix& B, const Matrix& iG ){
  return iG * B * grad_c;
}

template <typename Vector3, typename Vector>
Vector cartesian_to_irc(const Vector& x_c,
       const std::vector<connectivity::Bond<Vector3>>& bonds,
       const std::vector<connectivity::Angle<Vector3>>& angles,
       const std::vector<connectivity::Dihedral<Vector3>>& dihedrals){
  
  // Compute the number of internal redundant coordinates
  size_t n_irc{ bonds.size() + angles.size() + dihedrals.size() };
  
  // Allocate internal redundant coordinates
  Vector q_irc{ linalg::zeros<Vector>(n_irc) };
  
  size_t idx1, idx2, idx3, idx4;
  
  Vector3 p1, p2, p3, p4;
  
  size_t offset{0};
  
  // Compute bonds
  for(size_t i{0}; i < bonds.size(); i++){
    idx1 = bonds[i].i;
    idx2 = bonds[i].j;
    
    p1 = { x_c(3*idx1), x_c(3*idx1+1), x_c(3*idx1+2) };
    p2 = { x_c(3*idx2), x_c(3*idx2+1), x_c(3*idx2+2) };
    
    q_irc(i) = connectivity::distance(p1, p2);
  }
  
  offset = bonds.size();
  for(size_t i{0}; i < angles.size(); i++){
    idx1 = angles[i].i;
    idx2 = angles[i].j;
    idx3 = angles[i].k;
    
    p1 = { x_c(3*idx1), x_c(3*idx1+1), x_c(3*idx1+2) };
    p2 = { x_c(3*idx2), x_c(3*idx2+1), x_c(3*idx2+2) };
    p3 = { x_c(3*idx3), x_c(3*idx3+1), x_c(3*idx3+2) };
    
    q_irc(i + offset) = connectivity::angle(p1, p2, p3);
  }
  
  offset = bonds.size() + angles.size();
  for(size_t i{0}; i < dihedrals.size(); i++){
    idx1 = dihedrals[i].i;
    idx2 = dihedrals[i].j;
    idx3 = dihedrals[i].k;
    idx4 = dihedrals[i].l;
    
    p1 = { x_c(3*idx1), x_c(3*idx1+1), x_c(3*idx1+2) };
    p2 = { x_c(3*idx2), x_c(3*idx2+1), x_c(3*idx2+2) };
    p3 = { x_c(3*idx3), x_c(3*idx3+1), x_c(3*idx3+2) };
    p4 = { x_c(3*idx4), x_c(3*idx4+1), x_c(3*idx4+2) };
    
    q_irc(i + offset) = connectivity::dihedral(p1, p2, p3, p4);
  }
  
  return q_irc;
}

///
/// \tparam Vector3
/// \tparam Vector
/// \tparam Matrix
/// \param q_irc_old
/// \param dq_irc
/// \param x_c_old
/// \param bonds
/// \param angles
/// \param dihedrals
/// \param B
/// \param iG
/// \param max_iters
/// \param tolerance
/// \return
///
/// Convergence: V. Bakken and T. Helgaker, J. Chem. Phys. 117, 9160 (2002).
template <typename Vector3, typename Vector, typename Matrix>
Vector irc_to_cartesian(const Vector& q_irc_old,
           const Vector& dq_irc,
           const Vector& x_c_old,
           const std::vector<connectivity::Bond<Vector3>>& bonds,
           const std::vector<connectivity::Angle<Vector3>>& angles,
           const std::vector<connectivity::Dihedral<Vector3>>& dihedrals,
           const Matrix& B,
           const Matrix& iG,
           size_t max_iters = 25,
           double tolerance = 1e-6){
  
  // Convergence flag
  bool converged{false};
  
  // Cartesian coordinates
  Vector x_c{ x_c_old };
  
  // Store change in internal redundant coordinates
  Vector dq{ dq_irc };
  
  Vector dq1{ dq_irc };
  
  // Old internal coordinates
  Vector q_old{ q_irc_old };
  
  // New internal coordinates
  Vector q_new{ q_irc_old };
  
  // Change in cartesian coordinates
  Vector dx{ linalg::zeros<Vector>( linalg::size(x_c_old)) };
  
  // Compute the transpose of B
  Matrix Bt{ linalg::transpose(B) };
  
  // Start iterative search
  for(size_t i{0}; i < max_iters; i ++){
    
    // Compute displacement in cartesian coordinates
    dx = Bt * iG * dq;
    
    // Check for convergence
    if( rms<Vector>(dx) < tolerance ){
      converged = true;
      break;
    }
    
    // Update cartesian coordinates
    x_c += dx;
    
    // Store old internal coordinates
    q_old = q_new;
  
    // TODO: Check change in angles and dihedrals
    // Compute new internal coordinates
    q_new = cartesian_to_irc(x_c, bonds, angles, dihedrals);
    
    // New difference in internal coordinates
    dq = dq - (q_new - q_old);
  }
  
  // If iteration does not converge, use first estimate
  if( !converged ){
    x_c = x_c_old + linalg::transpose(B) * iG * dq_irc;
    
  // TODO: Something better?
    std::cerr << "WARNING: IRC_TO_CARTESIAN not converged." << std::endl;
  }
  
  return x_c;
}

}

#endif //IRC_TRANSFORMATION_H
