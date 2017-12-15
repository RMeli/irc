#ifndef IRC_TRANSFORMATION_H
#define IRC_TRANSFORMATION_H

#include "../linear_algebra/linalg.h"

namespace transformation{

template <typename Vector, typename Matrix>
Vector gradient_cartesian_to_irc(const Vector& grad_irc,
                                 const Matrix& B, const Matrix& iG ){
  return iG * B * grad;
}

template <typename Vector, typename Matrix>
Vector irc_to_cartesian(const Vector& q_irc_new,
                        const Vector& q_irc_old,
                        const Vector& x_cartesian_old,
                        const Matrix& B, const Matrix& iG,
                        size_t max_iters = 25){
  bool converged{false};
  
  Vector dq_irc{q_irc_new - q_irc_old};
  
  Vector x_cartesian{ x_cartesian_old };
  
  for(size_t i{0}; i < max_iters; i ++){
    
    // TODO: Check angles and dihedrals
  }
  
  // If iteration does not converge, use first estimate
  if( !converged ){
    x_cartesian = x_cartesian_old + linalg::transpose(B) * iG * dq_irc;
  }
  
  return x_cartesian;
}

template <typename Vector>
double rms(const Vector& v){
  size_t size{linalg::size(v)};
  
  double sum{0};
  
  // TODO: Use linalg::norm instead
  for(size_t i{0}; i < size; i++){
    sum += v(i)*v(i);
  }
  
  return std::sqrt( sum / size );
}

}

#endif //IRC_TRANSFORMATION_H
