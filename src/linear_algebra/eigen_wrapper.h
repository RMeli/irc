#ifndef IRC_EIGEN_WRAPPER_H
#define IRC_EIGEN_WRAPPER_H

#ifdef HAVE_EIGEN

#include <utility>

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/LU>

namespace linalg{

template<typename T>
size_t size(const T& a){
  return a.size();
}

template <typename T>
size_t n_rows(const T& a){
  return a.rows();
}

template <typename T>
size_t n_cols(const T& a){
  return a.cols();
}

template<typename T>
double norm(const T& a){
  return a.norm();
}

template<typename T>
double dot(const T& a, const T& b){
  return a.dot(b);
}

template <typename Vector3>
Vector3 cross(const Vector3& a, const Vector3& b){
  return a.cross(b);
}

template <typename Vector>
Vector zeros(size_t nelements){
  return Vector::Zero(nelements);
}

template <typename Matrix>
Matrix zeros(size_t nrows, size_t ncols){
  return Matrix::Zero(nrows, ncols);
}

template <typename Matrix>
Matrix identity(size_t n){
  return Matrix::Identity(n, n);
}

template <typename Matrix>
Matrix transpose(const Matrix& mat){
  return mat.transpose();
}

template <typename Matrix>
Matrix inv(const Matrix& mat){
  return mat.inverse();
}

template <typename Matrix>
Matrix pseudo_inverse(const Matrix& mat){
  return mat.pseudoInverse();
}

} // namespace linalg

#endif

#endif //IRC_EIGEN_WRAPPER_H
