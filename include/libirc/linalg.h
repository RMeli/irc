#ifndef IRC_LINALG_H
#define IRC_LINALG_H

#ifdef HAVE_ARMA
#include <armadillo>
#elif HAVE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/LU>
#else
#error
#endif

namespace irc {

namespace linalg {

/// Size of a given container
///
/// \tparam T Container type
/// \param a Container
/// \return Container size
template<typename T>
size_t size(const T &a){
  return a.size();
}

/// Number of rows of a given matrix
///
/// \tparam T Matrix type
/// \param a Matrix
/// \return Number of rows
template<typename T>
size_t n_rows(const T &a){
#ifdef HAVE_ARMA
  return a.n_rows;
#elif HAVE_EIGEN
  return a.rows();
#else
#error
#endif
}

/// Number of columns of a given matrix
///
/// \tparam T Matrix type
/// \param a Matrix
/// \return Number of columns
template<typename T>
size_t n_cols(const T &a){
#ifdef HAVE_ARMA
  return a.n_cols;
#elif HAVE_EIGEN
  return a.cols();
#else
#error
#endif
}

/// Norm of a given vector or matrix
///
/// \tparam T Vector or matrix type
/// \return Norm
template<typename T>
double norm(const T &a){
#ifdef HAVE_ARMA
  return arma::norm(a);
#elif HAVE_EIGEN
  return a.norm();
#else
#error
#endif
}

/// Dot product between two vectors
///
/// \tparam T Vector type
/// \param a Vector
/// \param b Vector
/// \return Dot product
template<typename T>
double dot(const T &a, const T &b){
#ifdef HAVE_ARMA
  return arma::dot(a,b);
#elif HAVE_EIGEN
  return a.dot(b);
#else
#error
#endif
}

template<typename Vector3>
Vector3 cross(const Vector3 &a, const Vector3 &b){
#ifdef HAVE_ARMA
  return arma::cross(a,b);
#elif HAVE_EIGEN
  return a.cross(b);
#else
#error
#endif
}

template<typename Vector>
Vector zeros(size_t nelements){
#ifdef HAVE_ARMA
  return arma::zeros<Vector>(nelements);
#elif HAVE_EIGEN
  return Vector::Zero(nelements);
#else
#error
#endif
}

template<typename Matrix>
Matrix zeros(size_t nrows, size_t ncols){
#ifdef HAVE_ARMA
  return arma::zeros<Matrix>(nrows, ncols);
#elif HAVE_EIGEN
  return Matrix::Zero(nrows, ncols);
#else
#error
#endif
}

template<typename Matrix>
Matrix ones(size_t nrows, size_t ncols){
#ifdef HAVE_ARMA
  return arma::ones<Matrix>(nrows, ncols);
#elif HAVE_EIGEN
#error
#else
#error
#endif
}

template<typename Matrix>
Matrix identity(size_t n){
#ifdef HAVE_ARMA
  return arma::eye(n, n);
#elif HAVE_EIGEN
  return Matrix::Identity(n, n);
#else
#error
#endif
}

template<typename Matrix>
Matrix transpose(const Matrix &mat){
#ifdef HAVE_ARMA
  return arma::trans(mat);
#elif HAVE_EIGEN
  return mat.transpose();
#else
#error
#endif
}

template<typename Matrix>
Matrix inv(const Matrix &mat){
#ifdef HAVE_ARMA
  return arma::inv(mat);
#elif HAVE_EIGEN
  return mat.inverse();
#else
#error
#endif
}

template<typename Matrix>
Matrix pseudo_inverse(const Matrix &mat){
#ifdef HAVE_ARMA
  return arma::pinv(mat);
#elif HAVE_EIGEN
  return mat.pseudoInverse();
#else
#error
#endif
}


} // namespace linalg

} // namespace irc

#endif //IRC_LINALG_H_H
