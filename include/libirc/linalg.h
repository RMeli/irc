#ifndef IRC_LINALG_H
#define IRC_LINALG_H

#ifdef HAVE_ARMA
#include <armadillo>
#elif HAVE_EIGEN3
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
/// \tparam T
/// \param a Container
/// \return Size of \param a
template<typename T>
std::size_t size(const T& a) {
  return a.size();
}

/// Number of rows of a given matrix
///
/// \tparam T
/// \param a Matrix
/// \return Number of rows of \param a
template<typename T>
std::size_t n_rows(const T& a) {
#ifdef HAVE_ARMA
  return a.n_rows;
#elif HAVE_EIGEN3
  return a.rows();
#else
#error
#endif
}

/// Number of columns of a given matrix
///
/// \tparam T
/// \param a Matrix
/// \return Number of columns of \param a
template<typename T>
std::size_t n_cols(const T& a) {
#ifdef HAVE_ARMA
  return a.n_cols;
#elif HAVE_EIGEN3
  return a.cols();
#else
#error
#endif
}

/// Norm of a given vector or matrix
///
/// \tparam T
/// \return Norm of \param a
template<typename T>
double norm(const T& a) {
#ifdef HAVE_ARMA
  return arma::norm(a);
#elif HAVE_EIGEN3
  return a.norm();
#else
#error
#endif
}

/// Dot product between two vectors
///
/// \tparam T
/// \param a Vector
/// \param b Vector
/// \return Dot product between \param a and \param b
template<typename T>
double dot(const T& a, const T& b) {
#ifdef HAVE_ARMA
  return arma::dot(a, b);
#elif HAVE_EIGEN3
  return a.dot(b);
#else
#error
#endif
}

/// Cross product between two vectors
///
/// \tparam Vector3
/// \param a Vector
/// \param b Vector
/// \return Cross product between \param a and \param b
template<typename Vector3>
Vector3 cross(const Vector3& a, const Vector3& b) {
#ifdef HAVE_ARMA
  return arma::cross(a, b);
#elif HAVE_EIGEN3
  return a.cross(b);
#else
#error
#endif
}

/// Allocate column vector of zeros
///
/// \tparam Vector
/// \param nelements Vector size
/// \return Column vector full of zeros
template<typename Vector>
Vector zeros(std::size_t nelements) {
#ifdef HAVE_ARMA
  return arma::zeros<Vector>(nelements);
#elif HAVE_EIGEN3
  return Vector::Zero(nelements);
#else
#error
#endif
}

/// Allocate matrix of zeros
///
/// \tparam Matrix
/// \param nrows Number of rows
/// \param ncols Number of columns
/// \return Matrix full of zeros
template<typename Matrix>
Matrix zeros(std::size_t nrows, std::size_t ncols) {
#ifdef HAVE_ARMA
  return arma::zeros<Matrix>(nrows, ncols);
#elif HAVE_EIGEN3
  return Matrix::Zero(nrows, ncols);
#else
#error
#endif
}

/// Allocate matrix of ones
/// \tparam Matrix
/// \param nrows Number of rows
/// \param ncols Number of columns
/// \return Matrix full of ones
template<typename Matrix>
Matrix ones(std::size_t nrows, std::size_t ncols) {
#ifdef HAVE_ARMA
  return arma::ones<Matrix>(nrows, ncols);
#elif HAVE_EIGEN3
  return Matrix::Ones(nrows, ncols);
#else
#error
#endif
}

/// Allocate identity matrix
///
/// \tparam Matrix
/// \param n Linear size of the identity matrix
/// \return Identity matrix
template<typename Matrix>
Matrix identity(std::size_t n) {
#ifdef HAVE_ARMA
  return arma::eye(n, n);
#elif HAVE_EIGEN3
  return Matrix::Identity(n, n);
#else
#error
#endif
}

/// Matrix transpose
///
/// \tparam Matrix
/// \param mat Matrix
/// \return Transpose of \param mat
template<typename Matrix>
Matrix transpose(const Matrix& mat) {
#ifdef HAVE_ARMA
  return arma::trans(mat);
#elif HAVE_EIGEN3
  return mat.transpose();
#else
#error
#endif
}

/// Inverse matrix
///
/// \tparam Matrix
/// \param mat Matrix
/// \return Inverse of \param mat
template<typename Matrix>
Matrix inv(const Matrix& mat) {
#ifdef HAVE_ARMA
  return arma::inv(mat);
#elif HAVE_EIGEN3
  return mat.inverse();
#else
#error
#endif
}

/// Pseudo-inverse matrix
///
/// \tparam Matrix
/// \param mat Matrix
/// \return Pseudo-inverse of \param mat
template<typename Matrix>
Matrix pseudo_inverse(const Matrix& mat) {
#ifdef HAVE_ARMA
  return arma::pinv(mat);
#elif HAVE_EIGEN3
  return mat.completeOrthogonalDecomposition().pseudoInverse();
#else
#error
#endif
}

} // namespace linalg

} // namespace irc

#endif // IRC_LINALG_H_H
