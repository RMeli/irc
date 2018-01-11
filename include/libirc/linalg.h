#ifndef IRC_LINALG_H
#define IRC_LINALG_H

// Armadillo wrapper
#include "libirc/wrapper/arma_wrapper.h"

// Eigen wrapper
#include "libirc/wrapper/eigen_wrapper.h"

namespace irc {

namespace linalg {

/// Size of a given container
///
/// \tparam T Container type
/// \param a Container
/// \return Container size
template<typename T>
size_t size(const T &a);

/// Number of rows of a given matrix
///
/// \tparam T Matrix type
/// \param a Matrix
/// \return Number of rows
template<typename T>
size_t n_rows(const T &a);

/// Number of columns of a given matrix
///
/// \tparam T Matrix type
/// \param a Matrix
/// \return Number of columns
template<typename T>
size_t n_cols(const T &a);

/// Norm of a given vector or matrix
///
/// \tparam T Vector or matrix type
/// \return Norm
template<typename T>
double norm(const T &);

/// Dot product between two vectors
///
/// \tparam T Vector type
/// \param a Vector
/// \param b Vector
/// \return Dot product
template<typename T>
double dot(const T &a, const T &b);

template<typename Vector3>
Vector3 cross(const Vector3 &a, const Vector3 &b);

template<typename Vector>
Vector zeros(size_t nelements);

template<typename Matrix>
Matrix zeros(size_t nrows, size_t ncols);

template<typename Matrix>
Matrix ones(size_t nrows, size_t ncols);

template<typename Matrix>
Matrix identity(size_t n);

template<typename Matrix>
Matrix transpose(const Matrix &mat);

template<typename Matrix>
Matrix inv(const Matrix &mat);

template<typename Matrix>
Matrix pseudo_inverse(const Matrix &mat);


} // namespace linalg

} // namespace irc

#endif //IRC_LINALG_H_H
