#ifndef IRC_LINALG_H
#define IRC_LINALG_H

// Armadillo wrapper
#include "arma_wrapper.h"

namespace linalg{

/// Size of a given container
///
/// \tparam T Container type
/// \param a Container
/// \return Container size
template<typename T>
size_t size(const T& a);

/// Number of rows of a given matrix
///
/// \tparam T Matrix type
/// \param a Matrix
/// \return Number of rows
template <typename T>
size_t n_rows(const T& a);

/// Number of columns of a given matrix
///
/// \tparam T Matrix type
/// \param a Matrix
/// \return Number of columns
template <typename T>
size_t n_cols(const T& a);

/// Norm of a given vector or matrix
///
/// \tparam T Vector or matrix type
/// \return Norm
template <typename T>
double norm(const T&);

/// Dot product between two vectors
///
/// \tparam T Vector type
/// \param a Vector
/// \param b Vector
/// \return Dot product
template<typename T>
double dot(const T& a, const T& b);

template <typename Vector>
Vector zeros(size_t nelements);

template <typename Matrix>
Matrix zeros(size_t nrows, size_t ncols);

template <typename Matrix>
Matrix transpose(const Matrix& mat);

template <typename Matrix>
Matrix pseudo_inverse(const Matrix& mat);


}

#endif //IRC_LINALG_H_H
