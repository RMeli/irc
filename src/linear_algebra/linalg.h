#ifndef IRC_LINALG_H
#define IRC_LINALG_H

// Armadillo wrapper
#include "arma_wrapper.h"

namespace linalg{

template<typename T>
size_t size(const T& a);

template <typename T>
size_t n_rows(const T& a);

template <typename T>
size_t n_cols(const T& a);

template <typename T>
double norm(const T&);

template<typename T>
double dot(const T& a, const T& b);

template <typename Vector>
Vector zeros(size_t nelements);

template <typename Matrix>
Matrix zeros(size_t nrows, size_t ncols);


}

#endif //IRC_LINALG_H_H
