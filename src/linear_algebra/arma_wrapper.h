#ifndef IRC_ARMA_H
#define IRC_ARMA_H

#include <utility>

#include <armadillo>

namespace linalg{

#ifdef HAVE_ARMA

template<typename T>
size_t size(const T& a){
  return a.size();
}

template <typename T>
size_t n_rows(const T& a){
  return a.n_rows;
}

template <typename T>
size_t n_cols(const T& a){
  return a.n_cols;
}

template<typename T>
double norm(const T& a){
  return arma::norm(a);
}

template<typename T>
double dot(const T& a, const T& b){
  return arma::dot(a, b);
}

template <typename Matrix>
Matrix zeros(size_t nrows, size_t ncols){
  return std::move(arma::zeros<Matrix>(nrows, ncols));
}

#endif

}

#endif //IRC_ARMA_H
