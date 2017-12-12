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

template <typename Vector3>
Vector3 cross(const Vector3& a, const Vector3& b){
  return arma::cross(a, b);
}

template <typename Vector>
Vector zeros(size_t nelements){
  return std::move(arma::zeros<Vector>(nelements));
}

template <typename Matrix>
Matrix zeros(size_t nrows, size_t ncols){
  return std::move(arma::zeros<Matrix>(nrows, ncols));
}

template <typename Matrix>
Matrix transpose(const Matrix& mat){
  return std::move(arma::trans(mat));
}

template <typename Matrix>
Matrix pseudo_inverse(const Matrix& mat){
  return std::move(arma::pinv(mat));
}

#endif

}

#endif //IRC_ARMA_H
