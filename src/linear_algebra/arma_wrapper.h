#ifndef IRC_ARMA_H
#define IRC_ARMA_H

#ifdef HAVE_ARMA

#include <utility>

#include <armadillo>

namespace irc {

namespace linalg {

template<typename T>
size_t size(const T &a) {
  return a.size();
}

template<typename T>
size_t n_rows(const T &a) {
  return a.n_rows;
}

template<typename T>
size_t n_cols(const T &a) {
  return a.n_cols;
}

template<typename T>
double norm(const T &a) {
  return arma::norm(a);
}

template<typename T>
double dot(const T &a, const T &b) {
  return arma::dot(a, b);
}

template<typename Vector3>
Vector3 cross(const Vector3 &a, const Vector3 &b) {
  return arma::cross(a, b);
}

template<typename Vector>
Vector zeros(size_t nelements) {
  return arma::zeros<Vector>(nelements);
}

template<typename Matrix>
Matrix zeros(size_t nrows, size_t ncols) {
  return arma::zeros<Matrix>(nrows, ncols);
}

template<typename Matrix>
Matrix ones(size_t nrows, size_t ncols){
  return arma::ones<Matrix>(nrows, ncols);
}

template<typename Matrix>
Matrix identity(size_t n) {
  return arma::eye(n, n);
}

template<typename Matrix>
Matrix transpose(const Matrix &mat) {
  return arma::trans(mat);
}

template<typename Matrix>
Matrix inv(const Matrix &mat) {
  return arma::inv(mat);
}

template<typename Matrix>
Matrix pseudo_inverse(const Matrix &mat) {
  return arma::pinv(mat);
}

} // namespace linalg

} // namespace irc

#endif

#endif //IRC_ARMA_H
