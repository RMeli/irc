#include "catch.hpp"

#include "libirc/linalg.h"

#ifdef HAVE_ARMA
#include <armadillo>
using vec3 = arma::vec3;
using vec = arma::vec;
using mat = arma::mat;
#elif HAVE_EIGEN3
#include <Eigen/Dense>
using vec3 = Eigen::Vector3d;
using vec = Eigen::VectorXd;
using mat = Eigen::MatrixXd;
#else
#error
#endif

using namespace irc;

TEST_CASE("Size", "[size]") {
  SECTION("Vector") {
    vec v = {1, 2, 3, 4, 5};

    CHECK(linalg::size(v) == 5);
  }

  SECTION("Matrix") {
    mat m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

    CHECK(linalg::size(m) == 9);
  }
}

TEST_CASE("Number of rows", "[nrows]") {
  SECTION("Vector") {
    vec v = {1, 2, 3, 4, 5};

    CHECK(linalg::n_rows(v) == 5);
  }

  SECTION("Matrix") {
    mat m = {{1, 2, 3}, {4, 5, 6}};

    CHECK(linalg::n_rows(m) == 2);
  }
}

TEST_CASE("Number of columns", "[ncols]") {
  SECTION("Vector") {
    vec v = {1, 2, 3, 4, 5};

    CHECK(linalg::n_cols(v) == 1);
  }

  SECTION("Matrix") {
    mat m = {{1, 2, 3}, {4, 5, 6}};

    CHECK(linalg::n_cols(m) == 3);
  }
}

TEST_CASE("Norm and normalization", "[norm]") {
  SECTION("Vector") {
    vec v = {1, 2, 3, 4};

    Approx t(std::sqrt(1 * 1 + 2 * 2 + 3 * 3 + 4 * 4));
    CHECK(linalg::norm(v) == t);

    v = linalg::normalize(v);

    CHECK(linalg::norm(v) == Approx(1.0));
  }

  SECTION("Matrix") {
    mat m = {{1, 2, 3}, {4, 5, 6}};

    Approx t(std::sqrt(91));
    CHECK(linalg::norm(m) == t);
  }
}

TEST_CASE("Dot product", "[dot]") {
  vec v1 = {1, 2, 3, 4};
  vec v2 = {4, 3, 2, 1};

  CHECK(linalg::dot(v1, v1) == Approx(std::pow(linalg::norm(v1), 2)));

  CHECK(linalg::dot(v2, v2) == Approx(std::pow(linalg::norm(v2), 2)));

  CHECK(linalg::dot(v1, v2) == Approx(1 * 4 + 2 * 3 + 3 * 2 + 4 * 1));

  CHECK(linalg::dot(v1, v2) == Approx(linalg::dot(v2, v1)));
}

TEST_CASE("Cross product", "[cross]") {
  vec3 v1 = {1, 2, 3};
  vec3 v2 = {3, 2, 1};

  vec3 v = linalg::cross(v1, v2);

  CHECK(v(0) == Approx(-4));
  CHECK(v(1) == Approx(8));
  CHECK(v(2) == Approx(-4));

  v = linalg::cross(v1, v1);

  CHECK(v(0) == Approx(0));
  CHECK(v(1) == Approx(0));
  CHECK(v(2) == Approx(0));

  v = linalg::cross(v2, v2);

  CHECK(v(0) == Approx(0));
  CHECK(v(1) == Approx(0));
  CHECK(v(2) == Approx(0));
}

TEST_CASE("Zeros", "[zeros]") {
  SECTION("Vector") {
    std::size_t n{100};

    vec v = linalg::zeros<vec>(n);

    REQUIRE(linalg::size(v) == n);
    for (std::size_t i{0}; i < n; i++) {
      CHECK(v(i) == Approx(0.));
    }
  }

  SECTION("Matrix") {
    std::size_t n_r{10};
    std::size_t n_c{20};

    std::size_t n{n_r * n_c};

    mat m = linalg::zeros<mat>(n_r, n_c);

    REQUIRE(linalg::size(m) == n);
    for (std::size_t i{0}; i < n; i++) {
      CHECK(m(i) == Approx(0.));
    }
  }
}

TEST_CASE("Ones", "[ones]") {
  SECTION("Matrix") {
    std::size_t n_r{10};
    std::size_t n_c{20};

    std::size_t n{n_r * n_c};

    mat m = linalg::ones<mat>(n_r, n_c);

    REQUIRE(linalg::size(m) == n);
    for (std::size_t i{0}; i < n; i++) {
      CHECK(m(i) == Approx(1.));
    }
  }
}

TEST_CASE("Identity", "[identity]") {
  SECTION("Matrix") {

    std::size_t n{10};

    mat m = linalg::identity<mat>(n);

    REQUIRE(linalg::size(m) == n * n);
    for (std::size_t i{0}; i < n; i++) {
      for (std::size_t j{0}; j < n; j++) {

        double m_ij = m(i, j);

        if (i == j) {
          CHECK(m_ij == Approx(1.));
        } else {
          CHECK(m_ij == Approx(0.));
        }
      }
    }
  }
}

TEST_CASE("Transpose", "[transpose]") {
  SECTION("Matrix") {

    mat m = {{1, 2, 3}, {4, 5, 6}};

    mat mt = linalg::transpose(m);

    std::size_t n_r{linalg::n_rows(m)};
    std::size_t n_c{linalg::n_cols(m)};

    REQUIRE(linalg::size(mt) == linalg::size(m));
    REQUIRE(linalg::n_rows(mt) == n_c);
    REQUIRE(linalg::n_cols(mt) == n_r);

    for (std::size_t i{0}; i < n_r; i++) {
      for (std::size_t j{0}; j < n_c; j++) {
        CHECK(mt(j, i) == Approx(m(i, j)));
      }
    }
  }
}

TEST_CASE("Inverse", "[inv]") {
  SECTION("diagonal") {

    std::size_t n{10};

    mat m = linalg::zeros<mat>(n, n);

    for (std::size_t i{0}; i < n; i++) {
      m(i, i) = i + 1;
    }

    mat inv = linalg::inv(m);

    for (std::size_t i{0}; i < n; i++) {
      for (std::size_t j{0}; j < n; j++) {
        if (i == j) {
          CHECK(inv(i, j) == Approx(1. / m(i, j)));
        } else {
          CHECK(inv(i, j) == Approx(0.));
        }
      }
    }
  }

  SECTION("2x2") {

    double a{1}, b{2}, c{3}, d{4};
    double det{a * d - b * c};

    mat m = {{1, 2}, {3, 4}};

    mat mi = {{d, -b}, {-c, a}};
    mi *= 1. / det;

    mat inv = linalg::inv(m);

    std::size_t n{linalg::size(m)};
    REQUIRE(linalg::size(inv) == n);

    for (std::size_t i{0}; i < n; i++) {
      CHECK(inv(i) == Approx(mi(i)));
    }
  }
}

TEST_CASE("Pseudo Inverse", "[pinv]") {
  mat m = {{1, 2, 3}, {4, 5, 6}};

  std::size_t n_r{linalg::n_rows(m)};
  std::size_t n_c{linalg::n_cols(m)};

  // WolframAlpha
  mat p = {{-17, 8}, {-2, 2}, {13, -4}};
  p /= 18.;

  mat pinv = linalg::pseudo_inverse(m);

  REQUIRE(linalg::n_rows(pinv) == n_c);
  REQUIRE(linalg::n_cols(pinv) == n_r);

  std::size_t n{n_r * n_c};

  for (std::size_t i{0}; i < n; i++) {
    CHECK(pinv(i) == Approx(p(i)));
  }
}
