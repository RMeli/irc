#include "catch.hpp"

#include "libirc/mathtools.h"

#include <cmath> // nextafter

#ifdef HAVE_ARMA
#include <armadillo>
using vec3 = arma::vec3;
using vec = arma::vec;
using mat = arma::mat;
template<typename T>
using Mat = arma::Mat<T>;
#elif HAVE_EIGEN3
#include <eigen3/Eigen/Dense>
using vec3 = Eigen::Vector3d;
using vec = Eigen::VectorXd;
using mat = Eigen::MatrixXd;
template<typename T>
using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
#else
#error
#endif

using namespace irc;

using tools::constants::pi;
using namespace tools::math;

TEST_CASE("pirange_rad produces angles in the interval (-pi,pi]") {

  const auto inclusive_angle = pi;
  CAPTURE(inclusive_angle);

  CHECK(pirange_rad(-pi) == Approx(inclusive_angle));
  CHECK(pirange_rad(0) == 0);
  CHECK(pirange_rad(pi) == Approx(inclusive_angle));
  CHECK(std::nextafter(pi, 0) < pi);
  CHECK(pirange_rad(std::nextafter(pi, 0)) == Approx(pi));
  CHECK(pirange_rad(std::nextafter(pi, 2 * pi)) == Approx(-pi));
  CHECK(pirange_rad(std::nextafter(-pi, -2 * pi)) == Approx(pi));
  CHECK(pirange_rad(std::nextafter(-pi, 0)) == Approx(-pi));

  CHECK(pirange_rad(2 * pi) == Approx(0));
  CHECK(pirange_rad(4 * pi) == Approx(0));
  CHECK(pirange_rad(6 * pi) == Approx(0));
  CHECK(pirange_rad(8 * pi) == Approx(0));
  CHECK(pirange_rad(1 * pi) == Approx(inclusive_angle));
  CHECK(pirange_rad(3 * pi) == Approx(inclusive_angle));
  CHECK(pirange_rad(5 * pi) == Approx(inclusive_angle));
  CHECK(pirange_rad(-1 * pi) == Approx(inclusive_angle));
  CHECK(pirange_rad(-3 * pi) == Approx(inclusive_angle));
  CHECK(pirange_rad(-5 * pi) == Approx(inclusive_angle));

  CHECK(pirange_rad(pi / 2.) == Approx(pi / 2.));
  CHECK(pirange_rad(-pi / 2.) == Approx(-pi / 2.));
  CHECK(pirange_rad(pi / 2. + pi) == Approx(-pi / 2.));
  CHECK(pirange_rad(-pi / 2. - pi) == Approx(pi / 2.));
}

TEST_CASE("pirange_deg produces angles in the interval (-180,180]") {

  const auto inclusive_angle = 180.0;
  CAPTURE(inclusive_angle);

  CHECK(pirange_deg(-180) == Approx(inclusive_angle));
  CHECK(pirange_deg(0) == 0);
  CHECK(pirange_deg(180) == Approx(inclusive_angle));
  CHECK(pirange_deg(std::nextafter(180, 0)) == Approx(180));
  CHECK(pirange_deg(std::nextafter(180, 2 * 180)) == Approx(-180));
  CHECK(pirange_deg(std::nextafter(-180, -2 * 180)) == Approx(180));
  CHECK(pirange_deg(std::nextafter(-180, 0)) == Approx(-180));

  CHECK(pirange_deg(2 * 180) == Approx(0));
  CHECK(pirange_deg(4 * 180) == Approx(0));
  CHECK(pirange_deg(6 * 180) == Approx(0));
  CHECK(pirange_deg(8 * 180) == Approx(0));
  CHECK(pirange_deg(1 * 180) == Approx(inclusive_angle));
  CHECK(pirange_deg(3 * 180) == Approx(inclusive_angle));
  CHECK(pirange_deg(5 * 180) == Approx(inclusive_angle));
  CHECK(pirange_deg(-1 * 180) == Approx(inclusive_angle));
  CHECK(pirange_deg(-3 * 180) == Approx(inclusive_angle));
  CHECK(pirange_deg(-5 * 180) == Approx(inclusive_angle));

  CHECK(pirange_deg(90) == Approx(90));
  CHECK(pirange_deg(-90) == Approx(-90));
  CHECK(pirange_deg(270) == Approx(-90));
  CHECK(pirange_deg(-270) == Approx(90));
}

TEST_CASE("collinear vectors") {

  vec3 v1 = {1., 0., 0.};
  vec3 v2 = {2., 0., 0.};
  vec3 v3 = {-0.5, 0., 0.};
  vec3 v4 = {1., 1., 1.};

  CHECK(collinear(v1, v2));
  CHECK(collinear(v1, v3));
  CHECK(collinear(v2, v3));

  CHECK(!collinear(v1, v4));
  CHECK(!collinear(v2, v4));
  CHECK(!collinear(v3, v4));
}