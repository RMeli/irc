#include "compare.h"

#include <cmath>

namespace tools {

namespace compare {

bool nearly_zero(double a, double epsilon) {
  return std::abs(a) < epsilon;
}

bool nearly_equal(double a, double b, double epsilon) {
  return nearly_zero(a - b, epsilon);
}

bool nearly_zero(arma::mat a, double epsilon) {
  bool nz{true};
  for (const auto &element : a) {
    if (!nearly_zero(element, epsilon)) {
      nz = false;
      break;
    }
  }
  
  return nz;
}

bool nearly_equal(arma::mat a, arma::mat b, double epsilon) {
  return nearly_zero(a - b, epsilon);
}

}

}