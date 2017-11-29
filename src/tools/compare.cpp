#include "compare.h"

#include <cmath>

namespace tools {

bool nearly_zero(double a, double epsilon){
  return std::abs(a) < epsilon;
}

bool nearly_equal(double a, double b, double epsilon){
  return nearly_zero(a - b, epsilon);
}

}