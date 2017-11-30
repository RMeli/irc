#ifndef IRC_COMPARE_H
#define IRC_COMPARE_H

#include <armadillo>

namespace tools{

namespace comparison {

bool nearly_zero(double a, double epsilon = 1e-12);

bool nearly_equal(double a, double b, double epsilon = 1e-12);

bool nearly_zero(arma::mat a, double epsilon = 1e-12);

bool nearly_equal(arma::mat a, arma::mat b, double epsilon = 1e-12);

}

}

#endif //IRC_COMPARE_H
