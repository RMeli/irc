#ifndef IRC_COMPARE_H
#define IRC_COMPARE_H

namespace tools{

bool nearly_zero(double a, double epsilon = 1e-12);

bool nearly_equal(double a, double b, double epsilon = 1e-12);

}

#endif //IRC_COMPARE_H
