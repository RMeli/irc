#undef NDEBUG

#include "comparison.h"

#include <cassert>

#include <armadillo>

int main(){
  using namespace tools::comparison;
  
  assert( nearly_zero(1.2e-14) );
  
  assert( !nearly_zero(1.2e-11) );
  
  assert( nearly_zero(1.2e-7, 1e-6) );
  
  assert( !nearly_zero(1.2e-5, 1e-6) );
  
  assert( nearly_equal(1+2.3e-13, 1-2.3e-13) );
  
  assert( !nearly_equal(1, 2) );
  
  assert( nearly_equal(1+2.3e-4, 1-2.3e-4, 1e-3) );
  
  assert( !nearly_equal(1+2.3e-3, 1-2.3e-3, 1e-3) );
  
  arma::mat a1{
      {1.1e-13, 2.2e-13},
      {3.3e-14, 3.3e-15}
  };
  
  assert( nearly_zero(a1) );
  
  arma::mat a2{
      {1.1e-4, 2.2e-5},
      {3.3e-5, 3.3e-6}
  };
  
  assert( nearly_zero(a2, 1e-3) );
  
  arma::mat a3{
      {1.1e-11, 2.2e-13},
      {3.3e-14, 3.3e-15}
  };
  
  assert( !nearly_zero(a3) );
  
  arma::mat a4{
      {1.1e-3, 2.2e-5},
      {3.3e-5, 3.3e-6}
  };
  
  assert( !nearly_zero(a4, 1e-3) );
  
  arma::mat b1{
      {2.2e-13, 2.2e-13},
      {3.3e-14, 3.3e-15}
  };
  
  assert( nearly_equal(a1, b1) );
  
  arma::mat b2{
      {2.2e-4, 2.2e-5},
      {3.3e-5, 3.3e-6}
  };
  
  assert( nearly_equal(a2, b2, 1e-3) );
  
  arma::mat b3{
      {2.2e-11, 2.2e-13},
      {3.3e-14, 3.3e-15}
  };
  
  assert( !nearly_equal(a3, b3) );
  
  arma::mat b4{
      {2.2e-3, 2.2e-5},
      {3.3e-5, 3.3e-6}
  };
  
  assert( !nearly_equal(a4, b4, 1e-3) );
  
  return 0;
}