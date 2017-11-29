#undef NDEBUG

#include "compare.h"

#include <cassert>

int main(){
  using namespace tools;
  
  assert( nearly_zero(1.2e-14) );
  
  assert( nearly_zero(1.2e-7, 1e-6) );
  
  assert( nearly_equal(1+2.3e-13, 1-2.3e-13) );
  
  assert( nearly_equal(1+2.3e-4, 1-2.3e-4, 1e-3) );
  
  return 0;
}