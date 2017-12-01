#undef NDEBUG

#include "periodic_table.h"

#include <cassert>
#include <iomanip>
#include <iostream>

int main(){
  using namespace std;
  using namespace periodic_table;
  
  bool valid{false};
  for(size_t i{0}; i < pt_size + 1; i++){
    valid = valid_atomic_number(i);
    
    if(valid){
      cout << std::setw(3) << i;
      cout << " is a valid atomic number corresponding to element ";
      cout << pt_symbols[i] << ".\n";
      
      assert( i < pt_size );
    }
    else{
      cout << std::setw(3) << i;
      cout << " is not a valid atomic number.\n";
      
      assert( i == 0 or i == pt_size );
    }
  }
  
  for(size_t i{0}; i < pt_size; i++){
    assert( atomic_number(pt_symbols[i]) == i );
  }
  
  return 0;
}