#undef NDEBUG

#include "periodic_table.h"

#include <cassert>
#include <iostream>
#include <string>

int main(){
  using namespace std;
  using namespace atoms;
  
  std::string symbol{""};
  for(size_t i{0}; i < pt_size; i++){
    symbol = pt_symbols[i];
    cout << '(' << symbol << ',' << atomic_number(symbol) << ") ";
    
    assert( i == atomic_number(symbol) );
  }
  cout << endl;
  
  for(size_t i{0}; i < pt_size; i++){
    print_atom_info(pt_symbols[i]);
  }
  
  return 0;
}