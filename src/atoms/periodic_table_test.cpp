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
    assert( i == atomic_number(symbol) );
  }
  cout << endl;
  
  return 0;
}