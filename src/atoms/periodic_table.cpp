#include "periodic_table.h"

#include <iomanip>

namespace atoms{

size_t atomic_number(const std::string& symbol){
  size_t idx{0};
  for(size_t i{0}; i < pt_size; i++){
    if(symbol == pt_symbols[i]){
      idx = i;
      break;
    }
  }
  
  return idx;
}

}