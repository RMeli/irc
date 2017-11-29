#include "periodic_table.h"

#include <iomanip>

namespace periodic_table{

bool valid_atomic_number(size_t an){
  return an > 0 && an < pt_size;
}

}