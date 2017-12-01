#undef NDEBUG

#include "atom.h"

#include "periodic_table.h"
#include "../tools/comparison.h"

#include <cassert>
#include <iostream>

int main(){
  using namespace std;
  using namespace atom;
  
  for(size_t i{1}; i < periodic_table::pt_size; i++){
    assert( periodic_table::valid_atomic_number(i) );
    
    Atom a{i,{0.,1.,2.}};
    
    assert( symbol(a.atomic_number) ==
            periodic_table::pt_symbols[i] );
    
    assert( mass(a.atomic_number) ==
            periodic_table::pt_masses[i] );
    
    assert( covalent_radius(a.atomic_number) ==
            periodic_table::pt_covalent_radii[i]);
    
    cout << a << endl;
  }
  
  for(size_t i{1}; i < periodic_table::pt_size; i++){
    assert( periodic_table::valid_atomic_number(i) );
    
    Atom a{periodic_table::pt_symbols[i],{0.,1.,2.}};
    
    assert( symbol(a.atomic_number) ==
            periodic_table::pt_symbols[i] );
    
    assert( mass(a.atomic_number) ==
            periodic_table::pt_masses[i] );
    
    assert( covalent_radius(a.atomic_number) ==
            periodic_table::pt_covalent_radii[i]);
    
    cout << a << endl;
  }

  return 0;
}