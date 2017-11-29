#undef NDEBUG

#include "atom.h"

#include <cassert>
#include <iostream>

int main(){
  using namespace std;
  using namespace atoms;
  
  for(size_t i{0}; i < pt_size; i++){
    Atom a{pt_symbols[i],{0.,1.,2.}};
    
    cout << a << endl;
  
    assert( a.get_atomic_number() == i);
    assert( a.get_symbol() == pt_symbols[i] );
    assert( a.get_mass() == pt_masses[i] );
    assert( a.get_covalent_radius() == pt_covalent_radii[i] );
  }
  
  return 0;
}