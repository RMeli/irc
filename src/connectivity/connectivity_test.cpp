#undef NDEBUG

#include "connectivity.h"
#include "../atoms/periodic_table.h"

#include <cassert>
#include <iostream>
#include <vector>

#include <armadillo>

int main(){
  
  using namespace std;
  
  size_t n_atoms{4};
  
  std::vector<std::string> CH2O_symbols = {"C", "O", "H", "H"};
  
  arma::mat CH2O_positions = {
      {0.000000,  0.000000,  -0.537500},
      {0.000000,  0.000000,   0.662500},
      {0.000000,  0.866025,  -1.037500},
      {0.000000, -0.866025,  -1.037500}
  };
  
  CH2O_positions *= atoms::angstrom_to_bohr;
  
  arma::vec bonds = connectivity::get_bonds(CH2O_symbols, CH2O_positions) /
                    atoms::angstrom_to_bohr;
  
  cout << "Bonds:\n" << endl;
  cout << bonds << endl;
  
  assert( arma::norm( bonds - arma::vec{1.2000, 1.0000, 1.0000}) < 1e-6);
  
  return 0;
}