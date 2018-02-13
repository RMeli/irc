#include "catch.hpp"

#include "libirc/periodic_table.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>

TEST_CASE("Test periodic table", "[periodic_table]") {
  using namespace std;

  using namespace irc;
  using namespace periodic_table;

  // Check periodic table size for supported elements
  REQUIRE(pt_size == 96);
  REQUIRE(pt_symbols.size() == pt_size);
  REQUIRE(pt_covalent_radii.size() == pt_size);
  REQUIRE(pt_vdv_radii.size() == pt_size);
  REQUIRE(pt_masses.size() == pt_size);

  // Check valid atomic numbers
  for (size_t atomic_num = 1; atomic_num < pt_size; atomic_num++) {
    CAPTURE(atomic_num);
    CHECK(valid_atomic_number(atomic_num));
  }

  // Check conversion from symbol to atomic number
  for (size_t atomic_num = 1; atomic_num < pt_size; atomic_num++) {
    CAPTURE(atomic_num);
    CAPTURE(pt_symbols[atomic_num]);
    REQUIRE(atomic_number(pt_symbols[atomic_num]) == atomic_num);
  }

  CHECK_THROWS_AS(atomic_number("ABC"), std::logic_error);
  CHECK_THROWS_AS(atomic_number(" H"), std::logic_error);
  CHECK_THROWS_AS(atomic_number("H "), std::logic_error);
  CHECK_THROWS_AS(atomic_number(" H "), std::logic_error);
  CHECK_THROWS_AS(atomic_number("H,"), std::logic_error);
  CHECK_THROWS_AS(atomic_number("H\n"), std::logic_error);
}
