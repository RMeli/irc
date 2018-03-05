#include "catch.hpp"

#include "libirc/periodic_table.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace irc {
namespace pt = periodic_table;
}

TEST_CASE("Test periodic table", "[periodic_table]") {

  // Check periodic table size for supported elements
  REQUIRE(irc::pt::pt_size == 96);
  REQUIRE(irc::pt::symbols.size() == irc::pt::pt_size);
  REQUIRE(irc::pt::covalent_radii.size() == irc::pt::pt_size);
  REQUIRE(irc::pt::vdw_radii.size() == irc::pt::pt_size);
  REQUIRE(irc::pt::masses.size() == irc::pt::pt_size);

  // Check valid atomic numbers
  for (size_t atomic_num = 1; atomic_num < irc::pt::pt_size; atomic_num++) {
    CAPTURE(atomic_num);
    CHECK(irc::pt::valid_atomic_number(atomic_num));
  }

  // Check conversion from symbol to atomic number
  for (size_t atomic_num = 1; atomic_num < irc::pt::pt_size; atomic_num++) {
    CAPTURE(atomic_num);
    CAPTURE(irc::pt::symbols[atomic_num]);
    REQUIRE(irc::pt::atomic_number(irc::pt::symbols[atomic_num]) == atomic_num);
  }

  CHECK_THROWS_AS(irc::pt::atomic_number("ABC"), std::logic_error);
  CHECK_THROWS_AS(irc::pt::atomic_number(" H"), std::logic_error);
  CHECK_THROWS_AS(irc::pt::atomic_number("H "), std::logic_error);
  CHECK_THROWS_AS(irc::pt::atomic_number(" H "), std::logic_error);
  CHECK_THROWS_AS(irc::pt::atomic_number("H,"), std::logic_error);
  CHECK_THROWS_AS(irc::pt::atomic_number("H\n"), std::logic_error);
}
