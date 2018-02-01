#include "../../include/catch/catch.hpp"

#include "libirc/periodic_table.h"

#include <iomanip>
#include <iostream>
#include <stdexcept>

TEST_CASE("Test periodic table", "[periodic_table]") {
  using namespace std;

  using namespace irc;
  using namespace periodic_table;

  bool verbose{false};

  // Check periodic table size for supported elements
  SECTION("size") { REQUIRE(pt_size == 96); }

  // Check valid atomic numbers
  SECTION("atomic numbers") {
    bool valid{false};
    for (size_t i{0}; i < pt_size + 1; i++) {
      valid = valid_atomic_number(i);

      if (valid) {
        REQUIRE(i > 0);
        REQUIRE(i < pt_size);
      } else {
        REQUIRE((i == 0 or i == pt_size));
      }
    }
  }

  // Check conversion from symbol to atomic number
  SECTION("atomic number from symbol") {
    for (size_t i{1}; i < pt_size; i++) {
      REQUIRE(atomic_number(pt_symbols[i]) == i);
    }
  }

  SECTION("invalid symbol") {
    bool exception{false};
    try {
      atomic_number("ABC");
    } catch (const std::logic_error &e) {
      exception = true;
    }

    REQUIRE(exception == true);
  }
}