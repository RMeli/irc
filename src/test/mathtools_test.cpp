#include "catch.hpp"

#include "libirc/mathtools.h"

using namespace irc;

TEST_CASE("Angle in the interval ]-pi,pi]") {

  using tools::constants::pi;
  using namespace tools::math;

  SECTION("Angle in ]-pi,pi]") {
    double angle{pi / 2.};

    Approx target{pi / 2.};

    target.margin(1e-6);

    REQUIRE(pirange_rad(angle) == target);
  }

  SECTION("Angle in ]pi,2*pi]") {
    double angle{pi / 2. + pi};

    Approx target{-pi / 2.};

    target.margin(1e-6);

    REQUIRE(pirange_rad(angle) == target);
  }

  SECTION("Angle in ]-2*pi,-pi]") {
    double angle{-pi / 2. - pi};

    Approx target{pi / 2.};

    target.margin(1e-6);

    REQUIRE(pirange_rad(angle) == target);
  }

  SECTION("Angle") {
    double angle{pi / 2. + 6 * pi};

    Approx target{pi / 2.};

    target.margin(1e-6);

    REQUIRE(pirange_rad(angle) == target);
  }

  SECTION("Angle") {
    double angle{pi / 2. - 5 * pi};

    Approx target{-pi / 2.};

    target.margin(1e-6);

    REQUIRE(pirange_rad(angle) == target);
  }
}

TEST_CASE("Angle in the interval ]-180,180]") {

  using tools::constants::pi;
  using namespace tools::math;

  SECTION("Angle in ]-180,180]") {
    double angle{90.};

    Approx target{90.};

    target.margin(1e-6);

    REQUIRE(pirange_deg(angle) == target);
  }

  SECTION("Angle in ]180,360]") {
    double angle{270.};

    Approx target{-90};

    target.margin(1e-6);

    REQUIRE(pirange_deg(angle) == target);
  }

  SECTION("Angle in ]-360,-180]") {
    double angle{-270.};

    Approx target{90.};

    target.margin(1e-6);

    REQUIRE(pirange_deg(angle) == target);
  }

  SECTION("Angle in ]360,540]") {
    double angle{450.};

    Approx target{90.};

    target.margin(1e-6);

    REQUIRE(pirange_deg(angle) == target);
  }

  SECTION("Angle in [-450,360[") {
    double angle{-450};

    Approx target{-90.};

    target.margin(1e-6);

    REQUIRE(pirange_deg(angle) == target);
  }
}
