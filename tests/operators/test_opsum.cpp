#include "../catch.hpp"

#include <iostream>
#include <xdiag/operators/opsum.hpp>

TEST_CASE("opsum", "[operators]") {
  using namespace xdiag;

  OpSum ops;
  ops += "J1" * Op("HB", {0, 1});
  ops += "J1" * Op("HB", {1, 2});
  ops += "J1" * Op("HB", {2, 3});
  ops += "J1" * Op("HB", {3, 4});
  ops += "J1" * Op("HB", {4, 5});
  ops += "J1" * Op("HB", {5, 0});
  ops += "J2" * Op("HB", {0, 2});
  ops += "J2" * Op("HB", {1, 3});
  ops += "J2" * Op("HB", {2, 4});
  ops += "J2" * Op("HB", {3, 5});
  ops += "J2" * Op("HB", {4, 0});
  ops += "J2" * Op("HB", {5, 1});

  double a = 1.0;
  complex b(2.0, 3.0);

  ops["J1"] = a;
  ops["J2"] = b;

  OpSum ops2;
  ops2 += a * Op("HB", {0, 1});
  ops2 += a * Op("HB", {1, 2});
  ops2 += a * Op("HB", {2, 3});
  ops2 += a * Op("HB", {3, 4});
  ops2 += a * Op("HB", {4, 5});
  ops2 += a * Op("HB", {5, 0});
  ops2 += b * Op("HB", {0, 2});
  ops2 += b * Op("HB", {1, 3});
  ops2 += b * Op("HB", {2, 4});
  ops2 += b * Op("HB", {3, 5});
  ops2 += b * Op("HB", {4, 0});
  ops2 += b * Op("HB", {5, 1});

  REQUIRE(ops == ops2);
}
