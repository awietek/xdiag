#include "../catch.hpp"

#include <iostream>
#include <xdiag/operators/opsum.hpp>

TEST_CASE("opsum", "[operators]") {
  using namespace xdiag;

  OpSum ops;
  ops += "J1" * Op("SDOTS", {0, 1});
  ops += "J1" * Op("SDOTS", {1, 2});
  ops += "J1" * Op("SDOTS", {2, 3});
  ops += "J1" * Op("SDOTS", {3, 4});
  ops += "J1" * Op("SDOTS", {4, 5});
  ops += "J1" * Op("SDOTS", {5, 0});
  ops += "J2" * Op("SDOTS", {0, 2});
  ops += "J2" * Op("SDOTS", {1, 3});
  ops += "J2" * Op("SDOTS", {2, 4});
  ops += "J2" * Op("SDOTS", {3, 5});
  ops += "J2" * Op("SDOTS", {4, 0});
  ops += "J2" * Op("SDOTS", {5, 1});

  double a = 1.0;
  complex b(2.0, 3.0);

  ops["J1"] = a;
  ops["J2"] = b;

  OpSum ops2;
  ops2 += a * Op("SDOTS", {0, 1});
  ops2 += a * Op("SDOTS", {1, 2});
  ops2 += a * Op("SDOTS", {2, 3});
  ops2 += a * Op("SDOTS", {3, 4});
  ops2 += a * Op("SDOTS", {4, 5});
  ops2 += a * Op("SDOTS", {5, 0});
  ops2 += b * Op("SDOTS", {0, 2});
  ops2 += b * Op("SDOTS", {1, 3});
  ops2 += b * Op("SDOTS", {2, 4});
  ops2 += b * Op("SDOTS", {3, 5});
  ops2 += b * Op("SDOTS", {4, 0});
  ops2 += b * Op("SDOTS", {5, 1});

  REQUIRE(ops.plain() == ops2);
}
