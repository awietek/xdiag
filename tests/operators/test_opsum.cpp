#include "../catch.hpp"

#include <iostream>
#include <xdiag/operators/opsum.hpp>

TEST_CASE("opsum", "[operators]") {
  using namespace xdiag;

  OpSum ops;
  ops += "J1" * Op("SdotS", {0, 1});
  ops += "J1" * Op("SdotS", {1, 2});
  ops += "J1" * Op("SdotS", {2, 3});
  ops += "J1" * Op("SdotS", {3, 4});
  ops += "J1" * Op("SdotS", {4, 5});
  ops += "J1" * Op("SdotS", {5, 0});
  ops += "J2" * Op("SdotS", {0, 2});
  ops += "J2" * Op("SdotS", {1, 3});
  ops += "J2" * Op("SdotS", {2, 4});
  ops += "J2" * Op("SdotS", {3, 5});
  ops += "J2" * Op("SdotS", {4, 0});
  ops += "J2" * Op("SdotS", {5, 1});

  double a = 1.0;
  complex b(2.0, 3.0);

  ops["J1"] = a;
  ops["J2"] = b;

  OpSum ops2;
  ops2 += a * Op("SdotS", {0, 1});
  ops2 += a * Op("SdotS", {1, 2});
  ops2 += a * Op("SdotS", {2, 3});
  ops2 += a * Op("SdotS", {3, 4});
  ops2 += a * Op("SdotS", {4, 5});
  ops2 += a * Op("SdotS", {5, 0});
  ops2 += b * Op("SdotS", {0, 2});
  ops2 += b * Op("SdotS", {1, 3});
  ops2 += b * Op("SdotS", {2, 4});
  ops2 += b * Op("SdotS", {3, 5});
  ops2 += b * Op("SdotS", {4, 0});
  ops2 += b * Op("SdotS", {5, 1});

  REQUIRE(ops.plain() == ops2);
}
