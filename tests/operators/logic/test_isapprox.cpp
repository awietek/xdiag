#include "../../catch.hpp"

#include <xdiag/operators/logic/isapprox.hpp>
using namespace xdiag;

TEST_CASE("isapprox", "[operators]") try {
  using namespace arma;

  Log("Testing isapprox of operators");
  auto o1 = Op("SdotS", {0, 1});
  auto o2 = Op("SdotS", {1, 0});
  REQUIRE(isapprox(o1, o2));

  o1 = Op("ScalarChirality", {0, 1, 2});
  o2 = Op("ScalarChirality", {1, 2, 0});
  REQUIRE(isapprox(o1, o2));

  o1 = Op("ScalarChirality", {0, 1, 2});
  o2 = Op("ScalarChirality", {2, 1, 0});
  REQUIRE(!isapprox(o1, o2));

  auto os1 = OpSum();
  os1 += 1.0 * o1;
  auto os2 = OpSum();
  os2 += -1.0 * o2;
  REQUIRE(isapprox(os1, os2));

  os1 = OpSum();
  os1 += 1.0 * o1;
  os2 = OpSum();
  os2 += 1.0 * o2;
  auto factor = isapprox_multiple(os1, os2);
  REQUIRE(factor);
  REQUIRE(*factor == -1.0);

  std::vector<std::string> types = {"Exchange", "Hop", "Hopup", "Hopdn"};
  for (auto type : types) {
    os1 = OpSum();
    os1 += complex(0, 1) * Op("Hop", {0, 1});
    os2 = OpSum();
    os2 += complex(0, 1) * Op("Hop", {1, 0});
    REQUIRE(!isapprox(os1, os2));
    
    os1 = OpSum();
    os1 += complex(0, 1) * Op(type, {0, 1});
    os2 = OpSum();
    os2 += complex(0, -1) * Op(type, {1, 0});
    REQUIRE(isapprox(os1, os2));
  }

  Log("done");
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
