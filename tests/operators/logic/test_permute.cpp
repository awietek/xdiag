#include "../../catch.hpp"

#include <xdiag/common.hpp>
#include <xdiag/operators/logic/permute.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

TEST_CASE("permute", "[operators]") try {
  Log("Testing permute of operator");

  Op op("A", {0, 1});
  Permutation p({1, 0});
  REQUIRE(permute(op, p) == Op("A", {1, 0}));

  op = Op("B", {0, 1, 2, 3, 4, 5});
  p = Permutation({5, 4, 3, 2, 1, 0});
  REQUIRE(permute(op, p) == Op("B", {5, 4, 3, 2, 1, 0}));

  op = Op("C", {0, 2, 4, 5});
  REQUIRE(permute(op, p) == Op("C", {5, 3, 1, 0}));

  OpSum ops1, ops2;
  for (int i=0; i<6; ++i){
    ops1 += Op("D", {i, (i+1)%6});
    ops2 += Op("D", {(5-i)%6, (5-i-1 + 6)%6});
  }
  auto opsp = permute(ops1, p);
  REQUIRE(opsp==ops2);
  // XDIAG_SHOW(ops1);
  // XDIAG_SHOW(opsp);
  // XDIAG_SHOW(ops2);
  
  // Error messages
  // op = Op("E", -1);
  // permute(op, p);

  // op = Op("F", 6);
  // permute(op, p);
  
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
