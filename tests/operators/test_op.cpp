#include "../catch.hpp"

#include <iostream>
#include <xdiag/operators/op.hpp>

TEST_CASE("Op", "[operators]") {
  using namespace xdiag;

  auto o1 = Op("HB", 1.0, {0, 1});
  auto o2 = Op("HB", {0, 1}) * 1.0;
  auto o3 = 1.0 * Op("HB", {0, 1});

  REQUIRE(o1 == o2);
  REQUIRE(o1 == o3);
  REQUIRE(o2 == o3);

  o1 = Op("HB", complex(1.0, 2.0), {0, 1});
  o2 = Op("HB", {0, 1}) * complex(1.0, 2.0);
  o3 = complex(1.0, 2.0) * Op("HB", {0, 1});

  REQUIRE(o1 == o2);
  REQUIRE(o1 == o3);
  REQUIRE(o2 == o3);

  o1 = Op("HB", "J", {0, 1});
  o2 = Op("HB", {0, 1}) * "J";
  o3 = "J" * Op("HB", {0, 1});

  REQUIRE(o1 == o2);
  REQUIRE(o1 == o3);
  REQUIRE(o2 == o3);

  arma::mat A(2, 2, arma::fill::randu);
  o1 = Op("SX", A, 1);
  o2 = 2.0 * o1;
  o3 = o1 * 2.0;

  // XDIAG_SHOW(o1);
  // XDIAG_SHOW(o2);
  // XDIAG_SHOW(o3);

  REQUIRE(o2 == o3);

  arma::cx_mat B(2, 2, arma::fill::randu);
  o1 = Op("SX", B, 1);
  o2 = 2.0 * o1;
  o3 = o1 * 2.0;

  // XDIAG_SHOW(o1);
  // XDIAG_SHOW(o2);
  // XDIAG_SHOW(o3);

  REQUIRE(o2 == o3);
}
