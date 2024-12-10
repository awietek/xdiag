#include "../catch.hpp"

#include <iostream>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

TEST_CASE("Op", "[operators]") {
  using namespace xdiag;

  auto o1 = Op("HB", {0, 1}) * 1.0;
  auto o2 = 1.0 * Op("HB", {0, 1});

  REQUIRE(o1 == o2);

  o1 = Op("HB", {0, 1}) * complex(1.0, 2.0);
  o2 = complex(1.0, 2.0) * Op("HB", {0, 1});

  REQUIRE(o1 == o2);

  o1 = Op("HB", {0, 1}) * "J";
  o2 = "J" * Op("HB", {0, 1});

  REQUIRE(o1 == o2);

  arma::mat A(2, 2, arma::fill::randu);
  auto op = Op("SX", 1, A);
  o1 = 2.0 * op;
  o2 = op * 2.0;

  XDIAG_SHOW(o1);
  XDIAG_SHOW(o2);

  REQUIRE(o1 == o2);

  arma::cx_mat B(2, 2, arma::fill::randu);
  op = Op("SX", 1, B);
  o1 = 2.0 * op;
  o2 = op * 2.0;

  XDIAG_SHOW(o1);
  XDIAG_SHOW(o2);

  REQUIRE(o1 == o2);
}
