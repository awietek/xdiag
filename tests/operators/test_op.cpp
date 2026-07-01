// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <iostream>
#include <vector>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>

TEST_CASE("Op construction and properties", "[operators]") try {
  using namespace xdiag;

  // --- Default construction ---
  {
    Op op;
    // No assertions needed; must not throw
  }

  // --- Type-only construction ---
  {
    Op op("HubbardU");
    REQUIRE(op.type() == "HubbardU");
    REQUIRE(op.hassites() == false);
    REQUIRE(op.hasmatrix() == false);
  }

  // --- Single-site construction ---
  {
    Op op("Sz", int64_t(3));
    REQUIRE(op.type() == "Sz");
    REQUIRE(op.size() == 1);
    REQUIRE(op[0] == 3);
    REQUIRE(op.hassites() == true);
    REQUIRE(op.hasmatrix() == false);
  }

  // --- Multi-site construction ---
  {
    std::vector<int64_t> sites = {0, 1, 2};
    Op op("SdotS", sites);
    REQUIRE(op.type() == "SdotS");
    REQUIRE(op.size() == 3);
    REQUIRE(op[0] == 0);
    REQUIRE(op[1] == 1);
    REQUIRE(op[2] == 2);
    REQUIRE(op.sites() == sites);
    REQUIRE(op.hassites() == true);
  }

  // --- Matrix construction (real arma::mat) ---
  {
    arma::mat M(2, 2, arma::fill::randu);
    Op op("Matrix", int64_t(0), M);
    REQUIRE(op.type() == "Matrix");
    REQUIRE(op.hassites() == true);
    REQUIRE(op.hasmatrix() == true);
  }

  // --- Matrix construction (complex arma::cx_mat) ---
  {
    arma::cx_mat M(2, 2, arma::fill::randu);
    Op op("Matrix", int64_t(1), M);
    REQUIRE(op.type() == "Matrix");
    REQUIRE(op.hassites() == true);
    REQUIRE(op.hasmatrix() == true);
  }

} catch (xdiag::Error e) {
  error_trace(e);
}

TEST_CASE("Op equality", "[operators]") try {
  using namespace xdiag;

  // Same type and site: equal
  REQUIRE(Op("Sz", int64_t(0)) == Op("Sz", int64_t(0)));

  // Different type: not equal
  REQUIRE(Op("Sz", int64_t(0)) != Op("Sp", int64_t(0)));

  // Different site: not equal
  REQUIRE(Op("Sz", int64_t(0)) != Op("Sz", int64_t(1)));

  // Same type+site, different matrix: not equal
  {
    arma::mat M1(2, 2, arma::fill::ones);
    arma::mat M2(2, 2, arma::fill::zeros);
    Op op1("Matrix", int64_t(0), M1);
    Op op2("Matrix", int64_t(0), M2);
    REQUIRE(op1 != op2);
  }

} catch (xdiag::Error e) {
  error_trace(e);
}

TEST_CASE("Op", "[operators]") try {
  using namespace xdiag;

  // Original test: *, "J"*, complex* constructors for OpSum

  auto o1 = Op("SdotS", {0, 1}) * 1.0;
  auto o2 = 1.0 * Op("SdotS", {0, 1});

  REQUIRE(o1 == o2);

  o1 = Op("SdotS", {0, 1}) * complex(1.0, 2.0);
  o2 = complex(1.0, 2.0) * Op("SdotS", {0, 1});

  REQUIRE(o1 == o2);

  o1 = Op("SdotS", {0, 1}) * "J";
  o2 = "J" * Op("SdotS", {0, 1});

  REQUIRE(o1 == o2);

  arma::mat A(2, 2, arma::fill::randu);
  auto op = Op("SX", 1, A);
  o1 = 2.0 * op;
  o2 = op * 2.0;

  REQUIRE(o1 == o2);

  arma::cx_mat B(2, 2, arma::fill::randu);
  op = Op("SX", 1, B);
  o1 = 2.0 * op;
  o2 = op * 2.0;

  REQUIRE(o1 == o2);

} catch (xdiag::Error e) {
  error_trace(e);
}
