// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <cmath>

#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/states/norm.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/error.hpp>

TEST_CASE("norm", "[states]") try {
  using namespace xdiag;

  auto block = Spinhalf(4, 2); // dimension C(4,2) = 6

  // --- Real state: norms of a known vector ---
  {
    arma::vec v = {0.0, 3.0, 0.0, 4.0, 0.0, 0.0};
    auto psi = State(block, v);
    REQUIRE(isreal(psi));
    REQUIRE(std::abs(norm(psi) - 5.0) < 1e-12);       // sqrt(9 + 16)
    REQUIRE(std::abs(norm1(psi) - 7.0) < 1e-12);      // |3| + |4|
    REQUIRE(std::abs(norminf(psi) - 4.0) < 1e-12);    // max |.|
  }

  // --- Complex state: same magnitudes distributed over re/im ---
  {
    arma::vec re = {3.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    arma::vec im = {0.0, 4.0, 0.0, 0.0, 0.0, 0.0};
    auto v = arma::cx_vec(re, im);
    auto psi = State(block, v);
    REQUIRE(!isreal(psi));
    REQUIRE(std::abs(norm(psi) - 5.0) < 1e-12);       // sqrt(3^2 + 4^2)
    REQUIRE(std::abs(norm1(psi) - 7.0) < 1e-12);      // |3| + |4i|
    REQUIRE(std::abs(norminf(psi) - 4.0) < 1e-12);
  }

  // --- Zero / invalid state returns 0 ---
  {
    State empty;
    REQUIRE(!isvalid(empty));
    REQUIRE(norm(empty) == 0.);
    REQUIRE(norm1(empty) == 0.);
    REQUIRE(norminf(empty) == 0.);
  }

  // --- Multi-column state: norms are ill-defined and must throw ---
  {
    auto psi = State(block, true, 2); // 2 columns
    REQUIRE(psi.ncols() == 2);
    REQUIRE_THROWS_AS(norm(psi), Error);
    REQUIRE_THROWS_AS(norm1(psi), Error);
    REQUIRE_THROWS_AS(norminf(psi), Error);
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
