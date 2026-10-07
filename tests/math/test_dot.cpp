// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/math/dot.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

TEST_CASE("dot_chunked", "[math]") {
  Log("dot_chunked test ...");

  // spans several chunks, so the in-order reduction of the chunk sums runs
  for (int64_t n : {1, 17, 4096, 10000}) {
    arma::vec v(n, arma::fill::randn);
    arma::vec w(n, arma::fill::randn);
    REQUIRE(std::abs(math::dot_chunked(v, w) - arma::dot(v, w)) < 1e-10);

    arma::cx_vec vc(n, arma::fill::randn);
    arma::cx_vec wc(n, arma::fill::randn);
    REQUIRE(std::abs(math::dot_chunked(vc, wc) - arma::cdot(vc, wc)) < 1e-10);
  }

  arma::vec v(7, arma::fill::randn);
  arma::vec w(5, arma::fill::randn);
  REQUIRE_THROWS_AS(math::dot_chunked(v, w), Error);

  arma::cx_vec vc(7, arma::fill::randn);
  arma::cx_vec wc(5, arma::fill::randn);
  REQUIRE_THROWS_AS(math::dot_chunked(vc, wc), Error);
}
