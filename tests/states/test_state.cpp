// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/xdiag_show.hpp>

TEST_CASE("state", "[states]") try {
  using namespace xdiag;
  using namespace arma;

  auto block = Spinhalf(4, 2);
  auto psir = State(block);
  REQUIRE(isreal(psir));
  // XDIAG_SHOW(psir);
  // XDIAG_SHOW(psir.matrix());
  REQUIRE(norm(psir.matrix()) < 1e-12);

  auto psic = State(block, false);
  REQUIRE(!isreal(psic));
  // XDIAG_SHOW(psic);
  // XDIAG_SHOW(psic.matrixC());
  REQUIRE(norm(psic.matrixC()) < 1e-12);

  arma::vec vr = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  auto psivr = State(block, vr);
  REQUIRE(isreal(psivr));
  // XDIAG_SHOW(psivr.vector());
  REQUIRE(norm(psivr.vector() - vr) < 1e-12);

  arma::vec vi = {6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
  auto vc = arma::cx_vec(vr, vi);
  auto psivc = State(block, vc);
  // XDIAG_SHOW(psivc.vectorC());
  REQUIRE(!isreal(psivc));
  REQUIRE(norm(psivc.vectorC() - vc) < 1e-12);
  REQUIRE(norm(psivc.real().vector() - vr) < 1e-12);
  REQUIRE(norm(psivc.imag().vector() - vi) < 1e-12);
    
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("state matrix constructors", "[states]") try {
  using namespace xdiag;

  auto block = Spinhalf(4, 2); // dim 6

  // --- real multi-column matrix constructor + accessors ---
  {
    arma::mat M(6, 3, arma::fill::randu);
    auto psi = State(block, M);
    REQUIRE(isreal(psi));
    REQUIRE(psi.ncols() == 3);
    REQUIRE(psi.nrows() == 6);
    REQUIRE(arma::norm(psi.matrix() - M) < 1e-12);

    // col(n) extracts a single-column State
    auto c1 = psi.col(1);
    REQUIRE(c1.ncols() == 1);
    REQUIRE(arma::norm(c1.vector() - M.col(1)) < 1e-12);
  }

  // --- complex multi-column matrix constructor ---
  {
    arma::cx_mat M(6, 2, arma::fill::randu);
    auto psi = State(block, M);
    REQUIRE(!isreal(psi));
    REQUIRE(psi.ncols() == 2);
    REQUIRE(arma::norm(psi.matrixC() - M) < 1e-12);
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("state constructor and accessor errors", "[states]") try {
  using namespace xdiag;

  auto block = Spinhalf(4, 2); // dim 6

  // --- constructor size mismatches throw ---
  REQUIRE_THROWS_AS(State(block, arma::vec(5, arma::fill::zeros)), Error);
  REQUIRE_THROWS_AS(State(block, arma::cx_vec(5, arma::fill::zeros)), Error);
  REQUIRE_THROWS_AS(State(block, arma::mat(5, 2, arma::fill::zeros)), Error);
  REQUIRE_THROWS_AS(State(block, arma::cx_mat(5, 2, arma::fill::zeros)), Error);

  // --- real state: complex accessors throw ---
  {
    auto psi = State(block, true);
    REQUIRE_THROWS_AS(psi.vectorC(), Error);
    REQUIRE_THROWS_AS(psi.matrixC(), Error);
    // out-of-range / negative column indices throw
    REQUIRE_THROWS_AS(psi.vector(5), Error);  // only 1 column
    REQUIRE_THROWS_AS(psi.vector(-1), Error);
  }

  // --- complex state: real accessors throw ---
  {
    auto psi = State(block, false);
    REQUIRE_THROWS_AS(psi.vector(), Error);
    REQUIRE_THROWS_AS(psi.matrix(), Error);
    REQUIRE_THROWS_AS(psi.vectorC(-1), Error);
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
