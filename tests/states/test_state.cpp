// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"
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
  XDIAG_SHOW(psivr.vector());
  REQUIRE(norm(psivr.vector() - vr) < 1e-12);

  arma::vec vi = {6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
  auto vc = arma::cx_vec(vr, vi);
  auto psivc = State(block, vc);
  XDIAG_SHOW(psivc.vectorC());
  REQUIRE(!isreal(psivc));
  REQUIRE(norm(psivc.vectorC() - vc) < 1e-12);
  REQUIRE(norm(psivc.real().vector() - vr) < 1e-12);
  REQUIRE(norm(psivc.imag().vector() - vi) < 1e-12);
    
} catch (xdiag::Error const &e) {
  error_trace(e);
}
