// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>
#include <xdiag/math/lapack.hpp>

TEST_CASE("lapack_sterf", "[math]") {
  using namespace xdiag;

  for (int64_t n : {1, 2, 5, 37}) {
    arma::vec diag(n, arma::fill::randn);
    arma::vec offdiag(std::max<int64_t>(n - 1, 1), arma::fill::randn);

    arma::mat tmat(n, n, arma::fill::zeros);
    tmat.diag() = diag;
    if (n > 1) {
      tmat.diag(1) = offdiag.head(n - 1);
      tmat.diag(-1) = offdiag.head(n - 1);
    }
    arma::vec expected = arma::eig_sym(tmat);

    arma::vec d = diag;
    arma::vec e = offdiag;
    math::lapack::sterf(n, d.memptr(), e.memptr());
    REQUIRE(arma::norm(d - expected) < 1e-12);

    arma::fvec fd = arma::conv_to<arma::fvec>::from(diag);
    arma::fvec fe = arma::conv_to<arma::fvec>::from(offdiag);
    math::lapack::sterf(n, fd.memptr(), fe.memptr());
    REQUIRE(arma::norm(arma::conv_to<arma::vec>::from(fd) - expected) < 1e-4);
  }
}
