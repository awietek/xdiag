// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>
#include <xdiag/math/lapack.hpp>
#include <xdiag/utils/error.hpp>

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

    arma::fvec fd0 = arma::conv_to<arma::fvec>::from(diag);
    arma::fvec fe0 = arma::conv_to<arma::fvec>::from(offdiag);

    arma::vec d = diag;
    arma::vec e = offdiag;
    math::lapack::sterf(n, d.memptr(), e.memptr());
    REQUIRE(arma::norm(d - expected) < 1e-12);

    arma::fvec fd = arma::conv_to<arma::fvec>::from(diag);
    arma::fvec fe = arma::conv_to<arma::fvec>::from(offdiag);
    math::lapack::sterf(n, fd.memptr(), fe.memptr());
    REQUIRE(arma::norm(arma::conv_to<arma::vec>::from(fd) - expected) < 1e-4);

    // stebz computes a range of eigenvalues and leaves its input alone
    int64_t k = std::min<int64_t>(n, 3);
    arma::vec lowest(k);
    math::lapack::stebz(n, diag.memptr(), offdiag.memptr(), 1, k,
                        lowest.memptr());
    REQUIRE(arma::norm(lowest - expected.head(k)) < 1e-12);

    arma::vec highest(1);
    math::lapack::stebz(n, diag.memptr(), offdiag.memptr(), n, n,
                        highest.memptr());
    REQUIRE(std::abs(highest(0) - expected(n - 1)) < 1e-12);

    // stedc, checked by residual and orthonormality so signs do not matter
    arma::vec sd = diag;
    arma::vec se = offdiag;
    arma::mat evecs(n, n);
    math::lapack::stedc(n, sd.memptr(), se.memptr(), evecs.memptr());
    REQUIRE(arma::norm(sd - expected) < 1e-12);
    REQUIRE(arma::norm(tmat * evecs - evecs * arma::diagmat(sd)) < 1e-10);
    REQUIRE(arma::norm(evecs.t() * evecs - arma::eye(n, n)) < 1e-10);

    arma::fvec flowest(k);
    math::lapack::stebz(n, fd0.memptr(), fe0.memptr(), 1, k, flowest.memptr());
    REQUIRE(arma::norm(arma::conv_to<arma::vec>::from(flowest) -
                       expected.head(k)) < 1e-4);
  }
}

TEST_CASE("lapack_stebz_range", "[math]") {
  using namespace xdiag;

  arma::vec diag(5, arma::fill::randn);
  arma::vec offdiag(4, arma::fill::randn);
  arma::vec out(5);

  REQUIRE_THROWS_AS(
      math::lapack::stebz(5, diag.memptr(), offdiag.memptr(), 0, 2,
                          out.memptr()),
      Error);
  REQUIRE_THROWS_AS(
      math::lapack::stebz(5, diag.memptr(), offdiag.memptr(), 3, 2,
                          out.memptr()),
      Error);
  REQUIRE_THROWS_AS(
      math::lapack::stebz(5, diag.memptr(), offdiag.memptr(), 1, 6,
                          out.memptr()),
      Error);
}
