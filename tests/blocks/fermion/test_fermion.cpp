// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/matrices/matrix.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/xdiag_show.hpp>

using namespace xdiag;

TEST_CASE("fermion", "[fermion]") try {
  Log("Fermion anti-commutation relations");
  int nsites = 4;
  auto b = Fermion(nsites);
  int dim = b.dim();

  for (int i = 0; i < nsites; ++i) {
    for (int j = 0; j < nsites; ++j) {
      arma::mat cdagi = matrix(Op("Cdag", i), b);
      arma::mat cdagj = matrix(Op("Cdag", j), b);
      arma::mat ci = matrix(Op("C", i), b);
      arma::mat cj = matrix(Op("C", j), b);
      arma::mat id = arma::eye(dim, dim);
      arma::mat zeros = arma::mat(dim, dim, arma::fill::zeros);
      arma::mat cdagcdag = cdagi * cdagj + cdagj * cdagi;
      arma::mat cdagc = cdagi * cj + cj * cdagi;
      arma::mat ccdag = ci * cdagj + cdagj * ci;
      arma::mat cc = ci * cj + cj * ci;

      if (i == j) {
        REQUIRE(isapprox(cdagc, id));
        REQUIRE(isapprox(ccdag, id));
        REQUIRE(isapprox(cdagcdag, zeros));
        REQUIRE(isapprox(cc, zeros));
      } else {
        REQUIRE(isapprox(cdagc, zeros));
        REQUIRE(isapprox(ccdag, zeros));
        REQUIRE(isapprox(cdagcdag, zeros));
        REQUIRE(isapprox(cc, zeros));
      }
    }
  }
} catch (xdiag::Error e) {
  error_trace(e);
}
