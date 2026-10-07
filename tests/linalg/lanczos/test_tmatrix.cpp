// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/linalg/lanczos/tmatrix.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

TEST_CASE("tmatrix_eigenvalues_lowest", "[lanczos]") {
  Log("tmatrix_eigenvalues_lowest test ...");

  std::vector<double> alphas{1.5, -0.5, 0.25, 2.0, -1.25};
  std::vector<double> betas{0.5, -0.75, 1.25, 0.5, 0.0};
  Tmatrix tmat(alphas, betas);
  arma::vec full = tmat.eigenvalues();
  int64_t n = tmat.size();

  // bisection path
  for (int64_t k = 1; k < n; ++k) {
    REQUIRE(arma::norm(tmat.eigenvalues_lowest(k) - full.head(k)) < 1e-12);
  }

  // a request covering the whole matrix falls back to the full spectrum
  REQUIRE(arma::norm(tmat.eigenvalues_lowest(n) - full) < 1e-12);
  REQUIRE(arma::norm(tmat.eigenvalues_lowest(n + 3) - full) < 1e-12);
  REQUIRE(Tmatrix().eigenvalues_lowest(1).n_elem == 0);
}
