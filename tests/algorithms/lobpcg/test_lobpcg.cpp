// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <iostream>

#include "../../blocks/electron/testcases_electron.hpp"

#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algorithms/lobpcg/lobpcg.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/utils/xdiag_show.hpp>

using namespace xdiag;

template <typename idx_t, typename coeff_t> void test_lobpcg_csr_matrix() {
  using namespace xdiag::testcases::electron;
  OpSum ops;
  ops["U"] = 5.0;

  int nsites = 6;
  int num_eigenvalue = 1;

  ops = freefermion_alltoall(nsites);
  int nup = 3;
  int ndn = 3;
  auto block = Electron(nsites, nup, ndn);

  auto mult = [&](arma::mat const &X, arma::mat &Y) {
    apply(ops, block, X, block, Y);
  };

  int64_t nrows = size(block);
  int64_t ncols = 5;
  auto X0 = arma::mat(nrows, ncols, arma::fill::randn);

  XDIAG_SHOW(eigval0(ops, block));
  Log.set_verbosity(1);

  lobpcg::lobpcg(mult, X0, 1e-6, 30, true, 20);
}

TEST_CASE("lobpcg", "[lobpcg]") try {
  
  Log("lobpcg_csr_matrix int32_t double");
  test_lobpcg_csr_matrix<int32_t, double>();
  // Log("lobpcg_csr_matrix int32_t complex");
  // test_lobpcg_csr_matrix<int32_t, complex>();
  // Log("lobpcg_csr_matrix int64_t double");
  // test_lobpcg_csr_matrix<int64_t, double>();
  // Log("lobpcg_csr_matrix int64_t complex");
  // test_lobpcg_csr_matrix<int64_t, complex>();
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
