// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <iostream>

#include "../blocks/spinhalf/testcases_spinhalf.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/time_evolution/exp_sym_v.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/common.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/logger.hpp>

TEST_CASE("exp_sym_v", "[algorithms]") try {
  using namespace xdiag;
  for (int64_t N = 2; N <= 10; ++N) {
    Log("Testing exp_sym_v on all-to-all HB model, N={}", N);

    auto block = Spinhalf(N);
    auto psi0 = random_state(block);
    auto psi0c = psi0;
    make_complex(psi0c);
    auto ops = xdiag::testcases::spinhalf::HB_alltoall(N);

    auto mult = [&](arma::vec const &v, arma::vec &w) {
      apply(ops, block, v, block, w);
    };
    auto dot_f = [&block](arma::vec const &v, arma::vec const &w) {
      return dot(block, v, w);
    };
    auto multC = [&](arma::cx_vec const &v, arma::cx_vec &w) {
      apply(ops, block, v, block, w);
    };
    auto dot_fC = [&block](arma::cx_vec const &v, arma::cx_vec const &w) {
      return dot(block, v, w);
    };
    {
      Log("real time");
      // Real time evolution
      double t = 1.2345;
      auto H = matrixC(ops, block);
      arma::cx_vec psi_ex = expmat(complex(0.0, -1.0 * t) * H) * psi0.vector();
      arma::cx_vec psi = psi0c.vectorC();
      exp_sym_v(multC, dot_fC, psi, complex(0, -t));

      // XDIAG_SHOW(norm(psi - psi_ex));
      REQUIRE(norm(psi - psi_ex) < 1e-6);
    }

    {
      Log("imaginary time");
      // Imaginary time evolution
      double t = 1.2345;
      auto H = matrix(ops, block);
      arma::vec psi_ex = expmat(-t * H) * psi0.vector();
      arma::vec psi = psi0.vector();
      exp_sym_v(mult, dot_f, psi, -t);
      // XDIAG_SHOW(norm(psi - psi_ex));
      REQUIRE(norm(psi - psi_ex) < 1e-6);
    }

    {
      Log("complex time");
      complex t(-1.2345, 3.2145);
      auto H = matrixC(ops, block);
      arma::cx_vec psi_ex = expmat(t * H) * psi0.vector();
      arma::cx_vec psi = psi0c.vectorC();
      exp_sym_v(multC, dot_fC, psi, t);
      // XDIAG_SHOW(norm(psi - psi_ex));
      REQUIRE(norm(psi - psi_ex) < 1e-6);
    }
  }

} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
