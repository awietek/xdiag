// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include "../blocks/electron/testcases_electron.hpp"
#include <xdiag/algebra/algebra.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/random_state.hpp>

#include <iostream>
#include <set>

using namespace xdiag;

TEST_CASE("random_state", "[states]") try {
  using namespace xdiag::testcases::electron;

  Log.out("random state Spinhalf distinction test");
  double first_r = 0.;
  complex first_c = 0.;

  // Test whether random states from different blocks are different
  for (int nsites = 6; nsites <= 8; ++nsites) {
    Log("N={}", nsites);
    auto irreps = get_cyclic_group_irreps(nsites);

    for (int nup = 0; nup <= nsites; ++nup) {
      for (auto irrep : irreps) {
        auto block = Spinhalf(nsites, nup, irrep);

        if (block.size() > 3) {
          // XDIAG_SHOW(irrep);
          // XDIAG_SHOW(block);

          auto state_real = State(block, true);
          auto state_cplx = State(block, false);
          fill(state_real, RandomState());
          fill(state_cplx, RandomState());
          // XDIAG_SHOW(state_real);
          // XDIAG_SHOW(state_real.vector());
          // XDIAG_SHOW(state_real.vector());
          // XDIAG_SHOW(state_cplx.vector());
          if (first_r == 0.) {
            first_r = state_real.vector(false)(0);
          } else {
            REQUIRE(std::abs(state_real.vector(false)(0) - first_r) > 1e-12);
          }
          if (first_c == 0.) {
            first_c = state_cplx.vectorC(false)(0);
          } else {
            REQUIRE(std::abs(state_cplx.vectorC(false)(0) - first_c) > 1e-12);
          }
        }
      }
    }
  }
#ifdef _OPENMP

  // Check whether result with multiple threads is the same as with a single
  // thread
  Log.out("random state Spinhalf omp test");

  auto block = Spinhalf(4);
  for (int seed = 0; seed < 10; ++seed) {
    auto state = State(block, true);
    fill(state, RandomState(seed));
    auto state_cplx = State(block, false);
    fill(state_cplx, RandomState(false));

    omp_set_num_threads(1);
    auto state2 = State(block, true);
    fill(state2, RandomState(seed));
    auto state2_cplx = State(block, false);
    fill(state2_cplx, RandomState(false));

    REQUIRE(arma::norm(state.vector() - state2.vector()) < 1e-12);
    REQUIRE(arma::norm(state_cplx.vectorC() - state2_cplx.vectorC()) < 1e-12);
  }
#endif
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
