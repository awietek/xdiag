// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/kernels/apply.hpp>
#include <xdiag/linalg/time_evolution/time_evolve.hpp>
#include <xdiag/linalg/time_evolution/time_evolve_expokit.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/blocks/distributed/tj_distributed.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/states/apply.hpp>
#include <xdiag/states/inner.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/utils/logger.hpp>

TEST_CASE("time_evolution_distributed", "[time_evolution]") try {
  using namespace xdiag;

  Log("Test time_evolution_distributed");
  int L = 3;
  int nsites = L * L;

  // Create square lattice t-J model
  OpSum ops;
  for (int x = 0; x < L; ++x) {
    for (int y = 0; y < L; ++y) {
      int nx = (x + 1) % L;
      int ny = (y + 1) % L;

      int site = y * L + x;
      int right = y * L + nx;
      int top = ny * L + x;
      ops += "T" * Op("Hop", {site, right});
      ops += "Tasym" * Op("HopAsym", {site, right});
      ops += "JZ" * Op("tJSzSz", {site, right});
      ops += "JEX" * Op("Exchange", {site, right});
      ops += "JEXasym" * Op("ExchangeAsym", {site, right});

      ops += "T" * Op("Hop", {site, top});
      ops += "Tasym" * Op("HopAsym", {site, top});
      ops += "JZ" * Op("tJSzSz", {site, top});
      ops += "JEX" * Op("Exchange", {site, top});
      ops += "JEXasym" * Op("ExchangeAsym", {site, top});
    }
  }
  ops["T"] = 1.0;
  ops["Tasym"] = complex(0.0, 0.2);
  ops["JZ"] = 0.4;
  ops["JEX"] = 0.3;
  ops["JEXasym"] = complex(0.0, 1.23);

  // Create initial state
  auto pstate = ProductState();
  for (int x = 0; x < L; ++x) {
    for (int y = 0; y < L; ++y) {
      if (((x + y) % 2) == 0) {
        pstate.push_back(2);
      } else {
        pstate.push_back(1);
      }
    }
  }
  pstate[nsites / 2] = 0;

  auto block = tJ(nsites, nsites / 2, nsites / 2);
  auto blockd = tJDistributed(nsites, nsites / 2, nsites / 2);

  auto psi_0 = State(block, false);
  auto psi_0d = State(blockd, false);
  fill(psi_0, pstate);
  fill(psi_0d, pstate);

  auto H_psi_0 = State(block, false);
  apply(ops, psi_0, H_psi_0);
  auto H_psi_0d = State(blockd, false);
  apply(ops, psi_0d, H_psi_0d);
  for (int s = 0; s < nsites; ++s) {
    auto n = innerC(Op("Ntot", s), H_psi_0);
    auto nd = innerC(Op("Ntot", s), H_psi_0d);
    // Log("i {} {} {}", s, n, nd);
    REQUIRE(isapprox(n, nd));
  }
  // Log("\n");

  arma::vec times = arma::logspace(-1, 1, 3);
  double tol = 1e-12;
  for (auto time : times) {
    // Log("time: {}", time);
    auto psi = time_evolve(ops, psi_0, time, tol);
    auto psid = time_evolve(ops, psi_0d, time, tol);
    for (int s = 0; s < nsites; ++s) {
      auto n = innerC(Op("Ntot", s), psi);
      auto nd = innerC(Op("Ntot", s), psid);
      // Log("{} {} {} {:.6f}", s, n, nd, time);
      REQUIRE(std::abs(n - nd) < 1e-6);
    }
    // Log("\n");
  }

} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}
