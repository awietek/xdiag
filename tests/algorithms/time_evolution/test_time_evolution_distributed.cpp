#include "../../catch.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algorithms/time_evolution/time_evolution.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/blocks/tj_distributed.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/product_state.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/close.hpp>

TEST_CASE("time_evolution_distributed", "[time_evolution]") try {
  using namespace xdiag;

  Log("Test time_evolution_distributed");
  int L = 3;
  int n_sites = L * L;

  // Create square lattice t-J model
  OpSum ops;
  for (int x = 0; x < L; ++x) {
    for (int y = 0; y < L; ++y) {
      int nx = (x + 1) % L;
      int ny = (y + 1) % L;

      int site = y * L + x;
      int right = y * L + nx;
      int top = ny * L + x;
      ops += "T" * Op("HOP", {site, right});
      ops += "JZ" * Op("TJSZSZ", {site, right});
      ops += "JEX" * Op("EXCHANGE", {site, right});

      ops += "T" * Op("HOP", {site, top});
      ops += "JZ" * Op("TJSZSZ", {site, top});
      ops += "JEX" * Op("EXCHANGE", {site, top});
    }
  }
  ops["T"] = 1.0 + 0.2i;
  ops["JZ"] = 0.4;
  ops["JEX"] = 0.3 + 1.23i;

  // Create initial state
  auto pstate = ProductState();
  for (int x = 0; x < L; ++x) {
    for (int y = 0; y < L; ++y) {
      if (((x + y) % 2) == 0) {
        pstate.push_back("Dn");
      } else {
        pstate.push_back("Up");
      }
    }
  }
  pstate[n_sites / 2] = "Emp";
  auto block = tJ(n_sites, n_sites / 2, n_sites / 2);
  auto blockd = tJDistributed(n_sites, n_sites / 2, n_sites / 2);

  auto psi_0 = State(block, false);
  auto psi_0d = State(blockd, false);
  fill(psi_0, pstate);
  fill(psi_0d, pstate);

  auto H_psi_0 = State(block, false);
  apply(ops, psi_0, H_psi_0);
  auto H_psi_0d = State(blockd, false);
  apply(ops, psi_0d, H_psi_0d);

  for (int s = 0; s < n_sites; ++s) {
    auto n = innerC(Op("NTOT", s), H_psi_0);
    auto nd = innerC(Op("NTOT", s), H_psi_0d);
    // Log("i {} {} {}", s, n, nd);
    REQUIRE(close(n, nd));
  }
  // Log("\n");

  arma::vec times = arma::logspace(-1, 1, 3);
  double tol = 1e-12;
  for (auto time : times) {
    auto psi = time_evolve(ops, psi_0, time, tol);
    auto psid = time_evolve(ops, psi_0d, time, tol);
    for (int s = 0; s < n_sites; ++s) {
      auto n = innerC(Op("NTOT", s), psi);
      auto nd = innerC(Op("NTOT", s), psid);
      Log("{} {} {} {:.6f}", s, n, nd, time);
      REQUIRE(std::abs(n - nd) < 1e-6);
    }
    Log("\n");
  }

} catch (xdiag::Error e) {
  error_trace(e);
}
