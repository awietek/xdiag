#include "../../catch.hpp"

#include <hydra/algebra/algebra.h>
#include <hydra/algebra/apply.h>
#include <hydra/algorithms/time_evolution/time_evolution.h>
#include <hydra/blocks/blocks.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/blocks/tj_distributed/tj_distributed.h>
#include <hydra/states/product_state.h>
#include <hydra/states/random_state.h>
#include <hydra/states/state.h>
#include <hydra/utils/close.h>

TEST_CASE("time_evolution_distributed", "[time_evolution]") try {
  using namespace hydra;

  int L = 3;
  int n_sites = L * L;

  // Create square lattice t-J model
  BondList bonds;
  for (int x = 0; x < L; ++x) {
    for (int y = 0; y < L; ++y) {
      int nx = (x + 1) % L;
      int ny = (y + 1) % L;

      int site = y * L + x;
      int right = y * L + nx;
      int top = ny * L + x;
      bonds << Bond("HOP", "T", {site, right});
      bonds << Bond("TJISING", "J", {site, right});
      bonds << Bond("HOP", "T", {site, top});
      bonds << Bond("TJISING", "J", {site, top});
    }
  }
  bonds["T"] = 1.0 + 0.2i;
  bonds["J"] = 0.4;

  // Create initial state
  auto pstate = ProductState();
  for (int x = 0; x < L; ++x) {
    for (int y = 0; y < L; ++y) {
      if (((x + y) % 2) == 0) {
        pstate << "Dn";
      } else {
        pstate << "Up";
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
  apply(bonds, psi_0, H_psi_0);
  auto H_psi_0d = State(blockd, false);
  apply(bonds, psi_0d, H_psi_0d);
  
  for (int s = 0; s < n_sites; ++s) {
    auto n = innerC(Bond("NUMBER", s), H_psi_0);
    auto nd = innerC(Bond("NUMBER", s), H_psi_0d);
    // Log("i {} {} {}", s, n, nd);
    REQUIRE(close(n, nd));
  }
  // Log("\n");

  arma::vec times = arma::logspace(-1, 1, 3);
  double tol = 1e-12;
  for (auto time : times) {
    auto psi = time_evolve(bonds, psi_0, time, tol);
    auto psid = time_evolve(bonds, psi_0d, -time, tol);
    for (int s = 0; s < n_sites; ++s) {
      auto n = innerC(Bond("NUMBER", s), psi);
      auto nd = innerC(Bond("NUMBER", s), psid);
      Log("{} {} {} {:.6f}", s, n, nd, time);
      REQUIRE(std::abs(n - nd) < 1e-6);
    }
    Log("\n");
  }

} catch (std::exception const &e) {
  hydra::traceback(e);
}
