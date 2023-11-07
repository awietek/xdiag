#include "../../catch.hpp"

#include <iostream>

#include "../spinhalf/testcases_spinhalf.h"
#include "../tj/testcases_tj.h"

#include <hydra/algorithms/sparse_diag.h>
#include <hydra/blocks/spinhalf/spinhalf_matrix.h>
#include <hydra/blocks/tj/tj_apply.h>
#include <hydra/blocks/tj/tj_matrix.h>
#include <hydra/blocks/tj_distributed/tj_distributed_apply.h>
#include <hydra/utils/close.h>
#include <hydra/utils/print_macro.h>

using namespace hydra;

TEST_CASE("tj_distributed_apply", "[tj_distributed]") try {
  using namespace hydra::testcases::tj;

  Log("tj_distributed_apply test");
  int N = 2;
  auto block = tJDistributed(N, N / 2, N / 2);
  BondList bonds;
  for (int i = 0; i < N - 1; ++i) {
    bonds << Bond("ISING", "Jz", {i, (i + 1) % N});
    bonds << Bond("EXCHANGE", "Jx", {i, (i + 1) % N});
  }
  bonds["Jz"] = 1.0;
  bonds["Jx"] = complex(1, 1);

  double e0 = eigval0(bonds, block);
  HydraPrint(e0);

  auto block2 = tJ(N, N / 2, N / 2);
  double e02 = eigval0(bonds, block2);
  REQUIRE(close(e0, e02));

  auto block3 = Spinhalf(N, N / 2);
  double e03 = eigval0(bonds, block3);
  REQUIRE(close(e0, e03));

  Log("tj_distributed: HB all-to-all comparison");
  for (int n_sites = 2; n_sites < 7; ++n_sites) {
    Log("N: {}", n_sites);
    auto bonds = testcases::spinhalf::HB_alltoall(n_sites);

    for (int nup = 0; nup <= n_sites; ++nup) {
      auto block = Spinhalf(n_sites, nup);
      auto block_tJ = tJDistributed(n_sites, nup, n_sites - nup);
      double e0_spinhalf = eigval0(bonds, block);
      double e0 = eigval0(bonds, block_tJ);
      // Log("{} {}", e0_spinhalf, e0);
      REQUIRE(close(e0_spinhalf, e0));
    }
  }

  Log.out("tj_distributed: HB triangular N=12 complex exchange");
  int n_sites = 12;
  std::vector<double> etas = {0.00, 0.01, 0.02,
                              0.03, 0.04, 0.05}; // dont change etas :-)
  for (auto eta : etas) {
    Log("eta: {}", eta);
    for (int nup = 0; nup <= n_sites; ++nup) {
      auto [bonds, e00] = testcases::spinhalf::triangular_12_complex(nup, eta);

      auto block = Spinhalf(n_sites, nup);
      auto block_tJ = tJDistributed(n_sites, nup, n_sites - nup);
      double e0_spinhalf = eigval0(bonds, block);
      double e0 = eigval0(bonds, block_tJ);
      // Log("{} {} {} {}", nup, e0_spinhalf, e0, e00);
      REQUIRE(close(e0_spinhalf, e0));
      if (nup == 6) {
        REQUIRE(std::abs(e0 - e00) < 1e-8);
      }
    }
  }

} catch (std::exception const &e) {
  traceback(e);
}
