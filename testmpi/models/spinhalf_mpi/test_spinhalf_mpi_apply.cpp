#include "../../catch.hpp"

#include <hydra/allmpi.h>

#include "../../../test/models/spinhalf/testcases_spinhalf.h"

using namespace hydra;

void test_e0_nompi(BondList bonds, Couplings couplings) {
  int N = bonds.n_sites();
  for (int nup = 0; nup <= N; ++nup) {
    auto block = Spinhalf<uint32>(N, nup);
    auto block_mpi = SpinhalfMPI<uint32>(N, nup);

    auto H = MatrixReal(bonds, couplings, block, block);
    REQUIRE(lila::close(H, lila::Herm(H)));

    auto evals_mat = lila::EigenvaluesSym(H);
    double e0_mat = evals_mat(0);
    double e0_app = E0Real(bonds, couplings, block_mpi);
    LogMPI.out("N: {}, n_up: {}, e0 mat: {:+.10f}, e0 mpi: {:+.10f}", N, nup, e0_mat,
               e0_app);
    REQUIRE(lila::close(e0_mat, e0_app, 1e-6, 1e-6));
  }
}

TEST_CASE("spinhalf_mpi_apply", "[spinhalf]") {

  {
    LogMPI.out("SpinhalfMPI: manual N=6 spin chain test");

    BondList bonds;
    std::string type = "ISING";
    bonds << Bond(type, "J", {0, 1});
    bonds << Bond(type, "J", {1, 2});
    bonds << Bond(type, "J", {2, 3});
    bonds << Bond(type, "J", {3, 4});
    bonds << Bond(type, "J", {4, 5});
    bonds << Bond(type, "J", {5, 0});

    // postfix bonds
    bonds << Bond("HB", "J", {0, 1});
    bonds << Bond("HB", "J", {1, 2});

    // mixed bonds
    bonds << Bond("HB", "J", {2, 3});
    bonds << Bond("HB", "J", {1, 4});
    bonds << Bond("HB", "J", {0, 3});

    // Prefix bonds
    bonds << Bond("HB", "J", {3, 4});
    bonds << Bond("HB", "J", {4, 5});

    Couplings couplings;
    couplings["J"] = 1;
    test_e0_nompi(bonds, couplings);
  }

  using namespace hydra::testcases::spinhalf;

  LogMPI.out("Spinhalf: Heisenberg chain apply test, J=1.0, N=2,..,6");
  for (int N = 2; N <= 6; ++N) {
    auto [bonds, couplings] = HBchain(N, 1.0);
    test_e0_nompi(bonds, couplings);
  }

  LogMPI.out("Spinhalf: Heisenberg alltoall apply test, N=2,..,6");
  for (int N = 2; N <= 6; ++N) {
    auto [bonds, couplings] = HB_alltoall(N);
    test_e0_nompi(bonds, couplings);
  }
}
