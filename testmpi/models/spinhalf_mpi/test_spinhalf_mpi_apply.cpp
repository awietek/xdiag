#include "../../catch.hpp"

#include <hydra/allmpi.h>

using namespace hydra;

void test_e0_nompi(BondList bonds, Couplings couplings) {
  int N = bonds.n_sites();
  for (int nup = 0; nup <= N; ++nup) {
    LogMPI.out("{} {}", N, nup);

    auto block = Spinhalf<uint32>(N, nup);
    auto block_mpi = SpinhalfMPI<uint32>(N, nup);

    auto H = MatrixReal(bonds, couplings, block, block);
    REQUIRE(lila::close(H, lila::Herm(H)));

    auto evals_mat = lila::EigenvaluesSym(H);
    double e0_mat = evals_mat(0);
    double e0_app = E0Real(bonds, couplings, block_mpi);
    LogMPI.out("{} {}", e0_mat, e0_app);
    REQUIRE(lila::close(e0_mat, e0_app, 1e-6, 1e-6));
  }
}

TEST_CASE("spinhalf_mpi_apply", "[spinhalf]") {
  LogMPI.out("SpinhalfMPI: Heisenberg chain apply test, J=1.0, N=2,..,6");
  BondList bonds;
  std::string type = "ISING";
  bonds << Bond(type, "J", {0, 1});
  bonds << Bond(type, "J", {1, 2});
  bonds << Bond(type, "J", {2, 3});
  bonds << Bond(type, "J", {3, 4});
  bonds << Bond(type, "J", {4, 5});
  bonds << Bond(type, "J", {5, 0});
  Couplings couplings;
  couplings["J"] = 1;
  test_e0_nompi(bonds, couplings);
}
