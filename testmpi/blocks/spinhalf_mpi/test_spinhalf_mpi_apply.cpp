#include "../../catch.hpp"

#include <hydra/allmpi.h>

#include "../../../test/blocks/spinhalf/testcases_spinhalf.h"

using namespace hydra;

void test_e0_nompi(BondList bonds, Couplings couplings) {
  int N = bonds.n_sites();
  for (int nup = 0; nup <= N; ++nup) {
    auto block = Spinhalf<uint32_t>(N, nup);
    auto block_mpi = SpinhalfMPI<uint32_t>(N, nup);

    auto H = MatrixReal(bonds, couplings, block, block);
    REQUIRE(lila::close(H, lila::Herm(H)));

    auto evals_mat = lila::EigenvaluesSym(H);
    double e0_mat = evals_mat(0);
    double e0_app = E0Real(bonds, couplings, block_mpi);
    double e0_app_cplx = E0Cplx(bonds, couplings, block_mpi);

    LogMPI.out("N: {}, n_up: {}, e0 mat: {:+.10f}, e0 mpi: {:+.10f}", N, nup,
               e0_mat, e0_app);
    CHECK(lila::close(e0_mat, e0_app, 1e-6, 1e-6));
    CHECK(lila::close(e0_mat, e0_app_cplx, 1e-6, 1e-6));
  }
}

void test_sz_sp_sm_energy(BondList bonds, Couplings couplings) {
  int N = bonds.n_sites();
  for (int nup = 0; nup <= N; ++nup) {
    auto block = Spinhalf<uint32_t>(N, nup);
    auto block_mpi = SpinhalfMPI<uint32_t>(N, nup);

    auto [e0_s, gs_s] = GroundstateReal(bonds, couplings, block);
    auto [e0_p, gs_p] = GroundstateReal(bonds, couplings, block_mpi);

    for (int i = 0; i < N; ++i) {
      auto bond = Bond("SZ", i);
      double exp_s = Inner(bond, gs_s);
      double exp_p = Inner(bond, gs_p);

      LogMPI.out("N: {}, n_up: {}, sz_s: {:+.10f}, sz_p: {:+.10f}", N, nup,
                 exp_s, exp_p);

      REQUIRE(lila::close(exp_s, exp_p));
      

      // bond = Bond("S+", i);
      // auto sz_i_gs_s = Apply(bond, gs_s);
      // double dot_s = Dot(sz_i_gs_s, sz_i_gs_s);
      // auto sz_i_gs_p = Apply(bond, gs_p);
      // double dot_p = Dot(sz_i_gs_p, sz_i_gs_p);

      // LogMPI.out("N: {}, n_up: {}, sz_s: {:+.10f}, sz_p: {:+.10f}", N, nup,
      //            dot_s, dot_p);
    }
  }
}

TEST_CASE("spinhalf_mpi_apply", "[spinhalf]") {

  using namespace hydra::testcases::spinhalf;

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
    bonds << Bond("HB", "J2", {4, 5});

    Couplings couplings;
    couplings["J"] = 1;
    couplings["J2"] = 0.1;
    test_e0_nompi(bonds, couplings);
  }

  LogMPI.out("Spinhalf: Heisenberg chain apply test, J=1.0, N=2,..,8");
  for (int N = 2; N <= 8; ++N) {
    auto [bonds, couplings] = HBchain(N, 1.0);
    test_e0_nompi(bonds, couplings);
  }

  LogMPI.out("Spinhalf: Heisenberg alltoall apply test, N=2,..,8");
  for (int N = 2; N <= 8; ++N) {
    auto [bonds, couplings] = HB_alltoall(N);
    test_e0_nompi(bonds, couplings);
  }

  // // Test S+, S-, Sz operators
  // LogMPI.out("Spinhalf: Heisenberg chain Sz,S+,S- test, N=2,..,8");
  // for (int N = 6; N <= 6; ++N) {
  //   auto [bonds, couplings] = HBchain(N, 1.0);
  //   test_sz_sp_sm_energy(bonds, couplings);
  // }

}
