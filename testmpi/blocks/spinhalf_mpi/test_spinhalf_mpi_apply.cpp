#include <mpi.h>
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
    REQUIRE(lila::close(e0_mat, e0_app, 1e-6, 1e-6));
    REQUIRE(lila::close(e0_mat, e0_app_cplx, 1e-6, 1e-6));
  }
}

void test_sz_sp_sm_energy(BondList bonds, Couplings couplings) {
  int N = bonds.n_sites();
  for (int nup = 0; nup <= N; ++nup) {
    auto block = Spinhalf<uint32_t>(N, nup);
    auto block_mpi = SpinhalfMPI<uint32_t>(N, nup);

    
    auto [e0_s, gs_s] = GroundstateReal(bonds, couplings, block);
    auto [e0_p, gs_p] = GroundstateReal(bonds, couplings, block_mpi);
    REQUIRE(std::abs(e0_s - e0_p) < 1e-8);
    
    for (int i = 0; i < N; ++i) {
      auto bond = Bond("SZ", i);
      double exp_s = Inner(bond, gs_s);
      double exp_p = Inner(bond, gs_p);

      // LogMPI.out("N: {}, n_up: {}, sz_s: {:+.10f}, sz_p: {:+.10f}", N, nup,
      //            exp_s, exp_p);

      REQUIRE(lila::close(exp_s, exp_p));

      if (nup < N - 1) {
        bond = Bond("S+", i);
        auto sz_i_gs_s = Apply(bond, gs_s);
        double dot_s = Dot(sz_i_gs_s, sz_i_gs_s);
        auto sz_i_gs_p = Apply(bond, gs_p);
        double dot_p = Dot(sz_i_gs_p, sz_i_gs_p);

        LogMPI.out("N: {}, n_up: {}, i: {} dot_s+: {:+.10f}, dot_p+: {:+.10f}",
                   N, nup, i, dot_s, dot_p);
        REQUIRE(lila::close(dot_s, dot_p));
      }

      if (nup > 0) {
        bond = Bond("S-", i);
        auto sz_i_gs_s = Apply(bond, gs_s);
        double dot_s = Dot(sz_i_gs_s, sz_i_gs_s);
        auto sz_i_gs_p = Apply(bond, gs_p);
        double dot_p = Dot(sz_i_gs_p, sz_i_gs_p);

        LogMPI.out("N: {}, n_up: {}, i: {} dot_s-: {:+.10f}, dot_p-: {:+.10f}",
                   N, nup, i, dot_s, dot_p);
        REQUIRE(lila::close(dot_s, dot_p));
      }
    }
  }
}

void test_sz_sp_sm_commutators(int n_sites) {
  for (int nup = 1; nup < n_sites; ++nup) {
    LogMPI.out("N: {}, n_up: {}", n_sites, nup);
    auto block = SpinhalfMPI(n_sites, nup);
    auto block_s = Spinhalf(n_sites, nup);

    for (int i = 0; i < n_sites; ++i)
      for (int j = 0; j < n_sites; ++j) {

        auto sp_i = Bond("S+", i);
        auto sm_i = Bond("S-", i);
        auto sp_j = Bond("S+", j);
        auto sm_j = Bond("S-", j);
        auto sz_i = Bond("SZ", i);

        auto rvec = RandomState<double>(block);
        rvec.vector() /= lila::Norm(rvec.vector());
        auto vec1 = Apply(sp_i, Apply(sm_j, rvec));
        auto vec2 = Apply(sm_j, Apply(sp_i, rvec));
        auto comm = vec1.vector() - vec2.vector();
        // LogMPI.out("i: {} j: {} comm: {}", i, j, lila::Norm(comm));

        auto rvec_s = RandomState<double>(block_s);
        rvec_s.vector() /= lila::Norm(rvec_s.vector());

        auto vec1_s = Apply(sp_i, Apply(sm_j, rvec_s));
        auto vec2_s = Apply(sm_j, Apply(sp_i, rvec_s));
        auto comm_s = vec1_s.vector() - vec2_s.vector();

        // LilaPrint(rvec.vector());
        // LilaPrint(vec1.vector());
        // LilaPrint(vec2.vector());

        // LilaPrint(rvec_s.vector());
        // LilaPrint(vec1_s.vector());
        // LilaPrint(vec2_s.vector());

        // check [S^+_i, S^-_j] = 2 S^z_i \delta_{ij}
        if (i == j) {
          auto vec3 = Apply(sz_i, rvec);
          REQUIRE(lila::close(2.0 * vec3.vector(), comm));
        } else {
          REQUIRE(lila::close(lila::Norm(comm), 0.0));
        }
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

  LogMPI.out("SpinhalfMPI: Heisenberg chain apply test, J=1.0, N=2,..,8");
  for (int N = 2; N <= 8; ++N) {
    auto [bonds, couplings] = HBchain(N, 1.0);
    test_e0_nompi(bonds, couplings);
  }

  LogMPI.out("SpinhalfMPI: Heisenberg alltoall apply test, N=2,..,8");
  for (int N = 2; N <= 8; ++N) {
    auto [bonds, couplings] = HB_alltoall(N);
    test_e0_nompi(bonds, couplings);
  }

  // Test S+, S-, Sz operators
  LogMPI.out("SpinhalfMPI: Heisenberg chain Sz,S+,S- test, N=2,..,8");
  for (int N = 6; N <= 6; ++N) {
    auto [bonds, couplings] = HBchain(N, 1.0);
    test_sz_sp_sm_energy(bonds, couplings);
  }

  // Test S+, S-, Sz operators
  LogMPI.out("SpinhalfMPI: Sz,S+,S- commutator test, N=2,..,8");
  for (int N = 2; N <= 8; ++N) {
    test_sz_sp_sm_commutators(N);
  }
}
