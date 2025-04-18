// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"
#include <mpi.h>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/spinhalf_distributed.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/utils/logger.hpp>

#include "../spinhalf/testcases_spinhalf.hpp"

using namespace xdiag;

static void test_onsite(std::string op1, std::string op12) {
  for (int nsites = 2; nsites < 5; ++nsites) {
    for (int nup = 0; nup <= nsites; ++nup) {
      auto b = SpinhalfDistributed(nsites, nup);
      for (int s = 0; s < nsites; ++s) {
        auto o1 = Op(op1, s);
        auto o12 = Op(op12, {s, s});
        auto r = random_state(b);
        auto v1 = apply(o1, apply(o1, r));
        auto v12 = apply(o12, r);
        REQUIRE(isapprox(v1, v12));
      }
    }
  }
}

static void test_e0_nompi(int N, OpSum ops) {
  for (int nup = 0; nup <= N; ++nup) {
    auto block = Spinhalf(N, nup);
    auto block_mpi = SpinhalfDistributed(N, nup);

    auto H = matrix(ops, block);
    REQUIRE(H.is_hermitian(1e-7));

    arma::vec evals_mat;
    arma::eig_sym(evals_mat, H);
    double e0_mat = evals_mat(0);

    double e0_app = eigval0(ops, block_mpi);

    // Log("N: {}, nup: {}, e0 mat: {:+.10f}, e0 mpi: {:+.10f}", N, nup,
    // e0_mat,
    //     e0_app);
    REQUIRE(isapprox(e0_mat, e0_app));
  }
}

static void test_sz_sp_sm_energy(int N, OpSum const &ops) {
  for (int nup = 0; nup <= N; ++nup) {
    auto block = Spinhalf(N, nup);
    auto block_mpi = SpinhalfDistributed(N, nup);

    auto [e0_s, gs_s] = eig0(ops, block);
    auto [e0_p, gs_p] = eig0(ops, block_mpi);
    REQUIRE(std::abs(e0_s - e0_p) < 1e-8);

    for (int i = 0; i < N; ++i) {
      // Log("N: {}, nup: {}, i: {}", N, nup, i);
      auto op = Op("Sz", i);
      double exp_s = inner(op, gs_s);
      double exp_p = inner(op, gs_p);
      REQUIRE(isapprox(exp_s, exp_p));

      if (nup < N - 1) {
        op = Op("S+", i);
        auto sz_i_gs_s = zero_state(Spinhalf(N, nup + 1));
        apply(op, gs_s, sz_i_gs_s);
        double dot_s = dot(sz_i_gs_s, sz_i_gs_s);
        auto sz_i_gs_p = zero_state(SpinhalfDistributed(N, nup + 1));
        apply(op, gs_p, sz_i_gs_p);
        double dot_p = dot(sz_i_gs_p, sz_i_gs_p);

        // Log("dot_s+: {:+.10f}, dot_p: {:+.10f}", dot_s, dot_p);
        REQUIRE(isapprox(dot_s, dot_p));
      }

      if (nup > 0) {
        op = Op("S-", i);
        auto sz_i_gs_s = zero_state(Spinhalf(N, nup - 1));
        apply(op, gs_s, sz_i_gs_s);
        double dot_s = dot(sz_i_gs_s, sz_i_gs_s);
        auto sz_i_gs_p = zero_state(SpinhalfDistributed(N, nup - 1));
        apply(op, gs_p, sz_i_gs_p);
        double dot_p = dot(sz_i_gs_p, sz_i_gs_p);

        // Log("dot_s-: {:+.10f}, dot_p: {:+.10f}", dot_s, dot_p);
        REQUIRE(isapprox(dot_s, dot_p));
      }
    }
  }
}

static void test_sz_sp_sm_commutators(int nsites) {
  for (int nup = 1; nup < nsites - 1; ++nup) {
    // Log("N: {}, nup: {}", nsites, nup);
    auto block = SpinhalfDistributed(nsites, nup);
    auto block_p = SpinhalfDistributed(nsites, nup + 1);
    auto block_m = SpinhalfDistributed(nsites, nup - 1);

    for (int i = 0; i < nsites; ++i)
      for (int j = 0; j < nsites; ++j) {

        auto sp_i = Op("S+", i);
        auto sm_i = Op("S-", i);
        auto sp_j = Op("S+", j);
        auto sm_j = Op("S-", j);
        auto sz_i = Op("Sz", i);

        auto rvec = random_state(block);
        auto sm_rvec = zero_state(block_m);
        auto sp_sm_rvec = zero_state(block);
        apply(sm_j, rvec, sm_rvec);
        apply(sp_i, sm_rvec, sp_sm_rvec);

        auto sp_rvec = zero_state(block_p);
        auto sm_sp_rvec = zero_state(block);
        apply(sp_i, rvec, sp_rvec);
        apply(sm_j, sp_rvec, sm_sp_rvec);

        double nrm = dot(rvec, sp_sm_rvec - sm_sp_rvec);

        // Check [S+_i, S-_j] = 2 Sz_i  delta_ij
        if (i == j) {
          double exp = 2 * inner(sz_i, rvec);
          REQUIRE(isapprox(nrm, exp));
        } else {
          REQUIRE(isapprox(nrm, 0.));
        }
      }
  }
}

TEST_CASE("spinhalf_distributed_apply", "[spinhalf_distributed]") try {

  using namespace xdiag::testcases::spinhalf;
  {
    Log("SpinhalfDistributed: manual N=6 spin chain test");

    OpSum ops;
    std::string type = "SzSz";
    ops += "J" * Op(type, {0, 1});
    ops += "J" * Op(type, {1, 2});
    ops += "J" * Op(type, {2, 3});
    ops += "J" * Op(type, {3, 4});
    ops += "J" * Op(type, {4, 5});
    ops += "J" * Op(type, {5, 0});

    // postfix ops
    ops += "J" * Op("Exchange", {0, 1});
    ops += "J" * Op("Exchange", {1, 2});

    // mixed ops
    ops += "J" * Op("SdotS", {2, 3});
    ops += "J" * Op("SdotS", {1, 4});
    ops += "J" * Op("SdotS", {0, 3});

    // Prefix ops
    ops += "J" * Op("SdotS", {3, 4});
    ops += "J" * Op("SdotS", {4, 5});
    ops += "J2" * Op("SdotS", {4, 5});

    ops["J"] = 1;
    ops["J2"] = 0.1;
    test_e0_nompi(6, ops);
  }

  Log("SpinhalfDistributed: Heisenberg chain apply test, J=1.0, N=2,..,8");
  for (int N = 2; N <= 8; ++N) {
    auto ops = HBchain(N, 1.0);
    test_e0_nompi(N, ops);
  }

  Log("SpinhalfDistributed: Heisenberg alltoall apply test, N=2,..,8");
  for (int N = 2; N <= 8; ++N) {
    auto ops = HB_alltoall(N);
    test_e0_nompi(N, ops);
  }

  // Test S+, S-, Sz operators
  Log("SpinhalfDistributed: Heisenberg chain Sz,S+,S- energy test, N=2,..,6");
  for (int N = 2; N <= 6; N += 2) {
    auto ops = HBchain(N, 1.0);
    test_sz_sp_sm_energy(N, ops);
  }

  // Test S+, S-, Sz operators
  Log("SpinhalfDistributed: Sz,S+,S- commutator test, N=2,..,8");
  for (int N = 2; N <= 8; ++N) {
    test_sz_sp_sm_commutators(N);
  }

  test_onsite("Sz", "SzSz");

  for (int nsites = 2; nsites < 6; ++nsites) {
    for (int nup = 1; nup < nsites; ++nup) {
      auto b = SpinhalfDistributed(nsites, nup);
      for (int s = 0; s < nsites; ++s) {

        auto r = random_state(b);

        // Exchange
        auto sp = Op("S+", s);
        auto sm = Op("S-", s);

        auto spsmr = apply(sp, apply(sm, r));
        auto smspr = apply(sm, apply(sp, r));
        auto r1 = 0.5 * (spsmr + smspr);
        auto r2 = apply(Op("Exchange", {s, s}), r);
        REQUIRE(isapprox(r1, r2));

        // SdotS
        auto szszr = apply(Op("SzSz", {s, s}), r);
        r1 = 0.5 * (spsmr + smspr) + szszr;
        r2 = apply(Op("SdotS", {s, s}), r);
        REQUIRE(isapprox(r1, r2));
      }
    }
  }

} catch (Error const &e) {
  error_trace(e);
}
