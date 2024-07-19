#include "../../catch.hpp"
#include <mpi.h>

#include <xdiag/all.hpp>

#include "../spinhalf/testcases_spinhalf.hpp"

using namespace xdiag;

void test_e0_nompi(int N, OpSum ops) {
  for (int nup = 0; nup <= N; ++nup) {
    auto block = Spinhalf(N, nup);
    auto block_mpi = SpinhalfDistributed(N, nup);

    auto H = matrix(ops, block);
    REQUIRE(H.is_hermitian(1e-7));
    
    arma::vec evals_mat;
    arma::eig_sym(evals_mat, H);
    double e0_mat = evals_mat(0);

    double e0_app = eigval0(ops, block_mpi);

    Log("N: {}, n_up: {}, e0 mat: {:+.10f}, e0 mpi: {:+.10f}", N, nup, e0_mat,
        e0_app);
    // REQUIRE(close(e0_mat, e0_app));
  }
}

// void test_sz_sp_sm_energy(OpSum ops, Couplings couplings) {
//   int N = ops.n_sites();
//   for (int nup = 0; nup <= N; ++nup) {
//     auto block = Spinhalf<uint32_t>(N, nup);
//     auto block_mpi = SpinhalfMPI<uint32_t>(N, nup);

//     auto [e0_s, gs_s] = GroundstateReal(ops, couplings, block);
//     auto [e0_p, gs_p] = GroundstateReal(ops, couplings, block_mpi);
//     REQUIRE(std::abs(e0_s - e0_p) < 1e-8);

//     for (int i = 0; i < N; ++i) {
//       auto op = Op("SZ", i);
//       double exp_s = Inner(op, gs_s);
//       double exp_p = Inner(op, gs_p);

//       // LogMPI.out("N: {}, n_up: {}, sz_s: {:+.10f}, sz_p: {:+.10f}", N,
//       nup,
//       //            exp_s, exp_p);

//       REQUIRE(lila::close(exp_s, exp_p));

//       if (nup < N - 1) {
//         op = Op("S+", i);
//         auto sz_i_gs_s = Apply(op, gs_s);
//         double dot_s = Dot(sz_i_gs_s, sz_i_gs_s);
//         auto sz_i_gs_p = Apply(op, gs_p);
//         double dot_p = Dot(sz_i_gs_p, sz_i_gs_p);

//         LogMPI.out("N: {}, n_up: {}, i: {} dot_s+: {:+.10f}, dot_p+:
//         {:+.10f}",
//                    N, nup, i, dot_s, dot_p);
//         REQUIRE(lila::close(dot_s, dot_p));
//       }

//       if (nup > 0) {
//         op = Op("S-", i);
//         auto sz_i_gs_s = Apply(op, gs_s);
//         double dot_s = Dot(sz_i_gs_s, sz_i_gs_s);
//         auto sz_i_gs_p = Apply(op, gs_p);
//         double dot_p = Dot(sz_i_gs_p, sz_i_gs_p);

//         LogMPI.out("N: {}, n_up: {}, i: {} dot_s-: {:+.10f}, dot_p-:
//         {:+.10f}",
//                    N, nup, i, dot_s, dot_p);
//         REQUIRE(lila::close(dot_s, dot_p));
//       }
//     }
//   }
// }

// void test_sz_sp_sm_commutators(int n_sites) {
//   for (int nup = 1; nup < n_sites; ++nup) {
//     LogMPI.out("N: {}, n_up: {}", n_sites, nup);
//     auto block = SpinhalfMPI(n_sites, nup);
//     auto block_s = Spinhalf(n_sites, nup);

//     for (int i = 0; i < n_sites; ++i)
//       for (int j = 0; j < n_sites; ++j) {

//         auto sp_i = Op("S+", i);
//         auto sm_i = Op("S-", i);
//         auto sp_j = Op("S+", j);
//         auto sm_j = Op("S-", j);
//         auto sz_i = Op("SZ", i);

//         auto rvec = RandomState<double>(block);
//         rvec.vector() /= lila::Norm(rvec.vector());
//         auto vec1 = Apply(sp_i, Apply(sm_j, rvec));
//         auto vec2 = Apply(sm_j, Apply(sp_i, rvec));
//         auto comm = vec1.vector() - vec2.vector();
//         // LogMPI.out("i: {} j: {} comm: {}", i, j, lila::Norm(comm));

//         auto rvec_s = RandomState<double>(block_s);
//         rvec_s.vector() /= lila::Norm(rvec_s.vector());

//         auto vec1_s = Apply(sp_i, Apply(sm_j, rvec_s));
//         auto vec2_s = Apply(sm_j, Apply(sp_i, rvec_s));
//         auto comm_s = vec1_s.vector() - vec2_s.vector();

//         // LilaPrint(rvec.vector());
//         // LilaPrint(vec1.vector());
//         // LilaPrint(vec2.vector());

//         // LilaPrint(rvec_s.vector());
//         // LilaPrint(vec1_s.vector());
//         // LilaPrint(vec2_s.vector());

//         // check [S^+_i, S^-_j] = 2 S^z_i \delta_{ij}
//         if (i == j) {
//           auto vec3 = Apply(sz_i, rvec);
//           REQUIRE(lila::close(2.0 * vec3.vector(), comm));
//         } else {
//           REQUIRE(lila::close(lila::Norm(comm), 0.0));
//         }
//       }
//   }
// }

TEST_CASE("spinhalf_distributed_apply", "[spinhalf_distributed]") {

  using namespace xdiag::testcases::spinhalf;
  {
    Log("SpinhalfDistributed: manual N=6 spin chain test");

    OpSum ops;
    std::string type = "ISING";
    ops += Op(type, "J", {0, 1});
    ops += Op(type, "J", {1, 2});
    ops += Op(type, "J", {2, 3});
    ops += Op(type, "J", {3, 4});
    ops += Op(type, "J", {4, 5});
    ops += Op(type, "J", {5, 0});

    // postfix ops
    ops += Op("HB", "J", {0, 1});
    ops += Op("HB", "J", {1, 2});

    // mixed ops
    ops += Op("HB", "J", {2, 3});
    ops += Op("HB", "J", {1, 4});
    ops += Op("HB", "J", {0, 3});

    // Prefix ops
    ops += Op("HB", "J", {3, 4});
    ops += Op("HB", "J", {4, 5});
    ops += Op("HB", "J2", {4, 5});

    ops["J"] = 1;
    ops["J2"] = 0.1;
    test_e0_nompi(6, ops);
  }

  // Log("SpinhalfDistributed: Heisenberg chain apply test, J=1.0, N=2,..,8");
  // for (int N = 2; N <= 8; ++N) {
  //   auto ops = HBchain(N, 1.0);
  //   test_e0_nompi(N, ops);
  // }

  // Log("SpinhalfDistributed: Heisenberg alltoall apply test, N=2,..,8");
  // for (int N = 2; N <= 8; ++N) {
  //   auto ops = HB_alltoall(N);
  //   test_e0_nompi(N, ops);
  // }

  // // Test S+, S-, Sz operators
  // LogMPI.out("SpinhalfMPI: Heisenberg chain Sz,S+,S- test, N=2,..,8");
  // for (int N = 6; N <= 6; ++N) {
  //   auto [ops, couplings] = HBchain(N, 1.0);
  //   test_sz_sp_sm_energy(ops, couplings);
  // }

  // // Test S+, S-, Sz operators
  // LogMPI.out("SpinhalfMPI: Sz,S+,S- commutator test, N=2,..,8");
  // for (int N = 2; N <= 8; ++N) {
  //   test_sz_sp_sm_commutators(N);
  // }
}
