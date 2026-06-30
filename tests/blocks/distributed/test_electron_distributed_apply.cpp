// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <mpi.h>

#include <tests/blocks/tj/testcases_tj.hpp>
#include <tests/blocks/spinhalf/testcases_spinhalf.hpp>
#include <tests/blocks/electron/testcases_electron.hpp>
#include <tests/catch.hpp>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/blocks/distributed/tj_distributed.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/kernels/apply.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/states/apply.hpp>
#include <xdiag/states/inner.hpp>
#include <xdiag/states/dot.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/utils/logger.hpp>

#include <tests/blocks/distributed/compare_observables.hpp>


using namespace xdiag;

static void test_onsite(std::string op1, std::string op12) {
  for (int nsites = 2; nsites < 5; ++nsites) {
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto b = ElectronDistributed(nsites, nup, ndn);
        for (int s = 0; s < nsites; ++s) {
          auto r = random_state(b);
          auto x = apply(Op(op1, s), apply(Op(op1, s), r));
          auto y = apply(Op(op12, {s, s}), r);
          REQUIRE(isapprox(x, y));
        }
      }
    }
  }
}

static void test_regular_block_e0(OpSum ops, int nsites, int nup, int ndn) {
  auto block = Electron(nsites, nup, ndn);
  auto blockd = ElectronDistributed(nsites, nup, ndn);
  double e0 = eigval0(ops, block);
  double e0d = eigval0(ops, blockd);
  REQUIRE(isapprox(e0, e0d));
}

TEST_CASE("electron_distributed_apply", "[electron_distributed]") {
  using namespace xdiag::testcases::electron;

  Log("electron_distributed_apply test");
  int N = 4;

  // Play-around test
  for (int nup = 0; nup <= N; ++nup) {
    for (int ndn = 0; ndn <= N; ++ndn) {
      auto block = ElectronDistributed(N, nup, ndn);
      OpSum ops;
      for (int i = 0; i < N; ++i) {
        ops += "Jz" * Op("SzSz", {i, (i + 1) % N});
        ops += "Jx" * Op("Exchange", {i, (i + 1) % N});
        ops += "Jxasym" * Op("ExchangeAsym", {i, (i + 1) % N});
        ops += "TDN" * Op("Hopdn", {i, (i + 1) % N});
        ops += "TUP" * Op("Hopup", {i, (i + 1) % N});
        ops += "TDNasym" * Op("HopdnAsym", {i, (i + 1) % N});
        ops += "TUPasym" * Op("HopupAsym", {i, (i + 1) % N});
      }
      ops["Jz"] = 1.32;
      // A complex Exchange/Hop is not Hermitian on its own; the Hermitian
      // combination is the real (symmetric) part plus the imaginary
      // (antisymmetric) part, matching the non-distributed Electron block.
      ops["Jx"] = 0.432;
      ops["Jxasym"] = complex(0., 0.576);
      ops["TDN"] = -0.1432;
      ops["TUP"] = -0.4321;
      ops["TDNasym"] = complex(0., .3576);
      ops["TUPasym"] = complex(0., .5763);

      ops += 8.34 * Op("HubbardU");

      double e0 = eigval0(ops, block);
      auto block2 = Electron(N, nup, ndn);
      double e02 = eigval0(ops, block2);
      // Log("{} {} {:.12f} {:.12f}", nup, ndn, e0, e02);
      REQUIRE(isapprox(e0, e02));
    }
  }

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (real)
  for (int nsites = 3; nsites < 7; ++nsites) {

    Log("electron_distributed_apply: Hubbard random all-to-all test (real), N: "
        "{}",
        nsites);
    OpSum ops = freefermion_alltoall(nsites);
    ops["U"] = 5.0;

    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        test_regular_block_e0(ops, nsites, nup, ndn);
      }
    }
  }

  // /////////////////
  // Test Fermion all to all, free fermions (cplx, up/dn different)
  for (int nsites = 3; nsites < 7; ++nsites) {
    Log("electron_apply: Hubbard random all-to-all test (cplx), N: {}", nsites);
    OpSum ops = freefermion_alltoall_complex_updn(nsites);
    ops["U"] = 5.0;

    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        test_regular_block_e0(ops, nsites, nup, ndn);
      }
    }
  }

  ////////////////////////////
  // Henry's MATLAB code test (tests Heisenberg terms)
  Log("electron_apply: U-hopping-HB apply of Henry's Matlab code");
  {
    auto [ops, eigs_correct] = randomAlltoAll4NoU();

    int nsites = 4;
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        test_regular_block_e0(ops, nsites, nup, ndn);
      }
    }
  }

  ////////////////////////////////
  // all-to-all complex
  for (int nsites = 3; nsites <= 5; ++nsites) {
    Log.out("electron_apply: random all-to-all complex exchange test N={}",
            nsites);
    auto ops = xdiag::testcases::tj::tj_alltoall_complex(nsites);
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        test_regular_block_e0(ops, nsites, nup, ndn);
      }
    }
  }

  test_onsite("Sz", "SzSz");
  test_onsite("Ntot", "NtotNtot");
  test_onsite("Nupdn", "NupdnNupdn");

  // Test corrs
  for (int nsites = 2; nsites < 5; ++nsites) {
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto b = ElectronDistributed(nsites, nup, ndn);
        auto r = random_state(b);

        auto basis = ElectronDistributed(nsites, nup, ndn);

        for (int i = 0; i < nsites; ++i) {
          auto a = apply(Op("Nup", i), apply(Op("Ndn", i), r));
          auto b = apply(Op("Nupdn", i), r);
          REQUIRE(isapprox(a, b));

          a = apply(Op("Nup", i), r) + apply(Op("Ndn", i), r);
          b = apply(Op("Ntot", i), r);
          REQUIRE(isapprox(a, b));

          a = 0.5 * (apply(Op("Nup", i), r) - apply(Op("Ndn", i), r));
          b = apply(Op("Sz", i), r);
          REQUIRE(isapprox(a, b));

          for (int j = 0; j < nsites; ++j) {


            a = apply(Op("SzSz", {i, j}), r);
            b = apply(Op("Sz", i), apply(Op("Sz", j), r));
            REQUIRE(isapprox(a, b));

            a = apply(Op("NtotNtot", {i, j}), r);
            b = apply(Op("Ntot", i), apply(Op("Ntot", j), r));
            REQUIRE(isapprox(a, b));

            a = apply(Op("NupdnNupdn", {i, j}), r);
            b = apply(Op("Nupdn", i), apply(Op("Nupdn", j), r));
            REQUIRE(isapprox(a, b));

            a = apply(Op("NupNdn", {i, j}), r);
            b = apply(Op("Nup", i), apply(Op("Ndn", j), r));
            REQUIRE(isapprox(a, b));

            a = apply(Op("NdnNup", {i, j}), r);
            b = apply(Op("Ndn", i), apply(Op("Nup", j), r));
            REQUIRE(isapprox(a, b));

            a = apply(Op("NupNup", {i, j}), r);
            b = apply(Op("Nup", i), apply(Op("Nup", j), r));
            REQUIRE(isapprox(a, b));

            a = apply(Op("NdnNdn", {i, j}), r);
            b = apply(Op("Ndn", i), apply(Op("Ndn", j), r));
            REQUIRE(isapprox(a, b));
          }
        }
      }
    }
  }
}

TEST_CASE("electron_distributed_observables", "[electron_distributed]") {
  Log("electron_distributed_observables: expect / correlation_matrix vs "
      "Electron");
  int N = 4;
  // Generic Hermitian electron Hamiltonian (real symmetric + imaginary
  // antisymmetric hop/exchange + Hubbard U) -> non-degenerate ground state.
  OpSum ops;
  for (int i = 0; i < N; ++i) {
    ops += 1.10 * Op("SzSz", {i, (i + 1) % N});
    ops += 0.43 * Op("Exchange", {i, (i + 1) % N});
    ops += complex(0.0, 0.31) * Op("ExchangeAsym", {i, (i + 1) % N});
    ops += -0.97 * Op("Hop", {i, (i + 1) % N});
    ops += complex(0.0, 0.22) * Op("HopAsym", {i, (i + 1) % N});
  }
  ops += 3.0 * Op("HubbardU");
  std::vector<std::string> onesite = {"Nup", "Ndn", "Ntot", "Sz"};
  std::vector<std::pair<std::string, std::string>> twosite = {
      {"Sz", "Sz"}, {"Ntot", "Ntot"}, {"Cdagup", "Cup"}, {"Cdagdn", "Cdn"}};
  std::vector<std::tuple<int, int>> nup_ndn = {{2, 1}, {2, 2}, {1, 2}};
  for (auto [nup, ndn] : nup_ndn) {
    compare_observables(ops, Electron(N, nup, ndn),
                        ElectronDistributed(N, nup, ndn), onesite, twosite);
  }
}
