#include "../../catch.hpp"

#include <iostream>

#include "../spinhalf/testcases_spinhalf.hpp"
#include "../tj/testcases_tj.hpp"
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/states/create_state.hpp>

using namespace xdiag;

static void test_onsite(std::string op1, std::string op12) {
  for (int nsites = 2; nsites < 5; ++nsites) {
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites - nup; ++ndn) {
        auto b = tJDistributed(nsites, nup, ndn);
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
}

static void test_tjdistributed_e0_real(OpSum ops, int nsites, int nup, int ndn,
                                       double e0) {
  auto block = tJDistributed(nsites, nup, ndn);
  double e0c = eigval0(ops, block);
  REQUIRE(std::abs(e0 - e0c) < 1e-6);
}

TEST_CASE("tj_distributed_apply", "[tj_distributed]") try {
  using namespace xdiag::testcases::tj;

  Log("tj_distributed_apply test");
  int N = 4;

  // Play-around test
  for (int nup = 0; nup <= N; ++nup) {
    for (int ndn = 0; ndn <= N - nup; ++ndn) {
      auto block = tJDistributed(N, nup, ndn);
      OpSum ops;
      for (int i = 0; i < N; ++i) {
        ops += "Jz" * Op("SzSz", {i, (i + 1) % N});
        ops += "Jx" * Op("Exchange", {i, (i + 1) % N});
        ops += "TDN" * Op("Hopdn", {i, (i + 1) % N});
        ops += "TUP" * Op("Hopup", {i, (i + 1) % N});
      }
      ops["Jz"] = 1.32;
      ops["Jx"] = complex(.432, .576);
      ops["TDN"] = complex(-0.1432, .3576);
      ops["TUP"] = complex(-0.4321, .5763); // 2.104;

      double e0 = eigval0(ops, block);
      auto block2 = tJ(N, nup, ndn);
      double e02 = eigval0(ops, block2);
      // Log("{} {} {:.12f} {:.12f}", nup, ndn, e0, e02);

      REQUIRE(isapprox(e0, e02));
    }
  }

  Log("tj_distributed: HB all-to-all comparison");
  for (int nsites = 2; nsites < 8; ++nsites) {
    // Log("N: {}", nsites);
    auto ops = testcases::spinhalf::HB_alltoall(nsites);
    // XDIAG_SHOW(ops);
    for (int nup = 0; nup <= nsites; ++nup) {
      auto block = Spinhalf(nsites, nup);
      auto block_tJ = tJDistributed(nsites, nup, nsites - nup);
      double e0_spinhalf = eigval0(ops, block);
      double e0 = eigval0(ops, block_tJ);
      // Log("{} {}", e0_spinhalf, e0);
      REQUIRE(isapprox(e0_spinhalf, e0));
    }
  }

  {
    Log("tj_distributed: TJModel: six-site chain test, t=1.0, J=1.0, N=6");
    auto ops = tJchain(6, 1.0, 1.0);
    std::vector<std::tuple<int, int, double>> nup_ndn_e0 = {
        {0, 0, 0.0},         {0, 1, -2.0},        {0, 2, -2.96081311},
        {0, 3, -3.79610527}, {0, 4, -2.46081311}, {0, 5, -0.99999999},
        {0, 6, 1.500000000}, {1, 0, -2.0},        {1, 1, -3.61222054},
        {1, 2, -4.04537829}, {1, 3, -4.10768318}, {1, 4, -2.42705097},
        {1, 5, -0.49999999}, {2, 0, -2.96081311}, {2, 1, -4.04537829},
        {2, 2, -4.16447847}, {2, 3, -3.52922048}, {2, 4, -2.11803398},
        {3, 0, -3.79610527}, {3, 1, -4.10768318}, {3, 2, -3.52922048},
        {3, 3, -2.80277563}, {4, 0, -2.46081311}, {4, 1, -2.42705097},
        {4, 2, -2.11803398}, {5, 0, -0.99999999}, {5, 1, -0.49999999},
        {6, 0, 1.500000000}};
    for (auto [nup, ndn, e0] : nup_ndn_e0) {
      test_tjdistributed_e0_real(ops, 6, nup, ndn, e0);
    }
  }

  {
    Log("tj_distributed: TJModel: six-site chain test, t=1.0, J=0.0, N=6");
    auto ops = tJchain(6, 1.0, 0.0);
    std::vector<std::tuple<int, int, double>> nup_ndn_e0 = {
        {0, 0, 0.0},         {0, 1, -2.0},        {0, 2, -3.00000000},
        {0, 3, -4.00000000}, {0, 4, -2.99999999}, {0, 5, -2.00000000},
        {0, 6, 0.000000000}, {1, 0, -2.0},        {1, 1, -3.46410161},
        {1, 2, -3.99999999}, {1, 3, -3.46410161}, {1, 4, -1.99999999},
        {1, 5, 0.000000000}, {2, 0, -3.00000000}, {2, 1, -3.99999999},
        {2, 2, -3.46410161}, {2, 3, -1.99999999}, {2, 4, 0.000000000},
        {3, 0, -4.00000000}, {3, 1, -3.46410161}, {3, 2, -1.99999999},
        {3, 3, 0.000000000}, {4, 0, -2.99999999}, {4, 1, -1.99999999},
        {4, 2, 0.000000000}, {5, 0, -2.00000000}, {5, 1, 0.000000000},
        {6, 0, 0.000000000}};
    for (auto [nup, ndn, e0] : nup_ndn_e0) {
      test_tjdistributed_e0_real(ops, 6, nup, ndn, e0);
    }
  }

  {
    for (int nsites = 1; nsites < 8; ++nsites) {
      Log("tj_distributed: tj_alltoall random (real) N={}", nsites);
      auto ops = tj_alltoall(nsites);
      for (int nup = 0; nup <= nsites; ++nup) {
        for (int ndn = 0; ndn <= nsites - nup; ++ndn) {
          auto block = tJ(nsites, nup, ndn);
          auto block2 = tJDistributed(nsites, nup, ndn);

          double e0 = eigval0(ops, block);
          double e02 = eigval0(ops, block2);
          REQUIRE(isapprox(e0, e02));
        }
      }
    }
  }
  {
    for (int nsites = 1; nsites < 8; ++nsites) {
      Log("tj_distributed: tj_alltoall random (complex) N={}", nsites);
      auto ops = tj_alltoall_complex(nsites);
      for (int nup = 0; nup <= nsites; ++nup) {
        for (int ndn = 0; ndn <= nsites - nup; ++ndn) {

          auto block = tJ(nsites, nup, ndn);
          auto block2 = tJDistributed(nsites, nup, ndn);

          double e0 = eigval0(ops, block);
          double e02 = eigval0(ops, block2);
          REQUIRE(isapprox(e0, e02));
        }
      }
    }
  }

  Log.out("tj_distributed: HB triangular N=12 complex exchange");
  int nsites = 12;
  std::vector<double> etas = {0.00, 0.01, 0.02,
                              0.03, 0.04, 0.05}; // dont change etas :-)
  for (auto eta : etas) {
    Log("eta: {}", eta);
    for (int nup = 0; nup <= nsites; ++nup) {
      auto [ops, e00] = testcases::spinhalf::triangular_12_complex(nup, eta);

      auto block = Spinhalf(nsites, nup);
      auto block_tJ = tJDistributed(nsites, nup, nsites - nup);
      double e0_spinhalf = eigval0(ops, block);
      double e0 = eigval0(ops, block_tJ);
      // Log("{} {} {} {}", nup, e0_spinhalf, e0, e00);
      REQUIRE(isapprox(e0_spinhalf, e0));
      if (nup == 6) {
        REQUIRE(std::abs(e0 - e00) < 1e-8);
      }
    }
  }

  test_onsite("Sz", "SzSz");
  test_onsite("Ntot", "NtotNtot");

  for (int nsites = 2; nsites < 6; ++nsites) {
    for (int nup = 1; nup < nsites; ++nup) {
      for (int ndn = 1; ndn < nsites - nup; ++ndn) {
        auto b = tJDistributed(nsites, nup, ndn);
        for (int s = 0; s < nsites; ++s) {

          auto r = random_state(b);

          // Exchange
          auto cdagup = Op("Cdagup", s);
          auto cup = Op("Cup", s);
          auto cdagdn = Op("Cdagdn", s);
          auto cdn = Op("Cdn", s);

          auto spsmr = apply(cdagup, apply(cdn, apply(cdagdn, apply(cup,
          r)))); auto smspr = apply(cdagdn, apply(cup, apply(cdagup,
          apply(cdn, r)))); auto r1 = 0.5 * (spsmr + smspr); auto r2 =
          apply(Op("Exchange", {s, s}), r); REQUIRE(isapprox(r1, r2));

          // SdotS
          auto szszr = apply(Op("SzSz", {s, s}), r);
          r1 = 0.5 * (spsmr + smspr) + szszr;
          r2 = apply(Op("SdotS", {s, s}), r);
          REQUIRE(isapprox(r1, r2));

          // tJSzSz
          r1 = szszr - 0.25 * apply(Op("NtotNtot", {s, s}), r);
          r2 = apply(Op("tJSzSz", {s, s}), r);
          REQUIRE(isapprox(r1, r2));

          // tJSdotS
          r1 = 0.5 * (spsmr + smspr) + szszr -
               0.25 * apply(Op("NtotNtot", {s, s}), r);
          r2 = apply(Op("tJSdotS", {s, s}), r);
          REQUIRE(isapprox(r1, r2));

          // Hopup
          r1 = -(apply(cdagup, apply(cup, r)) + apply(cdagup, apply(cup,
          r))); r2 = apply(Op("Hopup", {s, s}), r); REQUIRE(isapprox(r1,
          r2));

          // Hopdn
          r1 = -(apply(cdagdn, apply(cdn, r)) + apply(cdagdn, apply(cdn,
          r))); r2 = apply(Op("Hopdn", {s, s}), r); REQUIRE(isapprox(r1,
          r2));

          // Hop
          r1 = -(apply(cdagup, apply(cup, r)) + apply(cdagup, apply(cup, r))
          +
                 apply(cdagdn, apply(cdn, r)) + apply(cdagdn, apply(cdn,
                 r)));
          r2 = apply(Op("Hop", {s, s}), r);
          REQUIRE(isapprox(r1, r2));
        }
      }
    }
  }

  // Test corrs
  for (int nsites = 2; nsites < 5; ++nsites) {
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites - nup; ++ndn) {
        auto b = tJDistributed(nsites, nup, ndn);
        auto r = random_state(b);

        for (int i = 0; i < nsites; ++i) {
          auto a = apply(Op("Nup", i), r) + apply(Op("Ndn", i), r);
          auto b = apply(Op("Ntot", i), r);
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
          }
        }
      }
    }
  }

} catch (Error e) {
  error_trace(e);
}
