// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <iostream>

#include "../spinhalf/testcases_spinhalf.hpp"
#include "testcases_tj.hpp"
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/states/create_state.hpp>

using namespace xdiag;

static void test_onsite(std::string op1, std::string op12) {
  for (int nsites = 2; nsites < 5; ++nsites) {
    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites - nup; ++ndn) {
        auto b = tJ(nsites, nup, ndn);
        for (int s = 0; s < nsites; ++s) {
          arma::mat m1 = matrix(Op(op1, s), b);
          arma::mat m12 = matrix(Op(op12, {s, s}), b);
          REQUIRE(isapprox(m12, m1 * m1));
        }
      }
    }
  }
}

static void test_tjmodel_e0_real(OpSum ops, int nsites, int nup, int ndn,
                                 double e0) {
  auto block = tJ(nsites, nup, ndn);
  auto H = matrix(ops, block, block);
  arma::vec eigs;
  arma::eig_sym(eigs, H);

  REQUIRE(arma::norm(H - H.t()) < 1e-12);
  REQUIRE(std::abs(e0 - eigs(0)) < 1e-6);
}

static void test_tjmodel_fulleigs(OpSum ops, int nsites,
                                  arma::Col<double> exact_eigs) {

  std::vector<double> all_eigs;
  for (int ndn = 0; ndn <= nsites; ++ndn) {
    for (int nup = 0; nup <= nsites - ndn; ++nup) {

      auto block = tJ(nsites, nup, ndn);
      auto H = matrix(ops, block, block);
      // H.print();
      REQUIRE(arma::norm(H - H.t()) < 1e-12);

      arma::vec eigs;
      arma::eig_sym(eigs, H);

      for (auto eig : eigs)
        all_eigs.push_back(eig);
    }
  }
  std::sort(all_eigs.begin(), all_eigs.end());
  REQUIRE(all_eigs.size() == exact_eigs.size());
  // XDIAG_SHOW(all_eigs);
  // XDIAG_SHOW(exact_eigs);
  REQUIRE(isapprox(arma::vec(all_eigs), exact_eigs));
}

TEST_CASE("tj_matrix", "[tj]") try {
  using namespace xdiag::testcases::tj;

  {
    Log("tj_matrix: HB all-to-all comparison");
    for (int nsites = 2; nsites < 7; ++nsites) {
      Log("N: {}", nsites);
      int nup = nsites / 2;
      auto ops = testcases::spinhalf::HB_alltoall(nsites);
      auto block = Spinhalf(nsites, nup);
      auto block_tJ = tJ(nsites, nup, nsites - nup);
      auto H = matrix(ops, block, block);
      auto H_tJ = matrix(ops, block_tJ, block_tJ);

      REQUIRE(H.is_hermitian());
      REQUIRE(H_tJ.is_hermitian());

      // H.print();
      // H_tJ.print();

      arma::vec eigs;
      arma::eig_sym(eigs, H);

      arma::vec eigs_tJ;
      arma::eig_sym(eigs_tJ, H_tJ);

      // eigs.print();
      // eigs_tJ.print();

      REQUIRE(isapprox(eigs, eigs_tJ));
    }
  }

  {
    Log("tj_matrix: TJModel: six-site chain test, t=1.0, J=1.0, N=6");
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
    for (auto [nup, ndn, e0] : nup_ndn_e0)
      test_tjmodel_e0_real(ops, 6, nup, ndn, e0);
  }

  {
    Log.out("tj_matrix: TJModel: six-site chain test, t=1.0, J=0.0, N=6");
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
    for (auto [nup, ndn, e0] : nup_ndn_e0)
      test_tjmodel_e0_real(ops, 6, nup, ndn, e0);
  }

  for (int L = 3; L <= 6; ++L) {
    Log.out("tj_matrix: ALPS full spectrum test, chain N={}", L);
    auto [ops, eigs] = tJchain_fullspectrum_alps(L);
    test_tjmodel_fulleigs(ops, L, eigs);
  }

  {
    Log.out("tj_matrix: ALPS full spectrum test, square 2x2");
    auto [ops, eigs] = tj_square2x2_fullspectrum_alps();
    test_tjmodel_fulleigs(ops, 4, eigs);
  }

  for (int N = 3; N <= 6; ++N) {
    Log.out("tj_matrix:  random all-to-all complex exchange test, N={}", N);

    auto ops = tj_alltoall_complex(N);
    OpSum ops_hb;
    for (auto [cpl, op] : ops) {
      if (op.type() == "HB") {
        ops_hb += cpl * op;
        std::string name = cpl.string();
        ops_hb[name] = ops[name];
      }
    }

    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N - nup; ++ndn) {
        auto block = tJ(N, nup, ndn);
        auto H = matrixC(ops, block, block);
        REQUIRE(arma::norm(H - H.t()) < 1e-12);
      }

    // Check whether eigenvalues agree with HB model
    for (int nup = 0; nup <= N; ++nup) {
      int ndn = N - nup;
      auto block1 = tJ(N, nup, ndn);
      auto block2 = Spinhalf(N, nup);

      auto H1 = matrixC(ops_hb, block1, block1);
      auto H2 = matrixC(ops_hb, block2, block2);
      arma::vec eigs1;
      arma::eig_sym(eigs1, H1);
      arma::vec eigs2;
      arma::eig_sym(eigs2, H2);
      // Log("eigs1(0): {}, eigs2(0): {}", eigs1(0), eigs2(0));
      REQUIRE(isapprox(eigs1, eigs2));
    }
  }

  {
    Log.out("tj_matrix: Henry's Matlab test, random 3");
    auto [ops, eigs] = randomAlltoAll3();
    test_tjmodel_fulleigs(ops, 3, eigs);
  }

  {
    Log.out("tj_matrix: Henry's Matlab test, random 4");
    auto [ops, eigs] = randomAlltoAll4();
    test_tjmodel_fulleigs(ops, 4, eigs);
  }

  test_onsite("Sz", "SzSz");
  test_onsite("Ntot", "NtotNtot");

  for (int nsites = 2; nsites < 6; ++nsites) {
    for (int nup = 1; nup < nsites; ++nup) {
      for (int ndn = 1; ndn < nsites - nup; ++ndn) {
        auto b = tJ(nsites, nup, ndn);
        for (int s = 0; s < nsites; ++s) {

          auto r = random_state(b);

          // Exchange
          auto cdagup = Op("Cdagup", s);
          auto cup = Op("Cup", s);
          auto cdagdn = Op("Cdagdn", s);
          auto cdn = Op("Cdn", s);

          auto spsmr = apply(cdagup, apply(cdn, apply(cdagdn, apply(cup, r))));
          auto smspr = apply(cdagdn, apply(cup, apply(cdagup, apply(cdn, r))));
          auto r1 = 0.5 * (spsmr + smspr);
          auto r2 = apply(Op("Exchange", {s, s}), r);
          REQUIRE(isapprox(r1, r2));

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
          r1 = -(apply(cdagup, apply(cup, r)) + apply(cdagup, apply(cup, r)));
          r2 = apply(Op("Hopup", {s, s}), r);
          REQUIRE(isapprox(r1, r2));

          // Hopdn
          r1 = -(apply(cdagdn, apply(cdn, r)) + apply(cdagdn, apply(cdn, r)));
          r2 = apply(Op("Hopdn", {s, s}), r);
          REQUIRE(isapprox(r1, r2));

          // Hop
          r1 = -(apply(cdagup, apply(cup, r)) + apply(cdagup, apply(cup, r)) +
                 apply(cdagdn, apply(cdn, r)) + apply(cdagdn, apply(cdn, r)));
          r2 = apply(Op("Hop", {s, s}), r);
          REQUIRE(isapprox(r1, r2));
        }
      }
    }
  }

} catch (Error e) {
  xdiag::error_trace(e);
}
