#include "../../catch.hpp"

#include <iostream>

#include "../spinhalf/testcases_spinhalf.hpp"
#include "testcases_tj.hpp"
#include <xdiag/blocks/spinhalf/spinhalf.hpp>
#include <xdiag/blocks/spinhalf/matrix.hpp>
#include <xdiag/blocks/tj/matrix.hpp>
#include <xdiag/utils/close.hpp>
#include <xdiag/utils/print_macro.hpp>

using namespace xdiag;

void test_tjmodel_e0_real(OpSum ops, int n_sites, int nup, int ndn,
                          double e0) {
  auto block = tJ(n_sites, nup, ndn);
  auto H = matrix(ops, block, block);
  arma::vec eigs;
  arma::eig_sym(eigs, H);

  REQUIRE(arma::norm(H - H.t()) < 1e-12);
  REQUIRE(std::abs(e0 - eigs(0)) < 1e-6);
}

void test_tjmodel_fulleigs(OpSum ops, int n_sites,
                           arma::Col<double> exact_eigs) {

  std::vector<double> all_eigs;
  for (int ndn = 0; ndn <= n_sites; ++ndn) {
    for (int nup = 0; nup <= n_sites - ndn; ++nup) {

      auto block = tJ(n_sites, nup, ndn);
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
  REQUIRE(close(arma::vec(all_eigs), exact_eigs));
}

TEST_CASE("tj_matrix", "[tj]") try {
  using namespace xdiag::testcases::tj;

  {
    Log("tj_matrix: HB all-to-all comparison");
    for (int n_sites = 2; n_sites < 7; ++n_sites) {
      Log("N: {}", n_sites);
      int nup = n_sites / 2;
      auto ops = testcases::spinhalf::HB_alltoall(n_sites);
      auto block = Spinhalf(n_sites, nup);
      auto block_tJ = tJ(n_sites, nup, n_sites - nup);
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

      REQUIRE(close(eigs, eigs_tJ));
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
    for (auto op : ops) {
      if (op.type() == "HB") {
        ops_hb += op;
        std::string name = op.coupling().as<std::string>();
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
      REQUIRE(close(eigs1, eigs2));
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
} catch (Error e) {
  xdiag::error_trace(e);
}
