#include "../../catch.hpp"

#include <iostream>

#include "testcases_tj.h"
#include <hydra/all.h>

using namespace hydra;

void test_tjmodel_e0_real(BondList bonds, Couplings couplings, int nup, int ndn,
                          double e0) {
  int n_sites = bonds.n_sites();
  // if ((nup == 1) && (ndn == 2)) {
  auto block = tJ<uint32_t>(n_sites, nup, ndn);
  auto H = MatrixReal(bonds, couplings, block, block);

  auto eigs = lila::EigenvaluesSym(H);
  // Log("nup: {}, ndn: {}, e0: {}, eigs(0): {}", nup, ndn, e0, eigs(0));
  REQUIRE(lila::close(H, lila::Herm(H)));
  REQUIRE(std::abs(e0 - eigs(0)) < 1e-6);
  // }
}

void test_tjmodel_fulleigs(BondList bonds, Couplings couplings,
                           lila::Vector<double> exact_eigs) {
  int n_sites = bonds.n_sites();

  lila::Vector<double> all_eigs;
  for (int ndn = 0; ndn <= n_sites; ++ndn) {
    for (int nup = 0; nup <= n_sites - ndn; ++nup) {

      auto block = tJ<uint32_t>(n_sites, nup, ndn);
      auto H = MatrixReal(bonds, couplings, block, block);
      REQUIRE(lila::close(H, lila::Herm(H)));
      auto eigs = lila::EigenvaluesSym(H);
      // Log("nup {} ndn {}", nup, ndn);
      // LilaPrint(eigs);

      for (auto eig : eigs)
        all_eigs.push_back(eig);
    }
  }
  std::sort(all_eigs.begin(), all_eigs.end());
  REQUIRE(all_eigs.size() == exact_eigs.size());
  // LilaPrint(all_eigs);
  // LilaPrint(exact_eigs);
  REQUIRE(lila::close(all_eigs, exact_eigs));
}

TEST_CASE("tj_matrix", "[blocks][tj]") {
  using namespace hydra::testcases::tj;

  {
    Log("tj_matrix: TJModel: six-site chain test, t=1.0, J=1.0, N=6");
    auto [bonds, cpls] = tJchain(6, 1.0, 1.0);
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
      test_tjmodel_e0_real(bonds, cpls, nup, ndn, e0);
  }

  {
    Log.out("tj_matrix: TJModel: six-site chain test, t=1.0, J=0.0, N=6");
    auto [bonds, cpls] = tJchain(6, 1.0, 0.0);
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
      test_tjmodel_e0_real(bonds, cpls, nup, ndn, e0);
  }

  for (int L = 3; L <= 6; ++L) {
    Log.out("tj_matrix: ALPS full spectrum test, chain N={}", L);
    auto [bonds, cpls, eigs] = tJchain_fullspectrum_alps(L);
    test_tjmodel_fulleigs(bonds, cpls, eigs);
  }

  {
    Log.out("tj_matrix: ALPS full spectrum test, square 2x2");
    auto [bonds, cpls, eigs] = tj_square2x2_fullspectrum_alps();
    test_tjmodel_fulleigs(bonds, cpls, eigs);
  }

  for (int N = 3; N <= 6; ++N) {
    Log.out("tj_matrix:  random all-to-all complex exchange test, N={}", N);

    auto [bonds, cpls] = tj_alltoall_complex(N);
    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N - nup; ++ndn) {
        auto block = tJ<uint32_t>(N, nup, ndn);
        auto H = MatrixCplx(bonds, cpls, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));
      }

    // Check whether eigenvalues agree with HB model
    for (int nup = 0; nup <= N; ++nup) {
      int ndn = N - nup;
      auto block1 = tJ<uint32_t>(N, nup, ndn);
      auto block2 = Spinhalf<uint32_t>(N, nup);
      auto H1 = MatrixCplx(bonds, cpls, block1, block1);
      auto H2 = MatrixCplx(bonds, cpls, block2, block2);
      auto eigs1 = lila::EigenvaluesSym(H1);
      auto eigs2 = lila::EigenvaluesSym(H2);
      // Log("eigs1(0): {}, eigs2(0): {}", eigs1(0), eigs2(0));
      REQUIRE(lila::close(eigs1, eigs2));
    }
  }

  {
    Log.out("tj_matrix: Henry's Matlab test, random 3");
    auto [bonds, cpls, eigs] = randomAlltoAll3();
    test_tjmodel_fulleigs(bonds, cpls, eigs);
  }

  {
    Log.out("tj_matrix: Henry's Matlab test, random 4");
    auto [bonds, cpls, eigs] = randomAlltoAll4();
    test_tjmodel_fulleigs(bonds, cpls, eigs);
  }
}
