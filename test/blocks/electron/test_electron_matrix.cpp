#include "../../catch.hpp"

#include <iostream>

#include "../spinhalf/testcases_spinhalf.h"
#include "../tj/testcases_tj.h"
#include "testcases_electron.h"

#include <hydra/all.h>

using namespace hydra;

void test_electron_np_no_np_matrix(int n_sites, BondList bonds, Couplings cpls) {

  auto block_full = Electron(n_sites);
  auto H_full = MatrixCplx(bonds, cpls, block_full, block_full);
  REQUIRE(lila::close(H_full, lila::Herm(H_full)));
  auto eigs_full = lila::EigenvaluesSym(H_full);

  lila::Vector<double> all_eigs;
  for (int nup = 0; nup <= n_sites; ++nup)
    for (int ndn = 0; ndn <= n_sites; ++ndn) {

      auto block = Electron(n_sites, nup, ndn);
      auto H = MatrixCplx(bonds, cpls, block, block);
      REQUIRE(lila::close(H, lila::Herm(H)));
      auto eigs = lila::EigenvaluesSym(H);
      // LilaPrint(eigs);

      for (auto eig : eigs)
        all_eigs.push_back(eig);
    }
  std::sort(all_eigs.begin(), all_eigs.end());
  REQUIRE(lila::close(all_eigs, eigs_full));
}

TEST_CASE("electron_matrix", "[blocks][electron]") {
  using namespace hydra::testcases::electron;

  BondList bondlist;
  Couplings couplings;

  // Compare with matrix from Weisse and Fehske
  lila::Log("electron_matrix: Hubbard Weisse & Fehske");
  int n_sites = 4;
  int n_up = 3;
  int n_dn = 2;
  double t = 1.0;
  double U = 5.0;
  auto block = Electron<uint32_t>(n_sites, n_up, n_dn);

  for (int i = 0; i < n_sites; ++i)
    bondlist << Bond("HOP", "T", {i, (i + 1) % n_sites});
  couplings["T"] = 1.0;
  couplings["U"] = 5.0;

  auto H1 = MatrixReal(bondlist, couplings, block, block);

  double tp = t;
  double tm = -t;
  double UU = U;
  double U2 = 2 * U;
  lila::Matrix<double> H1_correct = {
      {U2, tm, 0., 0., tp, 0., tm, 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., tm, 0., 0., 0., 0., 0.}, // 1
      {tm, U2, tm, tm, 0., tp, 0., tm, 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., tm, 0., 0., 0., 0.}, // 2
      {0., tm, U2, 0., tm, 0., 0., 0., tm, 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., tm, 0., 0., 0.}, // 3
      {0., tm, 0., UU, tm, 0., 0., 0., 0., tm, 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., tm, 0., 0.}, // 4
      {tp, 0., tm, tm, UU, tm, 0., 0., 0., 0., tm, 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., tm, 0.}, // 5
      {0., tp, 0., 0., tm, UU, 0., 0., 0., 0., 0., tm,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., tm}, // 6
      {tm, 0., 0., 0., 0., 0., U2, tm, 0., 0., tp, 0.,
       tm, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, // 7
      {0., tm, 0., 0., 0., 0., tm, UU, tm, tm, 0., tp,
       0., tm, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, // 8
      {0., 0., tm, 0., 0., 0., 0., tm, UU, 0., tm, 0.,
       0., 0., tm, 0., 0., 0., 0., 0., 0., 0., 0., 0.}, // 9
      {0., 0., 0., tm, 0., 0., 0., tm, 0., U2, tm, 0.,
       0., 0., 0., tm, 0., 0., 0., 0., 0., 0., 0., 0.}, // 10
      {0., 0., 0., 0., tm, 0., tp, 0., tm, tm, U2, tm,
       0., 0., 0., 0., tm, 0., 0., 0., 0., 0., 0., 0.}, // 11
      {0., 0., 0., 0., 0., tm, 0., tp, 0., 0., tm, UU,
       0., 0., 0., 0., 0., tm, 0., 0., 0., 0., 0., 0.}, // 12

      {0., 0., 0., 0., 0., 0., tm, 0., 0., 0., 0., 0.,
       UU, tm, 0., 0., tp, 0., tm, 0., 0., 0., 0., 0.}, // 13
      {0., 0., 0., 0., 0., 0., 0., tm, 0., 0., 0., 0.,
       tm, U2, tm, tm, 0., tp, 0., tm, 0., 0., 0., 0.}, // 14
      {0., 0., 0., 0., 0., 0., 0., 0., tm, 0., 0., 0.,
       0., tm, UU, 0., tm, 0., 0., 0., tm, 0., 0., 0.}, // 15
      {0., 0., 0., 0., 0., 0., 0., 0., 0., tm, 0., 0.,
       0., tm, 0., U2, tm, 0., 0., 0., 0., tm, 0., 0.}, // 16
      {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., tm, 0.,
       tp, 0., tm, tm, UU, tm, 0., 0., 0., 0., tm, 0.}, // 17
      {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., tm,
       0., tp, 0., 0., tm, U2, 0., 0., 0., 0., 0., tm}, // 18
      {tm, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       tm, 0., 0., 0., 0., 0., UU, tm, 0., 0., tp, 0.}, // 19
      {0., tm, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., tm, 0., 0., 0., 0., tm, UU, tm, tm, 0., tp}, // 20
      {0., 0., tm, 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., tm, 0., 0., 0., 0., tm, U2, 0., tm, 0.}, // 21
      {0., 0., 0., tm, 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., tm, 0., 0., 0., tm, 0., UU, tm, 0.}, // 22
      {0., 0., 0., 0., tm, 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., tm, 0., tp, 0., tm, tm, U2, tm}, // 23
      {0., 0., 0., 0., 0., tm, 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., tm, 0., tp, 0., 0., tm, U2} // 24
  };

  // LilaPrint(H1);
  // LilaPrint(H1_correct);

  REQUIRE(lila::close(lila::Real(H1), H1_correct));
  for (int i = 0; i < 24; ++i)
    for (int j = 0; j < 24; ++j) {
      REQUIRE(lila::close(std::real(H1(i, j)), H1_correct(i, j)));
    }

  //////////////////////////////////
  // Test two site exact solution
  bondlist.clear();
  couplings.clear();
  bondlist << Bond("HOP", "T", {0, 1});
  auto block2 = Electron<uint32_t>(2, 1, 1);
  for (int i = 0; i < 20; ++i) {
    double U = 1.234 * i;
    lila::Log("electron_matrix: two-site exact solution test, U={}", U);
    couplings["T"] = 1.0;
    couplings["U"] = U;
    double e0_exact = 0.5 * (U - sqrt(U * U + 16));
    auto H = MatrixReal(bondlist, couplings, block2, block2);
    REQUIRE(lila::close(H, lila::Herm(H)));
    auto evecs = lila::EigenvaluesSym(H);
    double e0 = evecs(0);
    // printf("e0: %f, e0_exact: %f\n", e0, e0_exact);
    REQUIRE(lila::close(e0_exact, e0));
  }

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (real)
  for (int n_sites = 3; n_sites < 7; ++n_sites) {
    lila::Log(
        "electron_matrix: free fermion random all-to-all test, (real), N={}",
        n_sites);

    std::tie(bondlist, couplings) = freefermion_alltoall(n_sites);

    // Create single particle matrix
    auto Hs = lila::Zeros<double>(n_sites, n_sites);
    for (auto bond : bondlist) {
      assert(bond.size() == 2);
      int s1 = bond.site(0);
      int s2 = bond.site(1);
      auto name = bond.coupling();
      Hs(s1, s2) = -lila::real(couplings[name]);
      Hs(s2, s1) = -lila::real(couplings[name]);
    }
    auto seigs = lila::EigenvaluesSym(Hs);
    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {

        // Compute exact gs energy
        double e0_exact = 0;
        for (int i = 0; i < nup; ++i)
          e0_exact += seigs(i);
        for (int i = 0; i < ndn; ++i)
          e0_exact += seigs(i);

        auto block3 = Electron<uint32_t>(n_sites, nup, ndn);
        auto Hr = MatrixReal(bondlist, couplings, block3, block3);
        REQUIRE(lila::close(Hr, lila::Herm(Hr)));
        auto evecsr = lila::EigenvaluesSym(Hr);
        auto Hc = MatrixCplx(bondlist, couplings, block3, block3);
        REQUIRE(lila::close(Hc, lila::Herm(Hc)));
        auto evecsc = lila::EigenvaluesSym(Hc);
        REQUIRE(lila::close(evecsr, evecsc));
        double e0 = evecsr(0);
        // printf("nup: %d, ndn: %d, dim: %d, e0: %f, e0_exact: %f\n", nup, ndn,
        //        (int)evecs.size(), e0, e0_exact);
        REQUIRE(lila::close(e0_exact, e0));
      }
  }

  //////////////////////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (cplx, up/dn different)
  for (int n_sites = 3; n_sites < 7; ++n_sites) {
    lila::Log(
        "electron_matrix: free fermion random all-to-all test, (cplx), N={}",
        n_sites);
    std::tie(bondlist, couplings) = freefermion_alltoall_complex_updn(n_sites);

    // Create single particle matrix for upspins
    auto Hs_up = lila::Zeros<complex>(n_sites, n_sites);
    for (auto bond : bondlist.bonds_of_type("HOPUP")) {
      assert(bond.size() == 2);
      int s1 = bond.site(0);
      int s2 = bond.site(1);
      auto name = bond.coupling();
      Hs_up(s1, s2) = -couplings[name];
      Hs_up(s2, s1) = -couplings[name];
    }
    auto seigs_up = lila::EigenvaluesSym(Hs_up);

    // Create single particle matrix for dnspins
    auto Hs_dn = lila::Zeros<complex>(n_sites, n_sites);
    for (auto bond : bondlist.bonds_of_type("HOPDN")) {
      assert(bond.size() == 2);
      int s1 = bond.site(0);
      int s2 = bond.site(1);
      auto name = bond.coupling();
      Hs_dn(s1, s2) = -couplings[name];
      Hs_dn(s2, s1) = -lila::conj(couplings[name]);
    }
    auto seigs_dn = lila::EigenvaluesSym(Hs_dn);

    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {

        // Compute exact gs energy
        double e0_exact = 0;
        for (int i = 0; i < nup; ++i)
          e0_exact += seigs_up(i);
        for (int i = 0; i < ndn; ++i)
          e0_exact += seigs_dn(i);

        auto block3 = Electron<uint32_t>(n_sites, nup, ndn);
        auto H = MatrixCplx(bondlist, couplings, block3, block3);
        REQUIRE(lila::close(H, lila::Herm(H)));
        auto evecs = lila::EigenvaluesSym(H);
        double e0 = evecs(0);
        // printf("nup: %d, ndn: %d, dim: %d, e0: %f, e0_exact: %f\n", nup, ndn,
        //        (int)evecs.size(), e0, e0_exact);
        REQUIRE(lila::close(e0_exact, e0));
      }
  }

  // Test Heisenberg terms at half-filling
  for (int n_sites = 2; n_sites <= 6; ++n_sites) {
    lila::Log("electron_matrix: Heisenberg all-to-all comparison test, N={}",
              n_sites);
    int nup = n_sites / 2;
    int ndn = n_sites - nup;
    auto block_spinhalf = Spinhalf(n_sites, nup);
    auto block_electron = Electron(n_sites, nup, ndn);

    auto [bonds, couplings] = testcases::spinhalf::HB_alltoall(n_sites);
    couplings["U"] = 99999; // gap out doubly occupied sites
    auto H_spinhalf =
        MatrixReal(bonds, couplings, block_spinhalf, block_spinhalf);
    auto H_electron =
        MatrixReal(bonds, couplings, block_electron, block_electron);
    REQUIRE(lila::close(H_spinhalf, lila::Herm(H_spinhalf)));
    REQUIRE(lila::close(H_electron, lila::Herm(H_electron)));

    auto eigs_spinhalf = lila::EigenvaluesSym(H_spinhalf);
    auto eigs_electron = lila::EigenvaluesSym(H_electron);
    for (int idx = 0; idx < eigs_spinhalf.size(); ++idx)
      REQUIRE(std::abs(eigs_spinhalf(idx) - eigs_electron(idx)) < 1e-5);
  }

  /////////////////////////////////////////////////////
  // Henry's MATLAB code test (tests Heisenberg terms)
  lila::Log(
      "electron_matrix: U-hopping-HB full spectrum of Henry's Matlab code");

  lila::Vector<double> eigs_correct;
  std::tie(bondlist, couplings, eigs_correct) = randomAlltoAll4NoU();

  // Compute full spectrum in hydra
  n_sites = 4;
  lila::Vector<double> all_eigs;
  for (int nup = 0; nup <= n_sites; ++nup)
    for (int ndn = 0; ndn <= n_sites; ++ndn) {
      auto block = Electron(n_sites, nup, ndn);
      auto H = MatrixReal(bondlist, couplings, block, block);
      REQUIRE(lila::close(H, lila::Herm(H)));
      auto eigs = lila::EigenvaluesSym(H);
      for (auto eig : eigs)
        all_eigs.push_back(eig);
    }
  std::sort(all_eigs.begin(), all_eigs.end());
  REQUIRE(lila::close(all_eigs, eigs_correct));

  std::tie(bondlist, couplings, eigs_correct) = randomAlltoAll4();
  all_eigs.clear();
  for (int nup = 0; nup <= n_sites; ++nup)
    for (int ndn = 0; ndn <= n_sites; ++ndn) {

      auto block = Electron(n_sites, nup, ndn);
      auto H = MatrixReal(bondlist, couplings, block, block);
      REQUIRE(lila::close(H, lila::Herm(H)));
      auto eigs = lila::EigenvaluesSym(H);
      for (auto eig : eigs)
        all_eigs.push_back(eig);
    }
  std::sort(all_eigs.begin(), all_eigs.end());
  REQUIRE(lila::close(all_eigs, eigs_correct));

  for (int N = 3; N <= 6; ++N) {
    lila::Log.out(
        "electron_matrix: random all-to-all complex exchange test, N={}", N);

    auto [bonds, cpls] = hydra::testcases::tj::tj_alltoall_complex(N);
    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N - nup; ++ndn) {
        auto block = Electron<uint32_t>(N, nup, ndn);
        auto H = MatrixCplx(bonds, cpls, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));
      }

    // Set hoppings to zero
    for (auto cpl : cpls) {
      if (cpl.first[0] == 'T') {
        cpls[cpl.first] = 0.0;
        // std::cout << cpl.first << " " << cpl.first[0] << "\n";
      }
    }
    cpls["U"] = 1000;

    // Check whether eigenvalues agree with HB model
    for (int nup = 0; nup <= N; ++nup) {
      int ndn = N - nup;
      auto block1 = Electron<uint32_t>(N, nup, ndn);
      auto block2 = Spinhalf<uint32_t>(N, nup);
      auto H1 = MatrixCplx(bonds, cpls, block1, block1);
      auto H2 = MatrixCplx(bonds, cpls, block2, block2);
      auto eigs1 = lila::EigenvaluesSym(H1);
      auto eigs2 = lila::EigenvaluesSym(H2);
      // lila::Log("eigs1(0): {}, eigs2(0): {}", eigs1(0), eigs2(0));
      lila::Vector<double> eigs1_sub = eigs1({0, eigs2.size()});
      // LilaPrint(eigs1_sub);
      // LilaPrint(eigs2);
      REQUIRE(lila::close(eigs1_sub, eigs2, 1e-6, 1e-6));
    }
  }

  for (int N = 3; N <= 4; ++N) {
    lila::Log.out("electron_matrix: random all-to-all complex exchange test Np "
                  "<-> NoNp, N={}",
                  N);
    auto [bonds, cpls] = hydra::testcases::tj::tj_alltoall_complex(N);
    test_electron_np_no_np_matrix(N, bonds, cpls);
  }
}
