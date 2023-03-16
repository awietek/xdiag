#include "../../catch.hpp"

#include <iostream>

#include "../spinhalf/testcases_spinhalf.h"
#include "../tj/testcases_tj.h"
#include "testcases_electron.h"

#include <hydra/all.h>

using namespace hydra;

void test_electron_np_no_np_matrix(int n_sites, BondList bonds) {

  auto block_full = Electron(n_sites);
  auto H_full = matrix_cplx(bonds, block_full, block_full);
  REQUIRE(H_full.is_hermitian(1e-12));
  arma::Col<double> eigs_full;
  arma::eig_sym(eigs_full, H_full);

  std::vector<double> all_eigs;
  for (int nup = 0; nup <= n_sites; ++nup)
    for (int ndn = 0; ndn <= n_sites; ++ndn) {

      auto block = Electron(n_sites, nup, ndn);
      auto H = matrix_cplx(bonds, block, block);
      REQUIRE(H.is_hermitian(1e-12));

      arma::Col<double> eigs;
      arma::eig_sym(eigs, H);

      for (auto eig : eigs)
        all_eigs.push_back(eig);
    }
  std::sort(all_eigs.begin(), all_eigs.end());
  REQUIRE(close(arma::Col(all_eigs), eigs_full));
}

TEST_CASE("electron_matrix", "[electron]") {
  using namespace hydra::testcases::electron;

  BondList bondlist;

  // Compare with matrix from Weisse and Fehske
  Log("electron_matrix: Hubbard Weisse & Fehske");
  int n_sites = 4;
  int n_up = 3;
  int n_dn = 2;
  double t = 1.0;
  double U = 5.0;
  auto block = Electron(n_sites, n_up, n_dn);

  for (int i = 0; i < n_sites; ++i)
    bondlist << Bond("HOP", "T", {i, (i + 1) % n_sites});
  bondlist["T"] = 1.0;
  bondlist["U"] = 5.0;
  auto H1 = matrix_real(bondlist, block, block);
  double tp = t;
  double tm = -t;
  double UU = U;
  double U2 = 2 * U;
  arma::Mat<double> H1_correct = {
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

  REQUIRE(close(H1, H1_correct));
  for (int i = 0; i < 24; ++i)
    for (int j = 0; j < 24; ++j) {
      REQUIRE(close(std::real(H1(i, j)), H1_correct(i, j)));
    }

  //////////////////////////////////
  // Test two site exact solution
  bondlist.clear();
  bondlist << Bond("HOP", "T", {0, 1});
  auto block2 = Electron(2, 1, 1);
  for (int i = 0; i < 20; ++i) {
    double U = 1.234 * i;
    Log("electron_matrix: two-site exact solution test, U={}", U);
    bondlist["T"] = 1.0;
    bondlist["U"] = U;
    double e0_exact = 0.5 * (U - sqrt(U * U + 16));
    auto H = matrix_real(bondlist, block2, block2);
    REQUIRE(H.is_hermitian(1e-8));
    arma::Col<double> eigs;
    arma::eig_sym(eigs, H);
    double e0 = eigs(0);
    // printf("e0: %f, e0_exact: %f\n", e0, e0_exact);
    REQUIRE(close(e0_exact, e0));
  }

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (real)
  for (int n_sites = 3; n_sites < 7; ++n_sites) {
    Log("electron_matrix: free fermion random all-to-all test, (real), N={}",
        n_sites);

    bondlist = freefermion_alltoall(n_sites);

    // Create single particle matrix
    arma::Mat<double> Hs(n_sites, n_sites, arma::fill::zeros);
    for (auto bond : bondlist) {
      assert(bond.size() == 2);
      int s1 = bond.site(0);
      int s2 = bond.site(1);
      auto name = bond.coupling_name();
      Hs(s1, s2) = -bondlist[name].as<double>();
      Hs(s2, s1) = -bondlist[name].as<double>();
    }

    arma::vec seigs;
    arma::eig_sym(seigs, Hs);
    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {

        // Compute exact gs energy
        double e0_exact = 0;
        for (int i = 0; i < nup; ++i)
          e0_exact += seigs(i);
        for (int i = 0; i < ndn; ++i)
          e0_exact += seigs(i);

        auto block3 = Electron(n_sites, nup, ndn);
        auto Hr = matrix_real(bondlist, block3, block3);
        REQUIRE(Hr.is_hermitian(1e-8));
        arma::vec eigsr;
        arma::eig_sym(eigsr, Hr);
        auto Hc = matrix_cplx(bondlist, block3, block3);
        REQUIRE(Hc.is_hermitian(1e-8));
        arma::vec eigsc;
        arma::eig_sym(eigsc, Hc);
        REQUIRE(close(eigsr, eigsc));
        double e0 = eigsr(0);
        // printf("nup: %d, ndn: %d, dim: %d, e0: %f, e0_exact: %f\n", nup, ndn,
        //        (int)evecs.size(), e0, e0_exact);
        REQUIRE(close(e0_exact, e0));
      }
  }

  //////////////////////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (cplx, up/dn different)
  for (int n_sites = 3; n_sites < 7; ++n_sites) {
    Log("electron_matrix: free fermion random all-to-all test, (cplx), N={}",
        n_sites);
    bondlist = freefermion_alltoall_complex_updn(n_sites);

    // Create single particle matrix for upspins
    arma::cx_mat Hs_up(n_sites, n_sites, arma::fill::zeros);
    for (auto bond : bondlist.bonds_of_type("HOPUP")) {
      assert(bond.size() == 2);
      int s1 = bond.site(0);
      int s2 = bond.site(1);
      auto name = bond.coupling_name();
      Hs_up(s1, s2) = -bondlist[name].as<complex>();
      Hs_up(s2, s1) = -conj(bondlist[name].as<complex>());
    }
    arma::vec seigs_up;
    arma::eig_sym(seigs_up, Hs_up);

    // Create single particle matrix for dnspins
    arma::cx_mat Hs_dn(n_sites, n_sites, arma::fill::zeros);
    for (auto bond : bondlist.bonds_of_type("HOPDN")) {
      assert(bond.size() == 2);
      int s1 = bond.site(0);
      int s2 = bond.site(1);
      auto name = bond.coupling_name();
      Hs_dn(s1, s2) = -bondlist[name].as<complex>();
      Hs_dn(s2, s1) = -conj(bondlist[name].as<complex>());
    }
    arma::vec seigs_dn;
    arma::eig_sym(seigs_dn, Hs_dn);

    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {

        // Compute exact gs energy
        double e0_exact = 0;
        for (int i = 0; i < nup; ++i)
          e0_exact += seigs_up(i);
        for (int i = 0; i < ndn; ++i)
          e0_exact += seigs_dn(i);

        auto block3 = Electron(n_sites, nup, ndn);
        auto H = matrix_cplx(bondlist, block3, block3);
        REQUIRE(H.is_hermitian(1e-8));
        arma::vec evecs;
        arma::eig_sym(evecs, H);
        double e0 = evecs(0);
        // printf("nup: %d, ndn: %d, dim: %d, e0: %f, e0_exact: %f\n", nup, ndn,
        //        (int)evecs.size(), e0, e0_exact);
        REQUIRE(close(e0_exact, e0));
      }
  }

  // Test Heisenberg terms at half-filling
  for (int n_sites = 2; n_sites <= 6; ++n_sites) {
    Log("electron_matrix: Heisenberg all-to-all comparison test, N={}",
        n_sites);
    int nup = n_sites / 2;
    int ndn = n_sites - nup;
    auto block_spinhalf = Spinhalf(n_sites, nup);
    auto block_electron = Electron(n_sites, nup, ndn);

    auto bonds = testcases::spinhalf::HB_alltoall(n_sites);
    bonds["U"] = 999999; // gap out doubly occupied sites
    auto H_spinhalf = matrix_real(bonds, block_spinhalf, block_spinhalf);
    auto H_electron = matrix_real(bonds, block_electron, block_electron);
    REQUIRE(H_spinhalf.is_hermitian(1e-8));
    REQUIRE(H_electron.is_hermitian(1e-8));

    arma::vec eigs_spinhalf;
    arma::eig_sym(eigs_spinhalf, H_spinhalf);

    arma::vec eigs_electron;
    arma::eig_sym(eigs_electron, H_electron);

    for (uint64_t idx = 0; idx < eigs_spinhalf.size(); ++idx)
      REQUIRE(std::abs(eigs_spinhalf(idx) - eigs_electron(idx)) < 1e-5);
  }

  /////////////////////////////////////////////////////
  // Henry's MATLAB code test (tests Heisenberg terms)
  {
    Log("electron_matrix: U-hopping-HB full spectrum of Henry's Matlab code");

    auto [bondlist, eigs_correct] = randomAlltoAll4NoU();

    // Compute full spectrum in hydra
    n_sites = 4;
    std::vector<double> all_eigs;
    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {
        auto block = Electron(n_sites, nup, ndn);
        auto H = matrix_real(bondlist, block, block);
        REQUIRE(H.is_hermitian(1e-8));
        arma::vec eigs;
        arma::eig_sym(eigs, H);
        for (auto eig : eigs)
          all_eigs.push_back(eig);
      }
    std::sort(all_eigs.begin(), all_eigs.end());
    REQUIRE(close(arma::vec(all_eigs), eigs_correct));

    std::tie(bondlist, eigs_correct) = randomAlltoAll4();
    all_eigs.clear();
    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {

        auto block = Electron(n_sites, nup, ndn);
        auto H = matrix_real(bondlist, block, block);
        REQUIRE(H.is_hermitian(1e-8));
        arma::vec eigs;
        arma::eig_sym(eigs, H);

        for (auto eig : eigs)
          all_eigs.push_back(eig);
      }
    std::sort(all_eigs.begin(), all_eigs.end());
    REQUIRE(close(arma::vec(all_eigs), eigs_correct));
  }

  for (int N = 3; N <= 6; ++N) {
    Log.out("electron_matrix: random all-to-all complex exchange test, N={}",
            N);

    auto bonds = hydra::testcases::tj::tj_alltoall_complex(N);

    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N - nup; ++ndn) {
        auto block = Electron(N, nup, ndn);
        auto H = matrix_cplx(bonds, block, block);
        REQUIRE(H.is_hermitian(1e-8));
      }

    // Set hoppings to zero
    for (int s1 = 0; s1 < N; ++s1)
      for (int s2 = s1 + 1; s2 < N; ++s2) {
        std::stringstream ss;
        ss << "T" << s1 << "_" << s2;
        std::string name = ss.str();
        bonds[name] = 0.;
      }
    bonds["U"] = 1000;

    // Check whether eigenvalues agree with HB model
    for (int nup = 0; nup <= N; ++nup) {
      int ndn = N - nup;
      auto block1 = Electron(N, nup, ndn);
      auto block2 = Spinhalf(N, nup);
      auto H1 = matrix_cplx(bonds, block1, block1);
      auto H2 = matrix_cplx(bonds, block2, block2);
      arma::vec eigs1;
      arma::eig_sym(eigs1, H1);

      arma::vec eigs2;
      arma::eig_sym(eigs2, H2);

      // Log("eigs1(0): {}, eigs2(0): {}", eigs1(0), eigs2(0));
      // Log("eigs1.size(): {}, eigs2.size(): {}", eigs1.size(), eigs2.size());

      arma::vec eigs1_sub = eigs1.subvec(0, eigs2.size() - 1);
      REQUIRE(close(eigs1_sub, eigs2, 1e-6, 1e-6));
    }
  }

  for (int N = 3; N <= 4; ++N) {
    Log.out("electron_matrix: random all-to-all complex exchange test Np "
            "<-> NoNp, N={}",
            N);
    auto bonds = hydra::testcases::tj::tj_alltoall_complex(N);
    test_electron_np_no_np_matrix(N, bonds);
  }
}
