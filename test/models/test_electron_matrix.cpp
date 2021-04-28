#include "../catch.hpp"

#include <iostream>

#include "testcases_electron.h"
#include <hydra/all.h>

using namespace hydra;

TEST_CASE("electron_matrix", "[models]") {
  using namespace hydra::testcases::electron;

  BondList bondlist;
  Couplings couplings;

  // Compare with matrix from Weisse and Fehske
  int n_sites = 4;
  int n_up = 3;
  int n_dn = 2;
  double t = 1.0;
  double U = 5.0;
  auto block = Electron<uint32>(n_sites, n_up, n_dn);

  for (int i = 0; i < n_sites; ++i)
    bondlist << Bond("HOP", "T", {i, (i + 1) % n_sites});
  couplings["T"] = 1.0;
  couplings["U"] = 5.0;

  auto H1 = matrix_real(bondlist, couplings, block, block);

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
  auto block2 = Electron<uint32>(2, 1, 1);
  for (int i = 0; i < 20; ++i) {
    double U = 1.234 * i;
    printf("Electron: two-site exact solution test, U=%f\n", U);
    couplings["T"] = 1.0;
    couplings["U"] = U;
    double e0_exact = 0.5 * (U - sqrt(U * U + 16));
    auto H = matrix_real(bondlist, couplings, block2, block2);
    REQUIRE(lila::close(H, lila::Herm(H)));
    auto evecs = lila::EigenvaluesSym(H);
    double e0 = evecs(0);
    // printf("e0: %f, e0_exact: %f\n", e0, e0_exact);
    REQUIRE(lila::close(e0_exact, e0));
  }

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (real)
  for (int n_sites = 3; n_sites < 7; ++n_sites) {
    printf("HubbardModel: free fermion random all-to-all test, (real)\n");
    printf("N=%d\n", n_sites);
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

        auto block3 = Electron<uint32>(n_sites, nup, ndn);
        auto Hr = matrix_real(bondlist, couplings, block3, block3);
	REQUIRE(lila::close(Hr, lila::Herm(Hr)));
        auto evecsr = lila::EigenvaluesSym(Hr);
        auto Hc = matrix_cplx(bondlist, couplings, block3, block3);
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
    printf("HubbardModel: free fermion random all-to-all test, (cplx)\n");
    printf("N=%d\n", n_sites);
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

        auto block3 = Electron<uint32>(n_sites, nup, ndn);
        auto H = matrix_cplx(bondlist, couplings, block3, block3);
	REQUIRE(lila::close(H, lila::Herm(H)));
        auto evecs = lila::EigenvaluesSym(H);
        double e0 = evecs(0);
        // printf("nup: %d, ndn: %d, dim: %d, e0: %f, e0_exact: %f\n", nup, ndn,
        //        (int)evecs.size(), e0, e0_exact); 
        REQUIRE(lila::close(e0_exact, e0));
      }
  }

  // /////////////////////////////////////////////////////
  // // Henry's MATLAB code test (need Heisenberg terms)
  // printf("HubbardModel: testing full spectrum of Henry's Matlab code\n");

  // lila::Vector<double> eigs_correct;
  // std::tie(bondlist, couplings, eigs_correct) = randomAlltoAll4NoU();

  // // Compute full spectrum in hydra
  // n_sites = 4;
  // lila::Vector<double> all_eigs;
  // for (int nup = 0; nup <= n_sites; ++nup)
  //   for (int ndn = 0; ndn <= n_sites; ++ndn) {
  //     auto block = Electron<uint16>(n_sites, nup, ndn);
  //     auto H = matrix_real(bondlist, couplings, block, block);
  //     LilaPrint(H);
  //     REQUIRE(lila::close(H, lila::Herm(H)));
  //     auto eigs = lila::EigenvaluesSym(H);
  //     for (auto eig : eigs)
  //       all_eigs.push_back(eig);
  //   }
  // std::sort(all_eigs.begin(), all_eigs.end());  
  // LilaPrint(all_eigs);
  // LilaPrint(eigs_correct);
  // REQUIRE(lila::close(all_eigs, eigs_correct));

  // std::tie(bondlist, couplings, eigs_correct) = randomAlltoAll4();
  // all_eigs.clear();
  // for (int nup = 0; nup <= n_sites; ++nup)
  //   for (int ndn = 0; ndn <= n_sites; ++ndn) {
  //     auto block = Electron<uint16>(n_sites, nup, ndn);
  //     auto H = matrix_real(bondlist, couplings, block, block);
  //     REQUIRE(lila::close(H, lila::Herm(H)));
  //     auto eigs = lila::EigenvaluesSym(H);
  //     for (auto eig : eigs)
  //       all_eigs.push_back(eig);
  //   }
  // std::sort(all_eigs.begin(), all_eigs.end());  

  // REQUIRE(lila::close(all_eigs, eigs_correct));

}
