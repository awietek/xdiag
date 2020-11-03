#include "../catch.hpp"

#include <fstream>
#include <hydra/all.h>
#include <lila/all.h>
#include <random>

#include "testcases_hubbardmodel.h"

using namespace hydra;
using namespace lila;

// Test by comparing to exactly known solution for given quantumnumber
template <class coeff_t = double>
void test_hubbardmodel_e0(BondList bondlist, Couplings couplings,
                          qn_electron qn, double e0) {
  auto model = HubbardModel<coeff_t>(bondlist, couplings, qn);
  auto H = model.matrix();
  REQUIRE(lila::close(H, lila::Herm(H)));
  auto eigs = lila::EigenvaluesSym(H);
  // printf("eigs(0): %f, e0 %f\n", eigs(0), e0);
  REQUIRE(std::abs(e0 - eigs(0)) < 1e-6);
}

// Test by comparing to exactly known solution for given quantumnumber
template <class coeff_t = double>
void test_hubbardmodel_heisenberg(BondList bondlist, Couplings couplings,
                                  int max_nup = -1) {
  int n_sites = bondlist.n_sites();

  couplings["U"] = 999; // gap out doubly occupied sites
  for (int nup = 0; nup <= (max_nup < 0 ? n_sites : std::min(max_nup, n_sites));
       ++nup) {
    int ndn = n_sites - nup;
    qn_electron qn = {nup, ndn};

    auto HB = HeisenbergModel<>();
    auto HBM = HB.matrix(bondlist, couplings, nup);

    // Create Hubbard matrix
    auto Hubbard = HubbardModel<double>(bondlist, couplings, qn);
    auto HubbardM = Hubbard.matrix();

    auto ehubbard = lila::EigenvaluesSym(HubbardM);
    auto ehb = lila::EigenvaluesSym(HBM);

    REQUIRE(close(HBM, lila::Herm(HBM)));
    REQUIRE(close(HubbardM, lila::Herm(HubbardM)));

    for (int idx = 0; idx < ehb.size(); ++idx)
      REQUIRE(std::abs(ehubbard(idx) - ehb(idx)) < 1e-10);
  }
}

// Test by comparing to full spectrum of alps
void test_hubbardmodel_fullspectrum(BondList bondlist, Couplings couplings,
                                    std::string filename) {
  int n_sites = bondlist.n_sites();

  // Compute full spectrum in hydra
  lila::Vector<double> all_eigs;
  for (int nup = 0; nup <= n_sites; ++nup)
    for (int ndn = 0; ndn <= n_sites; ++ndn) {
      qn_electron qn = {nup, ndn};
      auto model = HubbardModel<double>(bondlist, couplings, qn);

      // Run Full ED
      auto H = model.matrix();

      if (nup + ndn == n_sites)
        REQUIRE(lila::close(H, lila::Herm(H)));
      auto eigs = lila::EigenvaluesSym(H);
      for (auto eig : eigs)
        all_eigs.push_back(eig);
    }
  std::sort(all_eigs.begin(), all_eigs.end());

  // Read alps_eigs
  lila::Vector<double> alps_eigs;
  std::ifstream in(filename.c_str());
  if (!in) {
    std::cerr << "test_tjmodel.cpp: Cannot open the File : " << filename
              << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string str;
  while (std::getline(in, str)) {
    if (str.size() > 0)
      alps_eigs.push_back(std::stod(str));
  }
  in.close();
  std::sort(alps_eigs.begin(), alps_eigs.end());
  REQUIRE(all_eigs.size() == alps_eigs.size());
  for (int i = 0; i < all_eigs.size(); ++i) {
    REQUIRE(close(all_eigs(i), alps_eigs(i)));
  }
}

// Test if complex matrix is hermitian
void test_hubbardmodel_hermitian(BondList bondlist, Couplings couplings) {
  int n_sites = bondlist.n_sites();
  for (int nup = 0; nup <= n_sites; ++nup)
    for (int ndn = 0; ndn <= n_sites - nup; ++ndn) {
      auto model = HubbardModel<complex>(bondlist, couplings, {nup, ndn});
      auto H = model.matrix();
      REQUIRE(lila::close(H, lila::Herm(H)));
    }
}

TEST_CASE("HubbardModel", "[HubbardModel]") {
  using namespace hydra::hubbardtestcases;
  BondList bondlist;
  Couplings couplings;

  //////////////////////////////////
  // Test two site exact solution
  bondlist << Bond("HUBBARDHOP", "T", {0, 1});
  couplings["T"] = 1.0;
  for (int i = 0; i < 20; ++i) {
    double U = 1.234 * i;
    couplings["U"] = U;
    printf("HubbardModel: two-site exact solution test, U=%f\n", U);
    double e0 = 0.5 * (U - sqrt(U * U + 16));
    test_hubbardmodel_e0(bondlist, couplings, {1, 1}, e0);
  }

  ////////////////////////////////////////////////
  // Cross-checks with various heisenberg models
  printf("Hubbardmodel: Heisenberg triangle test, N=3\n");
  std::tie(bondlist, couplings) = heisenberg_triangle();
  test_hubbardmodel_heisenberg(bondlist, couplings);

  ///////////////////////////////////////////////////
  // Test all-to-all random coupling Heisenberg
  for (int n_sites = 3; n_sites < 8; ++n_sites) {
    printf("Hubbardmodel: Heisenberg random all-to-all test, ");
    printf("N=%d\n", n_sites);
    std::tie(bondlist, couplings) = heisenberg_alltoall(n_sites);
    test_hubbardmodel_heisenberg(bondlist, couplings);
  }

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (real)
  for (int n_sites = 3; n_sites < 8; ++n_sites) {
    printf("HubbardModel: free fermion random all-to-all test, ");
    printf("N=%d\n", n_sites);
    std::tie(bondlist, couplings) = freefermion_alltoall(n_sites);

    // Create single particle matrix
    auto Hs = lila::Zeros<double>(n_sites, n_sites);
    for (auto bond : bondlist) {
      assert(bond.size() == 2);
      int s1 = bond.sites(0);
      int s2 = bond.sites(1);
      auto name = bond.coupling();
      Hs(s1, s2) = -lila::real(couplings[name]);
      Hs(s2, s1) = -lila::real(couplings[name]);
    }
    auto seigs = lila::EigenvaluesSym(Hs);
    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {
        double e0 = 0;
        for (int i = 0; i < nup; ++i)
          e0 += seigs(i);
        for (int i = 0; i < ndn; ++i)
          e0 += seigs(i);
        test_hubbardmodel_e0(bondlist, couplings, {nup, ndn}, e0);
      }
  }

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (complex)
  for (int n_sites = 3; n_sites < 8; ++n_sites) {
    printf("HubbardModel: free fermion random all-to-all test (cplx), ");
    printf("N=%d\n", n_sites);
    std::tie(bondlist, couplings) = freefermion_alltoall_complex(n_sites);

    // Create single particle matrix
    auto Hs = lila::Zeros<complex>(n_sites, n_sites);
    for (auto bond : bondlist) {
      assert(bond.size() == 2);
      int s1 = bond.sites(0);
      int s2 = bond.sites(1);
      auto name = bond.coupling();
      Hs(s1, s2) = -couplings[name];
      Hs(s2, s1) = -lila::conj(couplings[name]);
    }
    auto seigs = lila::EigenvaluesSym(Hs);

    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {
        double e0 = 0;
        for (int i = 0; i < nup; ++i)
          e0 += seigs(i);
        for (int i = 0; i < ndn; ++i)
          e0 += seigs(i);
        test_hubbardmodel_e0<complex>(bondlist, couplings, {nup, ndn}, e0);
      }
  }

  /////////////////////////////////////////////////////////////////
  // Test of full spectrum of random all-to-all interactions
  // by comparing to MATLAB results
  printf("HubbardModel: MATLAB full spectrum, chain N=4, no Hubbard U\n");
  std::tie(bondlist, couplings) = randomAlltoAll4NoU();
  test_hubbardmodel_fullspectrum(
      bondlist, couplings,
      "data/hubbardfullspectrum/spectrum.allToAll.N.4.U.0.txt");

  printf("HubbardModel: MATLAB full spectrum, chain N=4, Hubbard U\n");
  std::tie(bondlist, couplings) = randomAlltoAll4();
  test_hubbardmodel_fullspectrum(
      bondlist, couplings,
      "data/hubbardfullspectrum/spectrum.allToAll.N.4.U.5.txt");
}
