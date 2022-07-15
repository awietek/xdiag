#include "../../catch.hpp"

#include <iostream>

#include "testcases_electron.h"
#include "../tj/testcases_tj.h"
#include <hydra/all.h>

using namespace hydra;

void test_electron_np_no_np_apply(int n_sites, BondList bonds, Couplings cpls) {

  auto block_full = Electron(n_sites);
  auto e0_full = E0Cplx(bonds, cpls, block_full);
 

  lila::Vector<double> e0s;
  for (int nup = 0; nup <= n_sites; ++nup)
    for (int ndn = 0; ndn <= n_sites; ++ndn) {
      auto block = Electron(n_sites, nup, ndn);
      auto e0 = E0Cplx(bonds, cpls, block);
      e0s.push_back(e0);
    }
  auto e0_np = *std::min_element(e0s.begin(), e0s.end());
  REQUIRE(lila::close(e0_full, e0_np));
}

TEST_CASE("electron_apply", "[blocks][electron]") {
  using namespace hydra::testcases::electron;

  BondList bondlist;
  Couplings couplings;

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (real)
  for (int n_sites = 3; n_sites < 7; ++n_sites) {

    Log("electron_apply: Hubbard random all-to-all test (real), N: {}",
              n_sites);
    std::tie(bondlist, couplings) = freefermion_alltoall(n_sites);
    couplings["U"] = 5.0;

    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {

        // Create block and matrix for comparison
        auto block = Electron<uint32_t>(n_sites, nup, ndn);
        auto H = MatrixReal(bondlist, couplings, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));

        // Check whether apply gives the same as matrix multiplication
        auto v = lila::Random<double>(block.size());
        auto w1 = lila::Mult(H, v);
        auto w2 = lila::ZerosLike(v);
        Apply(bondlist, couplings, block, v, block, w2);
        REQUIRE(lila::close(w1, w2));

        // Compute eigenvalues and compare
        auto evals_mat = lila::EigenvaluesSym(H);
        double e0_mat = evals_mat(0);
        double e0_app = E0Real(bondlist, couplings, block);
        // Log.out("nup: {}, ndn: {}, e0_mat: {}, e0_app: {}", nup, ndn,
        // e0_mat, e0_app);
        CHECK(lila::close(e0_mat, e0_app, 1e-6, 1e-6));
      }
  }

  // /////////////////
  // Test Fermion all to all, free fermions (cplx, up/dn different)
  for (int n_sites = 3; n_sites < 7; ++n_sites) {
    Log("electron_apply: Hubbard random all-to-all test (cplx), N: {}",
              n_sites);
    std::tie(bondlist, couplings) = freefermion_alltoall_complex_updn(n_sites);
    couplings["U"] = 5.0;

    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {

        // Create block and matrix for comparison
        auto block = Electron<uint32_t>(n_sites, nup, ndn);
        auto H = MatrixCplx(bondlist, couplings, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));

        // Check whether apply gives the same as matrix multiplication
        auto v = lila::Random<complex>(block.size());
        auto w1 = lila::Mult(H, v);
        auto w2 = lila::ZerosLike(v);
        Apply(bondlist, couplings, block, v, block, w2);
        REQUIRE(lila::close(w1, w2));

        // Compute eigenvalues and compare
        auto evals_mat = lila::EigenvaluesSym(H);
        double e0_mat = evals_mat(0);
        double e0_app = E0Cplx(bondlist, couplings, block);
        // Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
        CHECK(lila::close(e0_mat, e0_app, 1e-6, 1e-6));
      }
  }

  /////////////////////////////////////////////////////
  // Henry's MATLAB code test (tests Heisenberg terms)
  Log("electron_apply: U-hopping-HB apply of Henry's Matlab code");
  {
    auto [bondlist, couplings, eigs_correct] = randomAlltoAll4NoU();

    int n_sites = 4;
    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {
        auto block = Electron(n_sites, nup, ndn);
        auto H = MatrixReal(bondlist, couplings, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));

        // Check whether apply gives the same as matrix multiplication
        auto v = lila::Random<double>(block.size());
        auto w1 = lila::Mult(H, v);
        auto w2 = lila::ZerosLike(v);
        Apply(bondlist, couplings, block, v, block, w2);
        REQUIRE(lila::close(w1, w2));

        // Compute eigenvalues and compare
        auto evals_mat = lila::EigenvaluesSym(H);
        double e0_mat = evals_mat(0);
        double e0_app = E0Cplx(bondlist, couplings, block);
        // Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
        CHECK(lila::close(e0_mat, e0_app, 1e-6, 1e-6));
      }

    std::tie(bondlist, couplings, eigs_correct) = randomAlltoAll4();
    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {
        auto block = Electron(n_sites, nup, ndn);
        auto H = MatrixReal(bondlist, couplings, block, block);
        REQUIRE(lila::close(H, lila::Herm(H)));

        // Check whether apply gives the same as matrix multiplication
        auto v = lila::Random<double>(block.size());
        auto w1 = lila::Mult(H, v);
        auto w2 = lila::ZerosLike(v);
        Apply(bondlist, couplings, block, v, block, w2);
        REQUIRE(lila::close(w1, w2));

        // Compute eigenvalues and compare
        auto evals_mat = lila::EigenvaluesSym(H);
        double e0_mat = evals_mat(0);
        double e0_app = E0Cplx(bondlist, couplings, block);
        // Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
        CHECK(lila::close(e0_mat, e0_app, 1e-6, 1e-6));
      }
  }

  for (int N = 3; N <= 5; ++N) {
    Log.out("electron_apply: random all-to-all complex exchange test Np "
                  "<-> NoNp, N={}",
                  N);
    auto [bonds, cpls] = hydra::testcases::tj::tj_alltoall_complex(N);
    test_electron_np_no_np_apply(N, bonds, cpls);
  }
}
