#include "../../catch.hpp"

#include <iostream>

#include "../tj/testcases_tj.h"
#include "testcases_electron.h"
#include <hydra/all.h>

using namespace hydra;

void test_electron_np_no_np_apply(int n_sites, BondList bonds) {

  auto block_full = Electron(n_sites);
  auto e0_full = e0_cplx(bonds, block_full);

  std::vector<double> e0s;
  for (int nup = 0; nup <= n_sites; ++nup)
    for (int ndn = 0; ndn <= n_sites; ++ndn) {
      auto block = Electron(n_sites, nup, ndn);
      auto e0 = e0_cplx(bonds, block);
      e0s.push_back(e0);
    }
  auto e0_np = *std::min_element(e0s.begin(), e0s.end());
  REQUIRE(close(e0_full, e0_np));
}

TEST_CASE("electron_apply", "[blocks][electron]") {
  using namespace hydra::testcases::electron;

  BondList bondlist;

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (real)
  for (int n_sites = 3; n_sites < 7; ++n_sites) {

    Log("electron_apply: Hubbard random all-to-all test (real), N: {}",
        n_sites);
    bondlist = freefermion_alltoall(n_sites);
    bondlist["U"] = 5.0;

    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {

        // Create block and matrix for comparison
        auto block = Electron<uint32_t>(n_sites, nup, ndn);
        auto H = matrix_real(bondlist, block, block);
        REQUIRE(H.is_hermitian(1e-8));

        // Check whether apply gives the same as matrix multiplication
        arma::cx_vec v(block.size(), arma::fill::randn);
        arma::cx_vec w1 = H * v;
        arma::cx_vec w2(block.size(), arma::fill::zeros);
        apply(bondlist, block, v, block, w2);
        REQUIRE(close(w1, w2));

        // Compute eigenvalues and compare
        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);
        double e0_mat = evals_mat(0);
        double e0_app = e0_real(bondlist, block);
        // Log.out("nup: {}, ndn: {}, e0_mat: {}, e0_app: {}", nup, ndn, e0_mat,
        //         e0_app);
        // REQUIRE(close(e0_mat, e0_app));
	REQUIRE(std::abs(e0_mat - e0_app)<1e-8);
      }
  }

  // /////////////////
  // Test Fermion all to all, free fermions (cplx, up/dn different)
  for (int n_sites = 3; n_sites < 7; ++n_sites) {
    Log("electron_apply: Hubbard random all-to-all test (cplx), N: {}",
        n_sites);
    bondlist = freefermion_alltoall_complex_updn(n_sites);
    bondlist["U"] = 5.0;

    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {

        // Create block and matrix for comparison
        auto block = Electron<uint32_t>(n_sites, nup, ndn);
        auto H = matrix_cplx(bondlist, block, block);
        REQUIRE(H.is_hermitian(1e-8));

        // Check whether apply gives the same as matrix multiplication
        arma::cx_vec v(block.size(), arma::fill::randn);
        arma::cx_vec w1 = H * v;
        arma::cx_vec w2(block.size(), arma::fill::zeros);
        apply(bondlist, block, v, block, w2);
        REQUIRE(close(w1, w2));

        // Compute eigenvalues and compare
        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);
        double e0_mat = evals_mat(0);
        double e0_app = e0_cplx(bondlist, block);
        // Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
        CHECK(close(e0_mat, e0_app));
      }
  }

  ////////////////////////////
  // Henry's MATLAB code test (tests Heisenberg terms)
  Log("electron_apply: U-hopping-HB apply of Henry's Matlab code");
  {
    auto [bondlist, eigs_correct] = randomAlltoAll4NoU();

    int n_sites = 4;
    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {
        auto block = Electron(n_sites, nup, ndn);
        auto H = matrix_real(bondlist, block, block);
        REQUIRE(H.is_hermitian(1e-8));

        // Check whether apply gives the same as matrix multiplication
        arma::vec v(block.size(), arma::fill::randn);
        arma::vec w1 = H * v;
        arma::vec w2(block.size(), arma::fill::zeros);
        apply(bondlist, block, v, block, w2);
        REQUIRE(close(w1, w2));

        // Compute eigenvalues and compare
        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);
        double e0_mat = evals_mat(0);
        double e0_app = e0_cplx(bondlist, block);
        // Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
        CHECK(close(e0_mat, e0_app));
      }

    std::tie(bondlist, eigs_correct) = randomAlltoAll4();
    for (int nup = 0; nup <= n_sites; ++nup)
      for (int ndn = 0; ndn <= n_sites; ++ndn) {
        auto block = Electron(n_sites, nup, ndn);
        auto H = matrix_real(bondlist, block, block);
        REQUIRE(H.is_hermitian(1e-8));

        // Check whether apply gives the same as matrix multiplication
        arma::vec v(block.size(), arma::fill::randn);
        arma::vec w1 = H * v;
        arma::vec w2(block.size(), arma::fill::zeros);
        apply(bondlist, block, v, block, w2);
        REQUIRE(close(w1, w2));

        // Compute eigenvalues and compare
        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);
        double e0_mat = evals_mat(0);
        double e0_app = e0_cplx(bondlist, block);
        // Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
        CHECK(close(e0_mat, e0_app, 1e-6, 1e-6));
      }
  }

  for (int N = 3; N <= 5; ++N) {
    Log.out("electron_apply: random all-to-all complex exchange test Np "
            "<-> NoNp, N={}",
            N);
    auto bonds = hydra::testcases::tj::tj_alltoall_complex(N);
    test_electron_np_no_np_apply(N, bonds);
  }
}
