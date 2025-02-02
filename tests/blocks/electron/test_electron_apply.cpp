#include "../../catch.hpp"

#include <iostream>

#include "../tj/testcases_tj.hpp"
#include "testcases_electron.hpp"
#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/algebra/isapprox.hpp>

using namespace xdiag;

void test_electron_np_no_np_apply(int nsites, OpSum ops) {

  auto block_full = Electron(nsites);
  auto e0_full = eigval0(ops, block_full);

  std::vector<double> e0s;
  for (int nup = 0; nup <= nsites; ++nup)
    for (int ndn = 0; ndn <= nsites; ++ndn) {
      auto block = Electron(nsites, nup, ndn);
      auto e0 = eigval0(ops, block);
      e0s.push_back(e0);
    }
  auto e0_np = *std::min_element(e0s.begin(), e0s.end());
  REQUIRE(isapprox(e0_full, e0_np));
}

TEST_CASE("electron_apply", "[electron]") {
  using namespace xdiag::testcases::electron;

  OpSum ops;

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (real)
  for (int nsites = 3; nsites < 7; ++nsites) {

    Log("electron_apply: Hubbard random all-to-all test (real), N: {}",
        nsites);
    ops = freefermion_alltoall(nsites);
    ops["U"] = 5.0;

    for (int nup = 0; nup <= nsites; ++nup)
      for (int ndn = 0; ndn <= nsites; ++ndn) {

        // Create block and matrix for comparison
        auto block = Electron(nsites, nup, ndn);
        auto H = matrix(ops, block, block);
        REQUIRE(H.is_hermitian(1e-8));

        // Check whether apply gives the same as matrix multiplication
        arma::cx_vec v(block.size(), arma::fill::randn);
        arma::cx_vec w1 = H * v;
        arma::cx_vec w2(block.size(), arma::fill::zeros);
        apply(ops, block, v, block, w2);
        REQUIRE(isapprox(w1, w2));
        arma::mat m(block.size(), 5, arma::fill::randn);
        arma::mat n1 = H * m;
        arma::mat n2(block.size(), 5, arma::fill::zeros);
        apply(ops, block, m, block, n2);
        REQUIRE(isapprox(n1, n2));

        // Compute eigenvalues and compare
        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);
        double e0_mat = evals_mat(0);
        double e0_app = eigval0(ops, block);
        // Log.out("nup: {}, ndn: {}, e0_mat: {}, e0_app: {}", nup, ndn, e0_mat,
        //         e0_app);
        // REQUIRE(isapprox(e0_mat, e0_app));
	REQUIRE(std::abs(e0_mat - e0_app)<1e-8);
      }
  }

  // /////////////////
  // Test Fermion all to all, free fermions (cplx, up/dn different)
  for (int nsites = 3; nsites < 7; ++nsites) {
    Log("electron_apply: Hubbard random all-to-all test (cplx), N: {}",
        nsites);
    ops = freefermion_alltoall_complex_updn(nsites);
    ops["U"] = 5.0;

    for (int nup = 0; nup <= nsites; ++nup)
      for (int ndn = 0; ndn <= nsites; ++ndn) {

        // Create block and matrix for comparison
        auto block = Electron(nsites, nup, ndn);
        auto H = matrixC(ops, block, block);
        REQUIRE(H.is_hermitian(1e-8));

        // Check whether apply gives the same as matrix multiplication
        arma::cx_vec v(block.size(), arma::fill::randn);
        arma::cx_vec w1 = H * v;
        arma::cx_vec w2(block.size(), arma::fill::zeros);
        apply(ops, block, v, block, w2);
        REQUIRE(isapprox(w1, w2));

        // Compute eigenvalues and compare
        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);
        double e0_mat = evals_mat(0);
        double e0_app = eigval0(ops, block);
        // Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
        CHECK(isapprox(e0_mat, e0_app));
      }
  }

  ////////////////////////////
  // Henry's MATLAB code test (tests Heisenberg terms)
  Log("electron_apply: U-hopping-HB apply of Henry's Matlab code");
  {
    auto [ops, eigs_correct] = randomAlltoAll4NoU();

    int nsites = 4;
    for (int nup = 0; nup <= nsites; ++nup)
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto block = Electron(nsites, nup, ndn);
        auto H = matrix(ops, block, block);
        REQUIRE(H.is_hermitian(1e-8));

        // Check whether apply gives the same as matrix multiplication
        arma::vec v(block.size(), arma::fill::randn);
        arma::vec w1 = H * v;
        arma::vec w2(block.size(), arma::fill::zeros);
        apply(ops, block, v, block, w2);
        REQUIRE(isapprox(w1, w2));

        // Compute eigenvalues and compare
        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);
        double e0_mat = evals_mat(0);
        double e0_app = eigval0(ops, block);
        // Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
        CHECK(isapprox(e0_mat, e0_app));
      }

    std::tie(ops, eigs_correct) = randomAlltoAll4();
    for (int nup = 0; nup <= nsites; ++nup)
      for (int ndn = 0; ndn <= nsites; ++ndn) {
        auto block = Electron(nsites, nup, ndn);
        auto H = matrix(ops, block, block);
        REQUIRE(H.is_hermitian(1e-8));

        // Check whether apply gives the same as matrix multiplication
        arma::vec v(block.size(), arma::fill::randn);
        arma::vec w1 = H * v;
        arma::vec w2(block.size(), arma::fill::zeros);
        apply(ops, block, v, block, w2);
        REQUIRE(isapprox(w1, w2));

        // Compute eigenvalues and compare
        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);
        double e0_mat = evals_mat(0);
        double e0_app = eigval0(ops, block);
        // Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
        CHECK(isapprox(e0_mat, e0_app, 1e-6, 1e-6));
      }
  }

  for (int N = 3; N <= 5; ++N) {
    Log.out("electron_apply: random all-to-all complex exchange test Np "
            "<-> NoNp, N={}",
            N);
    auto ops = xdiag::testcases::tj::tj_alltoall_complex(N);
    test_electron_np_no_np_apply(N, ops);
  }
}
