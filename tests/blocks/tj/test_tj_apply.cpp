#include "../../catch.hpp"

#include <iostream>

#include "testcases_tj.hpp"
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/blocks/spinhalf/spinhalf_matrix.hpp>
#include <xdiag/blocks/tj/tj_apply.hpp>
#include <xdiag/blocks/tj/tj_matrix.hpp>
#include <xdiag/utils/close.hpp>

using namespace xdiag;

TEST_CASE("tj_apply", "[tj]") {
  using namespace xdiag::testcases::tj;

  for (int N = 3; N <= 6; ++N) {
    Log("tj_apply: random all-to-all real apply=matrix, N={}", N);

    auto bonds = tj_alltoall(N);
    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N - nup; ++ndn) {

        auto block = tJ(N, nup, ndn);
        auto H = matrix(bonds, block, block);
        REQUIRE(arma::norm(H - H.t()) < 1e-12);

        arma::vec v(block.size(), arma::fill::randn);
        arma::vec w1 = H * v;
        arma::vec w2(block.size(), arma::fill::zeros);
        apply(bonds, block, v, block, w2);
        REQUIRE(close(w1, w2));

        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);
        double e0_mat = evals_mat(0);
        double e0_app = eigval0(bonds, block);
        // Log.out("nup: {}, ndn: {}, e0_mat: {}, e0_app: {}", nup, ndn,
        //              e0_mat, e0_app);
        REQUIRE(std::abs(e0_mat - e0_app) < 1e-10);
      }
  }

  for (int N = 3; N <= 6; ++N) {
    Log("tj_apply: random all-to-all complex apply=matrix, N={}", N);

    auto bonds = tj_alltoall_complex(N);
    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N - nup; ++ndn) {
        auto block = tJ(N, nup, ndn);
        auto H = matrixC(bonds, block, block);
        REQUIRE(arma::norm(H - H.t()) < 1e-12);
        arma::cx_vec v(block.size(), arma::fill::randn);
        arma::cx_vec w1 = H * v;
        arma::cx_vec w2(block.size(), arma::fill::zeros);
        apply(bonds, block, v, block, w2);
        REQUIRE(close(w1, w2));

        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);

        double e0_mat = evals_mat(0);
        double e0_app = eigval0(bonds, block);

        // Log.out("nup: {}, ndn: {}, e0_mat: {}, e0_app: {}", nup, ndn,
        //              e0_mat, e0_app);
        REQUIRE(close(e0_mat, e0_app));
      }
  }
}
