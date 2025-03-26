#include "../../catch.hpp"

#include <iostream>

#include "testcases_tj.hpp"
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/algebra/isapprox.hpp>

#include <xdiag/all.hpp>

using namespace xdiag;

TEST_CASE("tj_apply", "[tj]") {
  using namespace xdiag::testcases::tj;

  for (int N = 3; N <= 6; ++N) {
    Log("tj_apply: random all-to-all real apply=matrix, N={}", N);

    auto ops = tj_alltoall(N);
    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N - nup; ++ndn) {

        auto block = tJ(N, nup, ndn);
        auto H = matrix(ops, block, block);
        REQUIRE(arma::norm(H - H.t()) < 1e-12);

        arma::vec v(block.size(), arma::fill::randn);
        arma::vec w1 = H * v;
        arma::vec w2(block.size(), arma::fill::zeros);
        apply(ops, block, v, block, w2);
        REQUIRE(isapprox(w1, w2));
        arma::mat m(block.size(), 5, arma::fill::randn);
        arma::mat n1 = H * m;
        arma::mat n2(block.size(), 5, arma::fill::zeros);
        apply(ops, block, m, block, n2);
        REQUIRE(isapprox(n1, n2));

        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);
        double e0_mat = evals_mat(0);
        double e0_app = eigval0(ops, block);
        // Log.out("nup: {}, ndn: {}, e0_mat: {}, e0_app: {}", nup, ndn,
        //              e0_mat, e0_app);
        REQUIRE(std::abs(e0_mat - e0_app) < 1e-10);
      }
  }

  for (int N = 3; N <= 6; ++N) {
    Log("tj_apply: random all-to-all complex apply=matrix, N={}", N);
    auto ops = tj_alltoall_complex(N);
 
    for (int nup = 0; nup <= N; ++nup)
      for (int ndn = 0; ndn <= N - nup; ++ndn) {
        auto block = tJ(N, nup, ndn);
        auto H = matrixC(ops, block, block);
        REQUIRE(arma::norm(H - H.t()) < 1e-12);
        arma::cx_vec v(block.size(), arma::fill::randn);
        arma::cx_vec w1 = H * v;
        arma::cx_vec w2(block.size(), arma::fill::zeros);
        apply(ops, block, v, block, w2);
        REQUIRE(isapprox(w1, w2));

        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);

        double e0_mat = evals_mat(0);
        double e0_app = eigval0(ops, block);

        // Log.out("nup: {}, ndn: {}, e0_mat: {}, e0_app: {}", nup, ndn,
        //              e0_mat, e0_app);
        REQUIRE(isapprox(e0_mat, e0_app));
      }
  }

  
  // Test corrs
  for (int nsites = 2; nsites < 5; ++nsites) {
    Log("tj_apply: corrs, N={}", nsites);

    for (int nup = 0; nup <= nsites; ++nup) {
      for (int ndn = 0; ndn <= nsites - nup; ++ndn) {
        auto b = tJ(nsites, nup, ndn);
        auto r = random_state(b);

        for (int i = 0; i < nsites; ++i) {
          auto a = apply(Op("Nup", i), r) + apply(Op("Ndn", i), r);
          auto b = apply(Op("Ntot", i), r);
          REQUIRE(isapprox(a, b));

          a = 0.5 * (apply(Op("Nup", i), r) - apply(Op("Ndn", i), r));
          b = apply(Op("Sz", i), r);
          REQUIRE(isapprox(a, b));

          for (int j = 0; j < nsites; ++j) {

            a = apply(Op("SzSz", {i, j}), r);
            b = apply(Op("Sz", i), apply(Op("Sz", j), r));
            REQUIRE(isapprox(a, b));

            a = apply(Op("NtotNtot", {i, j}), r);
            b = apply(Op("Ntot", i), apply(Op("Ntot", j), r));
            REQUIRE(isapprox(a, b));
          }
        }
      }
    }
  }
}
