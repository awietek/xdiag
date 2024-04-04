#include "../../catch.hpp"

#include <iostream>

#include "testcases_spinhalf.h"
#include <xdiag/algebra/algebra.h>
#include <xdiag/algebra/matrix.h>
#include <xdiag/algorithms/sparse_diag.h>
#include <xdiag/blocks/spinhalf/spinhalf_apply.h>
#include <xdiag/utils/close.h>
#include <xdiag/utils/print_macro.h>

using namespace xdiag;

void test_apply(BondList bonds) {
  int N = bonds.n_sites();
  for (int nup = 1; nup <= N; ++nup) {
    auto block = Spinhalf(N, nup);
    auto H = matrix(bonds, block, block);
    REQUIRE(H.is_hermitian(1e-8));

    arma::vec v(block.size(), arma::fill::randn);
    arma::vec w1 = H * v;
    arma::vec w2(block.size(), arma::fill::zeros);
    apply(bonds, block, v, block, w2);

    arma::vec w3 = H * H * v;
    arma::vec w4(block.size(), arma::fill::zeros);
    apply(bonds, block, w2, block, w4);
    REQUIRE(close(w3, w4));

    arma::vec evals_mat;
    arma::eig_sym(evals_mat, H);

    double e0_mat = evals_mat(0);
    // double e0_app = eigval0(bonds, block);
    auto [e0_app, ev] = eig0(bonds, block);
    // Log("H: {}, nup: {}, mat: {:.5f} app: {:.5f}", N, nup, e0_mat, e0_app);
    REQUIRE(close(e0_mat, e0_app));
  }
}

TEST_CASE("spinhalf_apply", "[spinhalf]") {
  using namespace xdiag::testcases::spinhalf;

  Log.out("spinhalf_apply: Heisenberg chain apply test, J=1.0, N=2,..,6");
  for (int N = 2; N <= 6; ++N) {
    auto bonds = HBchain(N, 1.0);
    test_apply(bonds);
  }

  Log.out("spinhalf_apply: Heisenberg alltoall apply test, N=2,..,6");
  for (int N = 2; N <= 6; ++N) {
    auto bonds = HB_alltoall(N);
    test_apply(bonds);
  }

  Log.out("spinhalf_apply: Heisenberg all-to-all Sz <-> NoSz comparison");
  for (int n_sites = 2; n_sites <= 6; ++n_sites) {
    auto bonds = HB_alltoall(n_sites);
    auto block_no_sz = Spinhalf(n_sites);
    auto e0_no_sz = eigval0(bonds, block_no_sz);
    auto e0s_sz = std::vector<double>();

    for (int nup = 0; nup <= n_sites; ++nup) {
      auto block_sz = Spinhalf(n_sites, nup);
      auto e0_sz = eigval0(bonds, block_sz);
      e0s_sz.push_back(e0_sz);
    }
    auto e0_sz = *std::min_element(e0s_sz.begin(), e0s_sz.end());
    REQUIRE(close(e0_sz, e0_no_sz));
  }

  {
    Log("spinhalf_matrix: Triangular J1J2Jchi N=12");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/triangular.j1j2jch/"
                        "triangular.12.j1j2jch.sublattices.fsl.lat";

    auto bondlist = read_bondlist(lfile);
    bondlist["J1"] = 1.00;
    bondlist["J2"] = 0.15;
    bondlist["Jchi"] = 0.09;

    int n_sites = 12;
    int n_up = 6;
    auto spinhalf = Spinhalf(n_sites, n_up);
    auto e0 = eigval0(bondlist, spinhalf);
    double energy = -6.9456000700824329641;

    // Log("{:.18f} {:.18f}", e0, energy);

    REQUIRE(close(e0, energy));
  }
}
