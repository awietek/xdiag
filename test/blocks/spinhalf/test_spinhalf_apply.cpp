#include "../../catch.hpp"

#include <iostream>

#include "testcases_spinhalf.h"
#include <hydra/all.h>

using namespace hydra;

void test_apply(BondList bonds, Couplings couplings) {
  int N = bonds.n_sites();
  for (int nup = 0; nup <= N; ++nup) {
    auto block = Spinhalf<uint32_t>(N, nup);
    auto H = MatrixReal(bonds, couplings, block, block);
    REQUIRE(H.is_hermitian(1e-8));

    arma::vec v(block.size(), arma::fill::randn);
    arma::vec w1 = H * v;
    arma::vec w2(block.size(), arma::fill::zeros);
    Apply(bonds, couplings, block, v, block, w2);
    REQUIRE(close(w1, w2));

    arma::vec evals_mat;
    arma::eig_sym(evals_mat, H);

    double e0_mat = evals_mat(0);
    double e0_app = E0Real(bonds, couplings, block);
    REQUIRE(close(e0_mat, e0_app));
  }
}

TEST_CASE("spinhalf_apply", "[models][spinhalf]") {
  using namespace hydra::testcases::spinhalf;

  Log.out("spinhalf_apply: Heisenberg chain apply test, J=1.0, N=2,..,6");
  for (int N = 2; N <= 6; ++N) {
    auto [bonds, couplings] = HBchain(N, 1.0);
    test_apply(bonds, couplings);
  }

  Log.out("spinhalf_apply: Heisenberg alltoall apply test, N=2,..,6");
  for (int N = 2; N <= 6; ++N) {
    auto [bonds, couplings] = HB_alltoall(N);
    test_apply(bonds, couplings);
  }

  Log.out("spinhalf_apply: Heisenberg all-to-all Sz <-> NoSz comparison");
  for (int n_sites = 2; n_sites <= 6; ++n_sites) {
    auto [bonds, couplings] = HB_alltoall(n_sites);
    auto block_no_sz = Spinhalf(n_sites);
    auto e0_no_sz = E0Real(bonds, couplings, block_no_sz);
    auto e0s_sz = std::vector<double>();

    for (int nup = 0; nup <= n_sites; ++nup) {
      auto block_sz = Spinhalf(n_sites, nup);
      auto e0_sz = E0Real(bonds, couplings, block_sz);
      e0s_sz.push_back(e0_sz);
    }
    auto e0_sz = *std::min_element(e0s_sz.begin(), e0s_sz.end());
    REQUIRE(close(e0_sz, e0_no_sz));
  }

  {
    Log("spinhalf_matrix: Triangular J1J2Jchi N=12");
    std::string lfile = "data/triangular.j1j2jch/"
                        "triangular.12.j1j2jch.sublattices.fsl.lat";

    auto bondlist = read_bondlist(lfile);
    Couplings couplings;
    couplings["J1"] = 1.00;
    couplings["J2"] = 0.15;
    couplings["Jchi"] = 0.09;

    int n_sites = 12;
    int n_up = 6;
    auto spinhalf = Spinhalf<uint16_t>(n_sites, n_up);
    auto e0 = E0Cplx(bondlist, couplings, spinhalf);
    double energy = -6.9456000700824329641;

    // Log("{:.18f} {:.18f}", e0, energy);

    REQUIRE(close(e0, energy));
  }
}
