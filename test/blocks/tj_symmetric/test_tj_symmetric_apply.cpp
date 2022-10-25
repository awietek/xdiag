#include "../../catch.hpp"

#include <iostream>

#include "../electron/testcases_electron.h"
#include "../tj/testcases_tj.h"

#include <hydra/all.h>

using namespace hydra;

void test_apply_tj_symmetric(BondList bondlist, PermutationGroup space_group,
                             std::vector<Representation> irreps) {
  int n_sites = space_group.n_sites();

  for (int nup = 0; nup <= n_sites; ++nup) {
    for (int ndn = 0; ndn <= n_sites; ++ndn) {

      if (nup + ndn > n_sites)
        continue;

      for (int k = 0; k < (int)irreps.size(); ++k) {
        auto irrep = irreps[k];
        auto block = tJ(n_sites, nup, ndn, space_group, irrep);

        if (block.size() > 0) {
          auto H_sym = matrix_cplx(bondlist, block, block);
          REQUIRE(arma::norm(H_sym - H_sym.t()) < 1e-12);

          // Check whether apply gives the same as matrix multiplication
          arma::cx_vec v(block.size(), arma::fill::randn);
          arma::cx_vec w1 = H_sym * v;
          arma::cx_vec w2(block.size(), arma::fill::zeros);
          apply(bondlist, block, v, block, w2);
          REQUIRE(close(w1, w2));

          // Compute eigenvalues and compare
          arma::vec evals_mat;
          arma::eig_sym(evals_mat, H_sym);
          double e0_mat = evals_mat(0);
          double e0_app = eig0_cplx(bondlist, block);

          // Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
          REQUIRE(std::abs(e0_mat - e0_app) < 1e-7);

          // Compute eigenvalues with real arithmitic
          if (is_real(block.irrep()) && bondlist.is_real()) {
            auto H_real = matrix_real(bondlist, block, block);
            REQUIRE(arma::norm(H_real - H_real.t()) < 1e-12);

            arma::vec evals_mat_real;
            arma::eig_sym(evals_mat_real, H_real);
            REQUIRE(close(evals_mat_real, evals_mat));

            double e0_mat_real = evals_mat_real(0);
            double e0_app_real = eig0_real(bondlist, block);
            REQUIRE(std::abs(e0_mat_real - e0_app_real) < 1e-7);
          }
        }
      }
    }
  }
}

void test_tj_symmetric_apply_chains(int n_sites) {
  using namespace hydra::testcases::tj;
  using namespace hydra::testcases::electron;

  Log("tj_symmetric_apply: tJ chain, symmetric apply test, n_sites: {}",
      n_sites);
  auto bondlist = tJchain(n_sites, 1.0, 5.0);
  auto [space_group, irreps, multiplicities] =
      get_cyclic_group_irreps_mult(n_sites);
  (void)multiplicities;
  test_apply_tj_symmetric(bondlist, space_group, irreps);
}

TEST_CASE("tj_symmetric_apply", "[blocks][tj]") {
  using namespace hydra::testcases::tj;
  using namespace hydra::testcases::electron;

  // Test linear chains
  for (int n_sites = 2; n_sites < 8; ++n_sites) {
    test_tj_symmetric_apply_chains(n_sites);
  }
  {
    // test a 3x3 triangular lattice
    Log("tj_symmetric_apply: tJ 3x3 triangular, symmetric apply test");
    std::string lfile = "data/triangular.9.hop.sublattices.tsl.lat";

    auto bondlist = read_bondlist(lfile);
    bondlist["T"] = 1.0;
    bondlist["U"] = 5.0;
    auto permutations = hydra::read_permutations(lfile);
    auto space_group = PermutationGroup(permutations);

    std::vector<std::pair<std::string, int>> rep_name_mult = {
        {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
        {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
        {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
        {"Y.C1.A", 6}};

    std::vector<Representation> irreps;
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(read_represenation(lfile, name));
      (void)mult;
    }
    test_apply_tj_symmetric(bondlist, space_group, irreps);
  }

  {
    // test a 3x3 triangular lattice with complex flux
    Log("tj_symmetric_apply: tJ 3x3 triangular staggered, symmetric apply "
        "test, complex");
    std::string lfile =
        "data/triangular.9.tup.phi.tdn.nphi.sublattices.tsl.lat";

    auto bondlist = read_bondlist(lfile);
    std::vector<double> etas{0.0, 0.1, 0.2, 0.3};
    auto permutations = hydra::read_permutations(lfile);
    auto space_group = PermutationGroup(permutations);

    std::vector<std::pair<std::string, int>> rep_name_mult = {
        {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
        {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
        {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
        {"Y.C1.A", 6}};

    std::vector<Representation> irreps;
    std::vector<int> multiplicities;
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(read_represenation(lfile, name));
      multiplicities.push_back(mult);
    }

    for (auto eta : etas) {
      bondlist["TPHI"] = complex(cos(eta * M_PI), sin(eta * M_PI));
      bondlist["JPHI"] = complex(cos(2 * eta * M_PI), sin(2 * eta * M_PI));

      test_apply_tj_symmetric(bondlist, space_group, irreps);
    }
  }
}
