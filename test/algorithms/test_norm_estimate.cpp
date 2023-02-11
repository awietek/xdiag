#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

#include "../blocks/electron/testcases_electron.h"
#include "../blocks/spinhalf/testcases_spinhalf.h"

using namespace hydra;

void test_operator_norm_real(Block const &block, BondList const &bonds) {
  using namespace hydra;
  using namespace arma;

  auto H = matrix_real(bonds, block, block);
  double norm_exact = norm(H, 1);
  auto apply_H = [&H](vec const &v) { return vec(H * v); };
  double norm_est = norm_estimate_real(apply_H, apply_H, block.size());
  double ratio = norm_est / norm_exact;

  // Log("norm_exact: {}", norm_exact);
  // Log("norm_est: {}", norm_est);
  if (norm_exact > 1e-12) {
    // Log("ratio: {}", ratio);
    REQUIRE(((ratio <= 1.00001) && (ratio > 0.3)));
  }
}

void test_operator_norm_cplx(Block const &block, BondList const &bonds) {
  using namespace hydra;
  using namespace arma;

  auto H = matrix_cplx(bonds, block, block);
  double norm_exact = norm(H, 1);
  auto apply_H = [&H](cx_vec const &v) { return cx_vec(H * v); };
  double norm_est = norm_estimate_cplx(apply_H, apply_H, block.size());
  double ratio = norm_est / norm_exact;

  // Log("norm_exact: {}", norm_exact);
  // Log("norm_est: {}", norm_est);
  if (norm_exact > 1e-12) {
    // Log("ratio: {}", ratio);
    REQUIRE(((ratio <= 1.00001) && (ratio > 0.3)));
  }
}

TEST_CASE("norm_estimate", "[algorithms]") {

  using namespace hydra::testcases::spinhalf;
  using hydra::testcases::electron::get_cyclic_group_irreps_mult;

  using namespace hydra;
  using namespace arma;

  Log.set_verbosity(2);

  for (int n = 50; n <= 350; n += 50) {
    Log("norm_estimate for random real non-hermitian matrices, D={}", n);

    // Random real generic matrices
    for (int seed = 1; seed <= 10; ++seed) {
      arma_rng::set_seed(seed);
      mat A(n, n, fill::randn);
      double norm_exact = norm(A, 1);
      auto apply_A = [&A](arma::vec const &v) { return vec(A * v); };
      auto apply_A_T = [&A](arma::vec const &v) { return vec(A.t() * v); };
      double norm_est = norm_estimate_real(apply_A, apply_A_T, n);
      double ratio = norm_est / norm_exact;
      // Log("ratio: {}", ratio);
      // Log("norm_exact: {}", norm_exact);
      // Log("norm_est: {}", norm_est);
      REQUIRE(((ratio < 1.0001) && (ratio > 0.3)));
    }

    Log("norm_estimate for random cplx non-hermitian matrices, D={}", n);
    for (int seed = 1; seed <= 3; ++seed) {
      arma_rng::set_seed(seed);
      cx_mat A(n, n, fill::randn);
      double norm_exact = norm(A, 1);
      auto apply_A = [&A](arma::cx_vec const &v) { return cx_vec(A * v); };
      auto apply_A_T = [&A](arma::cx_vec const &v) {
        return cx_vec(A.t() * v);
      };
      double norm_est = norm_estimate_cplx(apply_A, apply_A_T, n);
      double ratio = norm_est / norm_exact;
      // Log("ratio: {}", ratio);
      // Log("norm_exact: {}", norm_exact);
      // Log("norm_est: {}", norm_est);
      REQUIRE(((ratio <= 1.00001) && (ratio > 0.3)));
    }
  }
  for (int n_sites = 2; n_sites < 12; ++n_sites) {
    HydraPrint(n_sites);

    Log("norm_estimate for Heisenberg all-to-all random, N={}", n_sites);

    // Random HB alltoall
    auto bonds = hydra::testcases::spinhalf::HB_alltoall(n_sites);
    auto block = Spinhalf(n_sites);
    test_operator_norm_real(block, bonds);
    for (int nup = 0; nup <= n_sites; ++nup) {
      auto block = Spinhalf(n_sites, nup);
      test_operator_norm_real(block, bonds);
    }

    Log("norm_estimate for Heisenberg chain symmetric, N={}", n_sites);

    // HB chain with lattice symmetries
    auto [group, irreps, multiplicities] =
        get_cyclic_group_irreps_mult(n_sites);
    auto bondlist = HBchain(n_sites, 3.21, 0.123);
    for (int nup = 0; nup <= n_sites; ++nup) {
      auto block = Spinhalf(n_sites, nup);
      test_operator_norm_real(block, bonds);
      for (auto irrep : irreps) {
        // HydraPrint(irrep);
        auto block = Spinhalf(n_sites, nup, group, irrep);
        if (hydra::is_real(irrep)) {
          test_operator_norm_real(block, bonds);
        } else {
          test_operator_norm_cplx(block, bonds);
        }
      }
    }
  }

  {
    Log("norm_estimate for tj_symmetric_matrix: tJ 3x3 triangular s");
    std::string lfile = "data/triangular.9.hop.sublattices.tsl.lat";
    int n_sites = 9;
    auto bonds = read_bondlist(lfile);
    bonds["T"] = 1.0;
    bonds["J"] = 0.4;
    auto permutations = hydra::read_permutations(lfile);
    auto group = PermutationGroup(permutations);

    std::vector<std::pair<std::string, int>> rep_name_mult = {
        {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
        {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
        {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
        {"Y.C1.A", 6}};

    std::vector<Representation> irreps;
    std::vector<int> multiplicities;
    for (int nup = 0; nup <= n_sites; ++nup) {
      for (int ndn = 0; ndn <= n_sites; ++ndn) {

        if (nup + ndn > n_sites)
          continue;

        auto block = tJ(n_sites, nup, ndn);
        test_operator_norm_real(block, bonds);

        for (auto [name, mult] : rep_name_mult) {
          auto irrep = read_representation(lfile, name);
          auto block = tJ(n_sites, nup, ndn, group, irrep);
          if (hydra::is_real(irrep)) {
            test_operator_norm_real(block, bonds);
          } else {
            test_operator_norm_cplx(block, bonds);
          }
        }
      }
    }
  } // 9 site tJ model triangular

  {
    int n_sites = 9;

    // test a 3x3 triangular lattice with complex flux
    Log("norm_estimate for tj_symmetric_matrix: tJ 3x3 triangular staggered "
        "flux, complex");
    std::string lfile =
        "data/triangular.9.tup.phi.tdn.nphi.sublattices.tsl.lat";

    auto bonds = read_bondlist(lfile);
    std::vector<double> etas{0.0, 0.1, 0.2, 0.3};
    auto permutations = hydra::read_permutations(lfile);
    auto group = PermutationGroup(permutations);

    std::vector<std::pair<std::string, int>> rep_name_mult = {
        {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
        {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
        {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
        {"Y.C1.A", 6}};

    std::vector<Representation> irreps;
    std::vector<int> multiplicities;

    for (auto eta : etas) {
      bonds["TPHI"] = complex(cos(eta * M_PI), sin(eta * M_PI));
      bonds["JPHI"] = complex(cos(2 * eta * M_PI), sin(2 * eta * M_PI));
      for (int nup = 0; nup <= n_sites; ++nup) {
        for (int ndn = 0; ndn <= n_sites; ++ndn) {

          if (nup + ndn > n_sites)
            continue;

          auto block = tJ(n_sites, nup, ndn);
          test_operator_norm_cplx(block, bonds);

          for (auto [name, mult] : rep_name_mult) {
            auto irrep = read_representation(lfile, name);
            auto block = tJ(n_sites, nup, ndn, group, irrep);
            test_operator_norm_cplx(block, bonds);
          }
        }
      }
    }
  }
}
