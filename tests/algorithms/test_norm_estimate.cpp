#include "../catch.hpp"

#include <iostream>

#include "../blocks/electron/testcases_electron.hpp"
#include "../blocks/spinhalf/testcases_spinhalf.hpp"
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/norm_estimate.hpp>
#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

template <typename block_t>
void test_operator_norm_real(block_t const &block, OpSum const &ops) {
  using namespace xdiag;
  using namespace arma;

  auto H = matrix(ops, block, block);
  double norm_exact = norm(H, 1);
  auto apply_H = [&H](vec const &v) { return vec(H * v); };
  double norm_est = norm_estimate(ops, block);
  double ratio = norm_est / norm_exact;

  // Log("norm_exact: {}", norm_exact);
  // Log("norm_est: {}", norm_est);
  if (norm_exact > 1e-12) {
    // Log("ratio: {}", ratio);
    REQUIRE(((ratio <= 1.00001) && (ratio > 0.2)));
  }
}

template <typename block_t>
void test_operator_norm_cplx(block_t const &block, OpSum const &ops) {
  using namespace xdiag;
  using namespace arma;

  auto H = matrixC(ops, block, block);
  double norm_exact = norm(H, 1);
  double norm_est = norm_estimate(ops, block);
  double ratio = norm_est / norm_exact;

  // Log("norm_exact: {}", norm_exact);
  // Log("norm_est: {}", norm_est);
  if (norm_exact > 1e-12) {
    // Log("ratio: {}", ratio);
    REQUIRE(((ratio <= 1.00001) && (ratio > 0.3)));
  }
}

TEST_CASE("norm_estimate", "[algorithms]") {

  using namespace xdiag::testcases::spinhalf;
  using xdiag::testcases::electron::get_cyclic_group_irreps_mult;

  using namespace xdiag;
  using namespace arma;

  Log.set_verbosity(1);

  for (int n = 50; n <= 350; n += 50) {
    Log("norm_estimate for random real non-hermitian matrices, D={}", n);

    // Random real generic matrices
    for (int seed = 1; seed <= 10; ++seed) {
      arma_rng::set_seed(seed);
      mat A(n, n, fill::randn);
      double norm_exact = norm(A, 1);
      double norm_est = norm_estimate(A);
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
      double norm_est = norm_estimate(A);
      double ratio = norm_est / norm_exact;
      // Log("ratio: {}", ratio);
      // Log("norm_exact: {}", norm_exact);
      // Log("norm_est: {}", norm_est);
      REQUIRE(((ratio <= 1.00001) && (ratio > 0.3)));
    }
  }
  for (int n_sites = 3; n_sites < 12; ++n_sites) {
    // XDIAG_SHOW(n_sites);
    {
      Log("norm_estimate for Heisenberg all-to-all random, N={}", n_sites);

      // Random HB alltoall
      auto ops = xdiag::testcases::spinhalf::HB_alltoall(n_sites);
      auto block = Spinhalf(n_sites);
      test_operator_norm_real(block, ops);
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto block = Spinhalf(n_sites, nup);
        test_operator_norm_real(block, ops);
      }
    }

    Log("norm_estimate for Heisenberg chain symmetric, N={}", n_sites);

    // HB chain with lattice symmetries
    auto [group, irreps, multiplicities] =
        get_cyclic_group_irreps_mult(n_sites);
    (void)multiplicities;
    auto ops = HBchain(n_sites, 3.21, 0.123);
    for (int nup = 0; nup <= n_sites; ++nup) {
      auto block = Spinhalf(n_sites, nup);
      test_operator_norm_real(block, ops);
      for (auto irrep : irreps) {
        // XDIAG_SHOW(irrep);
        auto block = Spinhalf(n_sites, nup, group, irrep);
        if (irrep.isreal()) {
          test_operator_norm_real(block, ops);
        } else {
          test_operator_norm_cplx(block, ops);
        }
      }
    }
  }

  {
    Log("norm_estimate for tj_symmetric_matrix: tJ 3x3 triangular s");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.lat";
    int n_sites = 9;
    auto ops = read_opsum(lfile);
    ops["T"] = 1.0;
    ops["J"] = 0.4;
    auto permutations = xdiag::read_permutations(lfile);
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
        test_operator_norm_real(block, ops);

        for (auto [name, mult] : rep_name_mult) {
          (void)mult;
          auto irrep = read_representation(lfile, name);
          auto block = tJ(n_sites, nup, ndn, group, irrep);
          if (irrep.isreal()) {
            test_operator_norm_real(block, ops);
          } else {
            test_operator_norm_cplx(block, ops);
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
    std::string lfile = XDIAG_DIRECTORY
        "/misc/data/triangular.9.tup.phi.tdn.nphi.sublattices.tsl.lat";

    auto ops = read_opsum(lfile);
    std::vector<double> etas{0.0, 0.1, 0.2, 0.3};
    auto permutations = xdiag::read_permutations(lfile);
    auto group = PermutationGroup(permutations);

    std::vector<std::pair<std::string, int>> rep_name_mult = {
        {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
        {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
        {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
        {"Y.C1.A", 6}};

    std::vector<Representation> irreps;
    std::vector<int> multiplicities;

    for (auto eta : etas) {
      ops["TPHI"] = complex(cos(eta * M_PI), sin(eta * M_PI));
      ops["JPHI"] = complex(cos(2 * eta * M_PI), sin(2 * eta * M_PI));
      for (int nup = 0; nup <= n_sites; ++nup) {
        for (int ndn = 0; ndn <= n_sites; ++ndn) {

          if (nup + ndn > n_sites)
            continue;

          auto block = tJ(n_sites, nup, ndn);
          test_operator_norm_cplx(block, ops);

          for (auto [name, mult] : rep_name_mult) {
            (void)mult;
            auto irrep = read_representation(lfile, name);
            auto block = tJ(n_sites, nup, ndn, group, irrep);
            test_operator_norm_cplx(block, ops);
          }
        }
      }
    }
  }
}
