// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <iostream>

#include "../blocks/electron/testcases_electron.hpp"
#include "../blocks/spinhalf/testcases_spinhalf.hpp"
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algebra/sparse/csr_matrix.hpp>
#include <xdiag/algorithms/norm_estimate.hpp>
#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/io/read.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/xdiag_show.hpp>

using namespace xdiag;

template <typename block_t, typename op_t>
void test_operator_norm_real(block_t const &block, op_t const &ops) try {
  using namespace xdiag;
  using namespace arma;

  mat H;
  if constexpr (std::is_same<op_t, OpSum>::value) {
    H = matrix(ops, block);
  } else {
    H = to_dense(ops);
  }
  double norm_exact = norm(H, 1);
  double norm_est = norm_estimate(ops, block);

  if (norm_exact > 1e-12) {
    double ratio = norm_est / norm_exact;
    // Log("norm_exact: {}", norm_exact);
    // Log("norm_est: {}", norm_est);
    // Log("ratio: {}", ratio);
    REQUIRE(((ratio <= 1.00001) && (ratio > 0.2)));
  } else {
    REQUIRE(norm_est < 1e-12);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename block_t, typename op_t>
void test_operator_norm_cplx(block_t const &block, op_t const &ops) {
  using namespace xdiag;
  using namespace arma;

  cx_mat H;
  if constexpr (std::is_same<op_t, OpSum>::value) {
    H = matrixC(ops, block);
  } else {
    H = to_dense(ops);
  }
  double norm_exact = norm(H, 1);
  double norm_est = norm_estimate(ops, block);

  if (norm_exact > 1e-12) {
    double ratio = norm_est / norm_exact;
    // Log("norm_exact: {}", norm_exact);
    // Log("norm_est: {}", norm_est);
    // Log("ratio: {}", ratio);
    REQUIRE(((ratio <= 1.00001) && (ratio > 0.3)));
  } else {
    REQUIRE(norm_est < 1e-12);
  }
}

TEST_CASE("norm_estimate", "[algorithms]") try {

  using namespace xdiag::testcases::spinhalf;
  using xdiag::testcases::electron::get_cyclic_group_irreps_mult;

  using namespace xdiag;
  using namespace arma;

  Log.set_verbosity(0);

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
  for (int nsites = 3; nsites < 8; ++nsites) {
    // XDIAG_SHOW(nsites);
    {
      Log("norm_estimate for Heisenberg all-to-all random, N={}", nsites);

      // Random HB alltoall
      auto ops = xdiag::testcases::spinhalf::HB_alltoall(nsites);
      auto block = Spinhalf(nsites);
      auto csr = csr_matrix(ops, block);
      test_operator_norm_real(block, ops);
      test_operator_norm_real(block, csr);
      for (int nup = 0; nup <= nsites; ++nup) {
        auto block = Spinhalf(nsites, nup);
        auto csr = csr_matrix(ops, block);
        test_operator_norm_real(block, ops);
        test_operator_norm_real(block, csr);
      }
    }

    Log("norm_estimate for Heisenberg chain symmetric, N={}", nsites);

    // HB chain with lattice symmetries
    auto [irreps, multiplicities] = get_cyclic_group_irreps_mult(nsites);
    (void)multiplicities;
    auto ops = HBchain(nsites, 3.21, 0.123);
    for (int nup = 0; nup <= nsites; ++nup) {
      auto block = Spinhalf(nsites, nup);
      auto csr = csr_matrix(ops, block);
      test_operator_norm_real(block, ops);
      test_operator_norm_real(block, csr);

      for (auto irrep : irreps) {
        auto block = Spinhalf(nsites, nup, irrep);
        if (irrep.isreal()) {
          auto csr = csr_matrix(ops, block);
          test_operator_norm_real(block, ops);
          test_operator_norm_real(block, csr);
        } else {
          auto csr = csr_matrixC(ops, block);
          test_operator_norm_cplx(block, ops);
          test_operator_norm_cplx(block, csr);
        }
      }
    }
  }

  // {
  //   Log("norm_estimate for tj_symmetric_matrix: tJ 3x3 triangular s");
  //   std::string lfile =
  //       XDIAG_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.toml";
  //   int nsites = 9;
  //   auto fl = FileToml(lfile);
  //   auto ops = fl["Interactions"].as<OpSum>();
  //   ops["T"] = 1.0;
  //   ops["J"] = 0.4;

  //   std::vector<std::pair<std::string, int>> rep_name_mult = {
  //       {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
  //       {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
  //       {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
  //       {"Y.C1.A", 6}};

  //   std::vector<Representation> irreps;
  //   std::vector<int> multiplicities;
  //   for (int nup = 0; nup <= nsites; ++nup) {
  //     for (int ndn = 0; ndn <= nsites; ++ndn) {

  //       if (nup + ndn > nsites)
  //         continue;

  //       auto block = tJ(nsites, nup, ndn);
  //       auto csr = csr_matrix(ops, block);
  //       test_operator_norm_real(block, ops);
  //       test_operator_norm_real(block, csr);

  //       for (auto [name, mult] : rep_name_mult) {
  //         (void)mult;
  //         auto irrep = read_representation(fl, name);
  //         auto block = tJ(nsites, nup, ndn, irrep);
  //         if (irrep.isreal()) {
  //           auto csr = csr_matrix(ops, block);
  //           test_operator_norm_real(block, ops);
  //           test_operator_norm_real(block, csr);
  //         } else {
  //           auto csr = csr_matrixC(ops, block);
  //           test_operator_norm_cplx(block, ops);
  //           test_operator_norm_cplx(block, csr);
  //         }
  //       }
  //     }
  //   }
  // } // 9 site tJ model triangular

  // {
  //   int nsites = 9;

  //   // test a 3x3 triangular lattice with complex flux
  //   Log("norm_estimate for tj_symmetric_matrix: tJ 3x3 triangular staggered "
  //       "flux, complex");
  //   std::string lfile = XDIAG_DIRECTORY
  //       "/misc/data/triangular.9.tup.phi.tdn.nphi.sublattices.tsl.toml";

  //   auto fl = FileToml(lfile);
  //   auto ops = fl["Interactions"].as<OpSum>();
  //   std::vector<double> etas{0.0, 0.1, 0.2, 0.3};

  //   std::vector<std::pair<std::string, int>> rep_name_mult = {
  //       {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
  //       {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
  //       {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
  //       {"Y.C1.A", 6}};

  //   std::vector<Representation> irreps;
  //   std::vector<int> multiplicities;

  //   for (auto eta : etas) {
  //     ops["TPHI"] = complex(cos(eta * M_PI), sin(eta * M_PI));
  //     ops["JPHI"] = complex(cos(2 * eta * M_PI), sin(2 * eta * M_PI));
  //     for (int nup = 0; nup <= nsites; ++nup) {
  //       for (int ndn = 0; ndn <= nsites; ++ndn) {

  //         if (nup + ndn > nsites)
  //           continue;

  //         auto block = tJ(nsites, nup, ndn);
  //         auto csr = csr_matrixC(ops, block);
  //         test_operator_norm_cplx(block, ops);
  //         test_operator_norm_cplx(block, csr);

  //         for (auto [name, mult] : rep_name_mult) {
  //           (void)mult;
  //           auto irrep = read_representation(fl, name);
  //           auto block = tJ(nsites, nup, ndn, irrep);
  //           auto csr = csr_matrixC(ops, block);
  //           test_operator_norm_cplx(block, ops);
  //           test_operator_norm_cplx(block, ops);
  //         }
  //       }
  //     }
  //   }
  // }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
