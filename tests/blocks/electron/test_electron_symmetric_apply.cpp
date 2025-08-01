// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <iostream>

#include "../electron/testcases_electron.hpp"
#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/operators/logic/real.hpp>

using namespace xdiag;

void test_electron_symmetric_apply(OpSum ops, int64_t nsites,
                                   std::vector<Representation> irreps) {

  for (int64_t nup = 0; nup <= nsites; ++nup) {
    for (int64_t ndn = 0; ndn <= nsites; ++ndn) {

      for (auto irrep : irreps) {
        // Log("nup: {}, ndn: {}", nup, ndn);
        // XDIAG_SHOW(irrep);

        // Create block and matrix for comparison
        // tic();
        auto block = Electron(nsites, nup, ndn, irrep);
        // Log("block.isreal {}", isreal(block));
        // toc("create block");

        if (block.size() > 0) {
          // tic();
          auto H = matrixC(ops, block, block);
          // toc("create matrix");

          // tic();
          // // REQUIRE(arma::norm(H - H.t()) < 1e-12);
          // toc("transpose");

          // REQUIRE(H.is_hermitian(1e-8));
          // Check whether apply gives the same as matrix multiplication
          // tic();
          arma::cx_vec v(block.size(), arma::fill::randn);
          arma::cx_vec w1 = H * v;
          arma::cx_vec w2(block.size(), arma::fill::randn);
          apply(ops, block, v, block, w2);
          REQUIRE(isapprox(w1, w2));
          arma::cx_mat m(block.size(), 5, arma::fill::randn);
          arma::cx_mat n1 = H * m;
          arma::cx_mat n2(block.size(), 5, arma::fill::zeros);
          apply(ops, block, m, block, n2);
          REQUIRE(isapprox(n1, n2));
          // toc("apply");

          // Compute eigenvalues and compare
          // tic();
          arma::vec evals_mat;
          arma::eig_sym(evals_mat, H);
          // toc("evals full");
          // tic();
          double e0_mat = evals_mat(0);
          double e0_app = eigval0(ops, block);
          // Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
          REQUIRE(std::abs(e0_mat - e0_app) < 1e-7);
          // toc("evals lcs");

          // Compute eigenvalues with real arithmitic
          // tic();
          if (isreal(block) && isreal(ops)) {
            auto H_real = matrix(ops, block, block);
            arma::vec evals_mat_real;
            arma::eig_sym(evals_mat_real, H_real);

            REQUIRE(isapprox(evals_mat_real, evals_mat));

            double e0_mat_real = evals_mat_real(0);
            double e0_app_real = eigval0(ops, block);
            REQUIRE(std::abs(e0_mat_real - e0_app_real) < 1e-7);
          }
          // toc("real");
        }
      }
    }
  }
}

void test_hubbard_symmetric_apply_chains(int64_t nsites) {
  using namespace xdiag::testcases::electron;

  // Without Heisenberg term
  Log.out("electron_symmetric_apply: Hubbard chain, nsites: {}", nsites);
  auto ops = get_linear_chain(nsites, 1.0, 5.0);
  auto irreps = get_cyclic_group_irreps(nsites);
  test_electron_symmetric_apply(ops, nsites, irreps);

  // With Heisenberg term
  Log.out("electron_symmetric_apply: Hubbard chain, nsites: {} (+ "
          "Heisenberg terms)",
          nsites);
  auto ops_hb = get_linear_chain_hb(nsites, 0.4);
  test_electron_symmetric_apply(ops_hb, nsites, irreps);
}

TEST_CASE("electron_symmetric_apply", "[electron]") try {

  // Test linear chains
  for (int64_t nsites = 2; nsites < 7; ++nsites) {
    test_hubbard_symmetric_apply_chains(nsites);
  }

  // test a 3x3 triangular lattice
  Log.out("electron_symmetric_apply: Hubbard 3x3 triangular");
  std::string lfile =
      XDIAG_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.toml";

  auto fl = FileToml(lfile);
  auto ops = fl["Interactions"].as<OpSum>();
  ops += "U" * Op("HubbardU");
  ops["T"] = 1.0;
  ops["U"] = 5.0;

  std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
      {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
      {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
      {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
      {"Y.C1.A", 6}};
  std::vector<Representation> irreps;
  for (auto [name, mult] : rep_name_mult) {
    irreps.push_back(read_representation(fl, name));
    (void)mult;
  }
  test_electron_symmetric_apply(ops, 9, irreps);

  // test a 3x3 triangular lattice with Heisenberg terms
  Log.out(
      "electron_symmetric_apply: Hubbard 3x3 triangular (+ Heisenberg terms)");
  auto ops_hb = ops;
  for (auto [cpl, op] : ops) {
    if (op.type() == "Hop") {
      ops_hb += "J" * Op("SdotS", {op[0], op[1]});
    }
  }
  ops_hb["J"] = 0.4;
  test_electron_symmetric_apply(ops_hb, 9, irreps);

  // test a 3x3 triangular lattice with complex hoppings
  {
    Log.out("electron_symmetric_apply: Hubbard 3x3 triangular (complex)");
    std::string lfile = XDIAG_DIRECTORY
        "/misc/data/triangular.9.tup.phi.tdn.nphi.sublattices.tsl.toml";
    auto fl = FileToml(lfile);
    auto ops = fl["Interactions"].as<OpSum>();
    ops += "U" * Op("HubbardU");
    ops["TPHI"] = complex(0.5, 0.5);
    ops["JPHI"] = 0.;
    ops["U"] = 5.0;

    std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
        {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
        {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
        {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
        {"Y.C1.A", 6}};
    irreps.clear();
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(read_representation(fl, name));
      (void)mult;
    }
    test_electron_symmetric_apply(ops, 9, irreps);
  }
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
