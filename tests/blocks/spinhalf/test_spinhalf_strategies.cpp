// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "test_spinhalf_strategies.hpp"

#include "../../catch.hpp"

#include <cassert>

#include <xdiag/armadillo.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/config.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/math/isapprox.hpp>
#include <xdiag/matrices/apply.hpp>
#include <xdiag/matrices/matrix.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

void test_apply(int64_t N, OpSum ops) {
  for (int64_t nup = 1; nup <= N; ++nup) {
    auto block = Spinhalf(N, nup);
    auto H = matrix(ops, block, block);
    REQUIRE(H.is_hermitian(1e-8));

    arma::vec v(block.size(), arma::fill::randn);
    arma::vec w2(block.size(), arma::fill::zeros);
    apply(ops, block, v, block, w2);

    arma::vec w3 = H * H * v;
    arma::vec w4(block.size(), arma::fill::zeros);
    apply(ops, block, w2, block, w4);
    REQUIRE(isapprox(w3, w4));

    arma::vec evals_mat;
    arma::eig_sym(evals_mat, H);
    double e0_mat = evals_mat(0);
    auto [e0_app, ev] = eig0(ops, block);
    REQUIRE(std::abs(e0_mat - e0_app) < 1e-12);
  }
}

void test_apply_mat(int64_t N, OpSum ops) {
  for (int64_t nup = 1; nup <= N; ++nup) {
    auto block = Spinhalf(N, nup);
    auto H = matrix(ops, block, block);
    REQUIRE(H.is_hermitian(1e-8));

    arma::mat v(block.size(), 5, arma::fill::randn);
    arma::mat w2(block.size(), 5, arma::fill::zeros);
    apply(ops, block, v, block, w2);

    arma::mat w3 = H * H * v;
    arma::mat w4(block.size(), 5, arma::fill::zeros);
    apply(ops, block, w2, block, w4);
    REQUIRE(isapprox(w3, w4));

    arma::vec evals_mat;
    arma::eig_sym(evals_mat, H);
    double e0_mat = evals_mat(0);
    auto [e0_app, ev] = eig0(ops, block);
    REQUIRE(std::abs(e0_mat - e0_app) < 1e-12);
  }
}

void test_onsite(std::string op1, std::string op12) {
  for (int64_t nsites = 2; nsites < 5; ++nsites) {
    for (int64_t nup = 0; nup <= nsites; ++nup) {
      auto b = Spinhalf(nsites, nup);
      for (int64_t s = 0; s < nsites; ++s) {
        arma::mat m1 = matrix(Op(op1, s), b);
        arma::mat m12 = matrix(Op(op12, {s, s}), b);
        arma::mat m1sq = m1 * m1;
        REQUIRE(isapprox(m12, m1sq));
      }
    }
  }
}

void test_kitaev_gamma(double K, double G,
                       std::vector<std::pair<std::string, double>> irrep_names_e0) {
  std::string lfile =
      XDIAG_DIRECTORY "/misc/data/kitaev_gamma/lattice-files/"
                      "honeycomb.8.HeisenbergKitaevGamma.fsl.toml";

  auto fl = FileToml(lfile);
  auto ops_read = read_opsum(fl, "Interactions");

  arma::cx_mat sx(arma::mat({{0., 0.5}, {0.5, 0.}}),
                  arma::mat({{0., 0.}, {0., 0.}}));
  arma::cx_mat sy(arma::mat({{0., 0.}, {0., 0.}}),
                  arma::mat({{0., -0.5}, {0.5, 0.}}));
  arma::cx_mat sz(arma::mat({{0.5, 0.0}, {0.0, -0.5}}),
                  arma::mat({{0., 0.}, {0., 0.0}}));

  arma::cx_mat sxsx = arma::kron(sx, sx);
  arma::cx_mat sysy = arma::kron(sy, sy);
  arma::cx_mat szsz = arma::kron(sz, sz);
  arma::cx_mat gsx = arma::kron(sy, sz) + arma::kron(sz, sy);
  arma::cx_mat gsy = arma::kron(sx, sz) + arma::kron(sz, sx);
  arma::cx_mat gsz = arma::kron(sx, sy) + arma::kron(sy, sx);

  auto ops = OpSum();
  for (auto term : ops_read) {
    Op const &op = term.monomial[0];
    std::string type = op.type();
    auto sites = op.sites();
    if (type == "KITAEVX") {
      ops += K * Op("Matrix", sites, sxsx);
    } else if (type == "KITAEVY") {
      ops += K * Op("Matrix", sites, sysy);
    } else if (type == "KITAEVZ") {
      ops += K * Op("Matrix", sites, szsz);
    } else if (type == "GAMMAX") {
      ops += G * Op("Matrix", sites, gsx);
    } else if (type == "GAMMAY") {
      ops += G * Op("Matrix", sites, gsy);
    } else if (type == "GAMMAZ") {
      ops += G * Op("Matrix", sites, gsz);
    }
  }

  for (auto [name, e0_reference] : irrep_names_e0) {
    auto irrep = read_representation(fl, name);
    auto block = Spinhalf(8, irrep);
    arma::cx_mat H = matrixC(ops, block);
    arma::vec eigs;
    arma::eig_sym(eigs, H);
    double e0_computed = eigs(0);
    REQUIRE(isapprox(e0_reference, e0_computed));
  }
}

void test_spinhalf_symmetric_apply(OpSum ops, int64_t nsites,
                                   std::vector<Representation> const &irreps) {
  for (int64_t nup = 0; nup <= nsites; ++nup) {
    for (auto irrep : irreps) {
      auto block = Spinhalf(nsites, nup, irrep);

      if (block.size() > 0) {
        auto H = matrixC(ops, block, block);
        REQUIRE(arma::norm(H - H.t()) < 1e-12);

        arma::cx_vec v(block.size(), arma::fill::randn);
        arma::cx_vec w1 = H * v;
        arma::cx_vec w2(block.size(), arma::fill::zeros);
        apply(ops, block, v, block, w2);
        REQUIRE(isapprox(w1, w2));

        arma::cx_mat m(block.size(), 5, arma::fill::randn);
        arma::cx_mat n1 = H * m;
        arma::cx_mat n2(block.size(), 5, arma::fill::zeros);
        apply(ops, block, m, block, n2);
        REQUIRE(isapprox(n1, n2));

        arma::vec evals_mat;
        arma::eig_sym(evals_mat, H);
        double e0_mat = evals_mat(0);
        double e0_app = eigval0(ops, block);
        REQUIRE(std::abs(e0_mat - e0_app) < 1e-7);

        if (isreal(block) && isreal(ops)) {
          auto H_real = matrix(ops, block, block);
          arma::vec evals_mat_real;
          arma::eig_sym(evals_mat_real, H_real);
          REQUIRE(isapprox(evals_mat_real, evals_mat));

          double e0_mat_real = evals_mat_real(0);
          double e0_app_real = eigval0(ops, block);
          REQUIRE(std::abs(e0_mat_real - e0_app_real) < 1e-7);
        }
      }
    }
  }
}

void test_spinhalf_symmetric_apply_no_sz(
    OpSum ops, int64_t nsites,
    std::vector<Representation> const &irreps) {
  for (auto irrep : irreps) {
    auto block = Spinhalf(nsites, irrep);

    if (block.size() > 0) {
      auto H = matrixC(ops, block, block);
      REQUIRE(arma::norm(H - H.t()) < 1e-12);

      arma::cx_vec v(block.size(), arma::fill::randn);
      arma::cx_vec w1 = H * v;
      arma::cx_vec w2(block.size(), arma::fill::zeros);
      apply(ops, block, v, block, w2);
      REQUIRE(isapprox(w1, w2));

      arma::cx_mat m(block.size(), 5, arma::fill::randn);
      arma::cx_mat n1 = H * m;
      arma::cx_mat n2(block.size(), 5, arma::fill::zeros);
      apply(ops, block, m, block, n2);
      REQUIRE(isapprox(n1, n2));

      arma::vec evals_mat;
      arma::eig_sym(evals_mat, H);
      double e0_mat = evals_mat(0);
      double e0_app = eigval0(ops, block);
      REQUIRE(std::abs(e0_mat - e0_app) < 1e-7);

      if (isreal(block) && isreal(ops)) {
        auto H_real = matrix(ops, block, block);
        arma::vec evals_mat_real;
        arma::eig_sym(evals_mat_real, H_real);
        REQUIRE(isapprox(evals_mat_real, evals_mat));

        double e0_mat_real = evals_mat_real(0);
        double e0_app_real = eigval0(ops, block);
        REQUIRE(std::abs(e0_mat_real - e0_app_real) < 1e-7);
      }
    }
  }
}

void test_spinhalf_symmetric_spectra(OpSum ops, int64_t nsites,
                                     std::vector<Representation> irreps,
                                     std::vector<int64_t> multiplicities) {
  assert(irreps.size() == multiplicities.size());

  for (int64_t nup = 0; nup <= nsites; ++nup) {
    auto spinhalf_nosym = Spinhalf(nsites, nup);

    if (spinhalf_nosym.size() < 1000) {
      std::vector<double> eigs_sym;

      auto H_nosym = matrixC(ops, spinhalf_nosym, spinhalf_nosym);
      REQUIRE(arma::norm(H_nosym - H_nosym.t()) < 1e-8);
      arma::vec eigs_nosym;
      arma::eig_sym(eigs_nosym, H_nosym);

      for (int64_t k = 0; k < (int64_t)irreps.size(); ++k) {
        auto irrep = irreps[k];
        int64_t multiplicity = multiplicities[k];
        auto spinhalf = Spinhalf(nsites, nup, irrep);

        if (spinhalf.size() > 0) {
          auto H_sym = matrixC(ops, spinhalf, spinhalf);
          REQUIRE(arma::norm(H_sym - H_sym.t()) < 1e-12);
          arma::vec eigs_sym_k;
          arma::eig_sym(eigs_sym_k, H_sym);

          if (isreal(spinhalf) && isreal(ops)) {
            auto H_sym_real = matrix(ops, spinhalf, spinhalf);
            arma::vec eigs_sym_k_real;
            arma::eig_sym(eigs_sym_k_real, H_sym_real);
            REQUIRE(isapprox(eigs_sym_k, eigs_sym_k_real));
          }

          for (auto eig : eigs_sym_k) {
            for (int64_t i = 0; i < multiplicity; ++i) {
              eigs_sym.push_back(eig);
            }
          }
        }
      }
      std::sort(eigs_sym.begin(), eigs_sym.end());
      REQUIRE(isapprox(arma::vec(eigs_sym), eigs_nosym));
    }
  }
}

void test_spinhalf_symmetric_spectra_no_sz(
    OpSum ops, int64_t nsites, std::vector<Representation> irreps,
    std::vector<int64_t> multiplicities) {
  assert(irreps.size() == multiplicities.size());

  auto spinhalf_nosym = Spinhalf(nsites);

  if (spinhalf_nosym.size() < 1000) {
    std::vector<double> eigs_sym;

    auto H_nosym = matrixC(ops, spinhalf_nosym, spinhalf_nosym);
    REQUIRE(H_nosym.is_hermitian());
    arma::vec eigs_nosym;
    arma::eig_sym(eigs_nosym, H_nosym);

    for (int64_t k = 0; k < (int64_t)irreps.size(); ++k) {
      auto irrep = irreps[k];
      int64_t multiplicity = multiplicities[k];
      auto spinhalf = Spinhalf(nsites, irrep);

      if (spinhalf.size() > 0) {
        auto H_sym = matrixC(ops, spinhalf, spinhalf);
        REQUIRE(arma::norm(H_sym - H_sym.t()) < 1e-12);
        arma::vec eigs_sym_k;
        arma::eig_sym(eigs_sym_k, H_sym);

        // cross-check: no-Sz irrep spectrum equals union of Sz-sector spectra
        auto eigs_sym_k_sz = std::vector<double>();
        for (int64_t nup = 0; nup <= nsites; ++nup) {
          auto spinhalf_sz = Spinhalf(nsites, nup, irrep);
          auto H_sym_sz = matrixC(ops, spinhalf_sz, spinhalf_sz);
          arma::vec es;
          arma::eig_sym(es, H_sym_sz);
          for (auto e : es) {
            eigs_sym_k_sz.push_back(e);
          }
        }
        std::sort(eigs_sym_k_sz.begin(), eigs_sym_k_sz.end());
        REQUIRE(isapprox(eigs_sym_k, arma::vec(eigs_sym_k_sz)));

        if (isreal(spinhalf) && isreal(ops)) {
          auto H_sym_real = matrix(ops, spinhalf, spinhalf);
          arma::vec eigs_sym_k_real;
          arma::eig_sym(eigs_sym_k_real, H_sym_real);
          REQUIRE(isapprox(eigs_sym_k, eigs_sym_k_real));
        }

        for (auto eig : eigs_sym_k) {
          for (int64_t i = 0; i < multiplicity; ++i) {
            eigs_sym.push_back(eig);
          }
        }
      }
    }
    std::sort(eigs_sym.begin(), eigs_sym.end());
    REQUIRE(isapprox(arma::vec(eigs_sym), eigs_nosym));
  }
}
