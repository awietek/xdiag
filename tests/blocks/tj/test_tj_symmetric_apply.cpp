#include "../../catch.hpp"

#include <iostream>

#include "../electron/testcases_electron.hpp"
#include "../tj/testcases_tj.hpp"
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/utils/close.hpp>

using namespace xdiag;

void test_apply_tj_symmetric(OpSum ops, PermutationGroup space_group,
                             std::vector<Representation> irreps) {
  int64_t n_sites = space_group.n_sites();

  for (int64_t nup = 0; nup <= n_sites; ++nup) {
    for (int64_t ndn = 0; ndn <= n_sites; ++ndn) {

      if (nup + ndn > n_sites)
        continue;

      for (int64_t k = 0; k < (int64_t)irreps.size(); ++k) {
        auto irrep = irreps[k];
        auto block = tJ(n_sites, nup, ndn, space_group, irrep);

        if (block.size() > 0) {
          auto H_sym = matrixC(ops, block, block);
          REQUIRE(arma::norm(H_sym - H_sym.t()) < 1e-12);

          // Check whether apply gives the same as matrix multiplication
          arma::cx_vec v(block.size(), arma::fill::randn);
          arma::cx_vec w1 = H_sym * v;
          arma::cx_vec w2(block.size(), arma::fill::zeros);
          apply(ops, block, v, block, w2);
          REQUIRE(close(w1, w2));
          arma::cx_mat m(block.size(), 5, arma::fill::randn);
          arma::cx_mat n1 = H_sym * m;
          arma::cx_mat n2(block.size(), 5, arma::fill::zeros);
          apply(ops, block, m, block, n2);
          REQUIRE(close(n1, n2));

          // Compute eigenvalues and compare
          arma::vec evals_mat;
          arma::eig_sym(evals_mat, H_sym);
          double e0_mat = evals_mat(0);
          double e0_app = eigval0(ops, block);

          // Log.out("e0_mat: {}, e0_app: {}", e0_mat, e0_app);
          REQUIRE(std::abs(e0_mat - e0_app) < 1e-7);

          // Compute eigenvalues with real arithmitic
          if (block.irrep().isreal() && isreal(ops)) {
            auto H_real = matrix(ops, block, block);
            REQUIRE(arma::norm(H_real - H_real.t()) < 1e-12);

            arma::vec evals_mat_real;
            arma::eig_sym(evals_mat_real, H_real);
            REQUIRE(close(evals_mat_real, evals_mat));

            double e0_mat_real = evals_mat_real(0);
            double e0_app_real = eigval0(ops, block);
            REQUIRE(std::abs(e0_mat_real - e0_app_real) < 1e-7);
          }
        }
      }
    }
  }
}

void test_tj_symmetric_apply_chains(int64_t n_sites) {
  using namespace xdiag::testcases::tj;
  using namespace xdiag::testcases::electron;

  Log("tj_symmetric_apply: tJ chain, symmetric apply test, n_sites: {}",
      n_sites);
  auto ops = tJchain(n_sites, 1.0, 5.0);
  auto [group, irreps, multiplicities] =
      get_cyclic_group_irreps_mult(n_sites);
  (void)multiplicities;
  test_apply_tj_symmetric(ops, group, irreps);
}

TEST_CASE("tj_symmetric_apply", "[tj]") {
  using namespace xdiag::testcases::tj;
  using namespace xdiag::testcases::electron;

  // Test linear chains
  for (int64_t n_sites = 2; n_sites < 8; ++n_sites) {
    test_tj_symmetric_apply_chains(n_sites);
  }
  {
    // test a 3x3 triangular lattice
    Log("tj_symmetric_apply: tJ 3x3 triangular, symmetric apply test");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.toml";

    auto fl = FileToml(lfile);
    auto ops = fl["Interactions"].as<OpSum>();
    ops["T"] = 1.0;
    ops["U"] = 5.0;
    auto group = fl["Symmetries"].as<PermutationGroup>();

    std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
        {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
        {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
        {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
        {"Y.C1.A", 6}};

    std::vector<Representation> irreps;
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(fl[name].as<Representation>());
      (void)mult;
    }
    test_apply_tj_symmetric(ops, group, irreps);
  }

  {
    // test a 3x3 triangular lattice with complex flux
    Log("tj_symmetric_apply: tJ 3x3 triangular staggered, symmetric apply "
        "test, complex");
    std::string lfile = XDIAG_DIRECTORY
        "/misc/data/triangular.9.tup.phi.tdn.nphi.sublattices.tsl.toml";

    auto fl = FileToml(lfile);
    auto ops = fl["Interactions"].as<OpSum>();
    std::vector<double> etas{0.0, 0.1, 0.2, 0.3};
    auto group = fl["Symmetries"].as<PermutationGroup>();

    std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
        {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
        {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
        {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
        {"Y.C1.A", 6}};

    std::vector<Representation> irreps;
    std::vector<int64_t> multiplicities;
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(fl[name].as<Representation>());
      multiplicities.push_back(mult);
    }

    for (auto eta : etas) {
      ops["TPHI"] = complex(cos(eta * M_PI), sin(eta * M_PI));
      ops["JPHI"] = complex(cos(2 * eta * M_PI), sin(2 * eta * M_PI));

      test_apply_tj_symmetric(ops, group, irreps);
    }
  }
}