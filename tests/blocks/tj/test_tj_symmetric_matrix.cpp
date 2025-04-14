// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <iostream>

#include "../electron/testcases_electron.hpp"
#include "../tj/testcases_tj.hpp"

#include <xdiag/algebra/matrix.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/algebra/isapprox.hpp>

using namespace xdiag;

void test_spectra_tj_symmetric(OpSum ops, int64_t nsites,
                               std::vector<Representation> irreps,
                               std::vector<int64_t> multiplicities) {
  // XDIAG_SHOW(ops);
  assert(irreps.size() == multiplicities.size());

  for (int64_t nup = 1; nup <= nsites; ++nup) {
    for (int64_t ndn = 1; ndn <= nsites; ++ndn) {

      if (nup + ndn > nsites)
        continue;
      // Compute the full spectrum from non-symmetrized block

      auto tj_nosym = tJ(nsites, nup, ndn);
      if (tj_nosym.size() < 1000) {
        auto H_nosym = matrixC(ops, tj_nosym, tj_nosym);

        REQUIRE(arma::norm(H_nosym - H_nosym.t()) < 1e-12);
        arma::vec eigs_nosym;
        arma::eig_sym(eigs_nosym, H_nosym);

        std::vector<double> eigs_sym;
        for (int64_t k = 0; k < (int64_t)irreps.size(); ++k) {
          auto irrep = irreps[k];
          int64_t multiplicity = multiplicities[k];
          auto tj = tJ(nsites, nup, ndn, irrep);
          if (tj.size() > 0) {

            // Compute partial spectrum from symmetrized block
            auto H_sym = matrixC(ops, tj, tj);

            // Log("nsites: {}, nup: {}, ndn: {}, k: {}", nsites, nup, ndn,
            // k); XDIAG_SHOW(irrep); XDIAG_SHOW(H_sym);

            REQUIRE(arma::norm(H_sym - H_sym.t()) < 1e-12);

            arma::vec eigs_sym_k;
            arma::eig_sym(eigs_sym_k, H_sym);

            // XDIAG_SHOW(eigs_sym_k);

            // Check whether results are the same for real blocks
            if (isreal(tj) && isreal(ops)) {
              auto H_sym_real = matrix(ops, tj, tj);
              arma::vec eigs_sym_k_real;
              arma::eig_sym(eigs_sym_k_real, H_sym_real);
              REQUIRE(isapprox(eigs_sym_k, eigs_sym_k_real));
            }
            // append all the eigenvalues with multiplicity
            for (auto eig : eigs_sym_k)
              for (int64_t i = 0; i < multiplicity; ++i)
                eigs_sym.push_back(eig);
          }
        }
        std::sort(eigs_sym.begin(), eigs_sym.end());

        // Check if all eigenvalues agree
        // XDIAG_SHOW(eigs_sym);
        // XDIAG_SHOW(eigs_nosym);
        REQUIRE(isapprox(arma::vec(eigs_sym), eigs_nosym));
      }
    }
  }
}

void test_tj_symmetric_spectrum_chains(int64_t nsites) {
  using namespace xdiag::testcases::tj;
  using namespace xdiag::testcases::electron;

  Log.out("tj_symmetric_matrix: tJ chain, symmetric spectra test, nsites: {}",
          nsites);
  auto ops = tJchain(nsites, 0.0, 1.0);
  auto [irreps, multiplicities] = get_cyclic_group_irreps_mult(nsites);
  test_spectra_tj_symmetric(ops, nsites, irreps, multiplicities);
}

TEST_CASE("tj_symmetric_matrix", "[tj]") try {
  using namespace xdiag::testcases::tj;
  using namespace xdiag::testcases::electron;

  for (int64_t nsites = 2; nsites < 7; ++nsites) {
    Log.out(
        "tj_symmetric_matrix: HB chain, symmetric spectra test, nsites: {}",
        nsites);
    OpSum ops;
    for (int64_t s = 0; s < nsites; ++s) {
      ops += Op("tJSdotS", {s, (s + 1) % nsites});
    }
    auto [irreps, multiplicities] = get_cyclic_group_irreps_mult(nsites);
    test_spectra_tj_symmetric(ops, nsites, irreps, multiplicities);
  }

  // Test linear chains
  for (int64_t nsites = 2; nsites < 7; ++nsites) {
    test_tj_symmetric_spectrum_chains(nsites);
  }

  {
    // test a 8 site square lattice Heisenberg model
    Log("tj_symmetric_matrix: 8 site square lattice HB model");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/square.8.heisenberg.2sl.toml";

    auto fl = FileToml(lfile);
    auto ops = fl["Interactions"].as<OpSum>();
    std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
        {"Gamma.D4.A1", 1}, {"Gamma.D4.A2", 1}, {"Gamma.D4.B1", 1},
        {"Gamma.D4.B2", 1}, {"Gamma.D4.E", 2},  {"M.D4.A1", 1},
        {"M.D4.A2", 1},     {"M.D4.B1", 1},     {"M.D4.B2", 1},
        {"M.D4.E", 2},      {"Sigma.D1.A", 4},  {"Sigma.D1.B", 4},
        {"X.D2.A1", 2},     {"X.D2.A2", 2},     {"X.D2.B1", 2},
        {"X.D2.B2", 2}};

    std::vector<Representation> irreps;
    std::vector<int64_t> multiplicities;
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(read_representation(fl, name));
      multiplicities.push_back(mult);
    }

    ops["J"] = 1.0;
    test_spectra_tj_symmetric(ops, 8, irreps, multiplicities);
  }

  {
    // test a 3x3 triangular lattice
    Log("tj_symmetric_matrix: tJ 3x3 triangular, symmetric spectra test");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.toml";

    auto fl = FileToml(lfile);
    auto ops = fl["Interactions"].as<OpSum>();
    ops["T"] = 1.0;
    ops["J"] = 0.4;
    auto group = fl["Symmetries"].as<PermutationGroup>();

    std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
        {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
        {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
        {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
        {"Y.C1.A", 6}};

    std::vector<Representation> irreps;
    std::vector<int64_t> multiplicities;
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(read_representation(fl, name));
      multiplicities.push_back(mult);
    }
    test_spectra_tj_symmetric(ops, 9, irreps, multiplicities);
  }

  {
    // test a 3x3 triangular lattice with complex flux
    Log("tj_symmetric_matrix: tJ 3x3 triangular staggered flux, "
        "symmetric spectra test, complex");
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
      irreps.push_back(read_representation(fl, name));
      multiplicities.push_back(mult);
    }

    for (auto eta : etas) {
      Log("eta: {:.2f}", eta);
      ops["TPHI"] = 1.0; // complex(cos(eta * M_PI), sin(eta * M_PI));
      ops["JPHI"] = 0.4; // complex(cos(2 * eta * M_PI), sin(2 * eta * M_PI));
      test_spectra_tj_symmetric(ops, 9, irreps, multiplicities);
    }
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
}
