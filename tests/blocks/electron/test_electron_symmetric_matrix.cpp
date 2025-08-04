// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../../catch.hpp"

#include <iostream>

#include "../electron/testcases_electron.hpp"
#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/algebra/isapprox.hpp>

using namespace xdiag;

void test_electron_symmetric_spectra_no_np(
    OpSum opsum, int64_t nsites, std::vector<Representation> irreps,
    std::vector<int64_t> multiplicities) try {
  (void)multiplicities;
  auto block_total = Electron(nsites);
  if (block_total.size() < 1000) {

    auto H_total = matrixC(opsum, block_total, block_total);
    REQUIRE(H_total.is_hermitian(1e-8));

    arma::vec eigs_total;
    arma::eig_sym(eigs_total, H_total);

    std::vector<double> eigs_all;
    for (auto irrep : irreps) {
      auto block_no_np = Electron(nsites, irrep);
      auto H_no_np = matrixC(opsum, block_no_np, block_no_np);

      REQUIRE(H_no_np.is_hermitian(1e-8));

      arma::vec eigs_no_np;
      arma::eig_sym(eigs_no_np, H_no_np);

      std::vector<double> eigs_np_all;
      for (int64_t nup = 0; nup <= nsites; ++nup) {
        for (int64_t ndn = 0; ndn <= nsites; ++ndn) {
          auto block_np = Electron(nsites, nup, ndn, irrep);
          auto H_np = matrixC(opsum, block_np, block_np);
          REQUIRE(H_np.is_hermitian(1e-8));

          arma::vec eigs_np;
          arma::eig_sym(eigs_np, H_np);

          for (auto e : eigs_np) {
            eigs_np_all.push_back(e);
            eigs_all.push_back(e);
          }
        }
      }
      std::sort(eigs_np_all.begin(), eigs_np_all.end());
      // XDIAG_SHOW(eigs_no_np);
      // XDIAG_SHOW(arma::vec(eigs_np_all));

      REQUIRE(isapprox(eigs_no_np, arma::vec(eigs_np_all)));
    }
    std::sort(eigs_all.begin(), eigs_all.end());
    REQUIRE(isapprox(eigs_total, arma::vec(eigs_all)));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void test_electron_symmetric_spectra(OpSum opsum, int64_t nsites,
                                     std::vector<Representation> irreps,
                                     std::vector<int64_t> multiplicities) try {
  assert(irreps.size() == multiplicities.size());

  for (int64_t nup = 0; nup <= nsites; ++nup) {
    for (int64_t ndn = 0; ndn <= nsites; ++ndn) {

      // Compute the full spectrum from non-symmetrized block
      auto electron_nosym = Electron(nsites, nup, ndn);
      if (electron_nosym.size() < 1000) {

        auto H_nosym = matrixC(opsum, electron_nosym, electron_nosym);
        REQUIRE(H_nosym.is_hermitian(1e-8));
        arma::vec eigs_nosym;
        arma::eig_sym(eigs_nosym, H_nosym);

        std::vector<double> eigs_sym;
        for (int64_t k = 0; k < (int64_t)irreps.size(); ++k) {
          auto irrep = irreps[k];
          int64_t multiplicity = multiplicities[k];

          auto electron = Electron(nsites, nup, ndn, irrep);
          // Log.out(
          //     "nup: {}, ndn: {}, k: {}, mult: {}, dim_nosym: {}, dim_sym:"
          //     "{} ",
          //     nup, ndn, k, multiplicity, electron_nosym.size(),
          //     electron.size());

          if (electron.size() > 0) {

            // Compute partial spectrum from symmetrized block
            auto H_sym = matrixC(opsum, electron, electron);
            REQUIRE(arma::norm(H_sym - H_sym.t()) < 1e-12);

            // REQUIRE(H_sym.is_hermitian(1e-7));
            arma::vec eigs_sym_k;
            arma::eig_sym(eigs_sym_k, H_sym);

            // Check whether results are the same for real blocks
            if (isreal(electron) && isreal(opsum)) {
              auto H_sym_real = matrix(opsum, electron, electron);
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
        // Log.out("{} {} {} {}", nup, ndn, eigs_sym(0), eigs_nosym(0));

        // if (!isapprox(arma::vec(eigs_sym), eigs_nosym)){
        //   XDIAG_SHOW(arma::norm(arma::vec(eigs_sym) - eigs_nosym));
        //   XDIAG_SHOW(eigs_sym);
        //   XDIAG_SHOW(eigs_nosym);
        // }
        REQUIRE(isapprox(arma::vec(eigs_sym), eigs_nosym));
      }
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void test_hubbard_symmetric_spectrum_chains(int64_t nsites) try {
  using namespace xdiag::testcases::electron;

  auto [irreps, multiplicities] = get_cyclic_group_irreps_mult(nsites);

  // Without Heisenberg term
  Log.out("electron_symmetric_matrix: Hubbard chain, nsites: {}", nsites);
  auto opsum = get_linear_chain(nsites, 1.0, 5.0);
  test_electron_symmetric_spectra(opsum, nsites, irreps, multiplicities);
  test_electron_symmetric_spectra_no_np(opsum, nsites, irreps, multiplicities);

  // With Heisenberg term
  Log("electron_symmetric_matrix: Hubbard chain, nsites: {} (+ "
      "Heisenberg terms)",
      nsites);
  auto opsum_hb = get_linear_chain_hb(nsites, 0.4);
  test_electron_symmetric_spectra(opsum_hb, nsites, irreps, multiplicities);
  test_electron_symmetric_spectra_no_np(opsum_hb, nsites, irreps,
                                        multiplicities);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

TEST_CASE("electron_symmetric_matrix", "[electron]") try {
  using namespace xdiag::testcases::electron;
  //////////////////////////////////////////////////////////////////////////////////////

  // Check matrices agains Weisse & Fehske
  Log("electron_symmetric_matrix: Weisse & Fehske matrix");
  int64_t nsites = 4;
  int64_t nup = 3;
  int64_t ndn = 2;
  double t = 1.0;
  double U = 5.0;
  auto opsum = get_linear_chain(nsites, t, U);
  auto [irreps, multiplicities] = get_cyclic_group_irreps_mult(nsites);

  for (int64_t k = 0; k < (int64_t)irreps.size(); ++k) {
    auto irrep = irreps[k];
    auto electron = Electron(nsites, nup, ndn, irrep);
    auto H_sym = matrixC(opsum, electron, electron);
    complex U2 = 2 * U;
    complex UU = U;
    complex tp = t;
    complex tm = -t;
    complex it = complex(0, t);
    if (k == 0) {
      arma::Mat<complex> H_correct = {
          {U2, tm, tm, tp, tp, 0.}, {tm, U2, tm, tm, 0., tp},
          {tm, tm, U2, 0., tm, tm}, {tp, tm, 0., UU, tm, tp},
          {tp, 0., tm, tm, UU, tm}, {0., tp, tm, tp, tm, UU}};
      REQUIRE(isapprox(H_correct, H_sym));
    }
    if (k == 1) {
      arma::Mat<complex> H_correct = {
          {U2, tm, it, it, tp, 0.},       {tm, U2, tm, tm, 2. * it, tp},
          {-it, tm, U2, 0., tm, it},      {-it, tm, 0., UU, tm, it},
          {tp, -2. * it, tm, tm, UU, tm}, {0., tp, -it, -it, tm, UU}};
      REQUIRE(isapprox(H_correct, H_sym));
    }
    if (k == 2) {
      arma::Mat<complex> H_correct = {
          {U2, tm, tp, tm, tp, 0.}, {tm, U2, tm, tm, 0., tp},
          {tp, tm, U2, 0., tm, tp}, {tm, tm, 0., UU, tm, tm},
          {tp, 0., tm, tm, UU, tm}, {0., tp, tp, tm, tm, UU}};
      REQUIRE(isapprox(H_correct, H_sym));
    }
    if (k == 3) {
      arma::Mat<complex> H_correct = {
          {U2, tm, -it, -it, tp, 0.},    {tm, U2, tm, tm, -2. * it, tp},
          {it, tm, U2, 0., tm, -it},     {it, tm, 0., UU, tm, -it},
          {tp, 2. * it, tm, tm, UU, tm}, {0., tp, it, it, tm, UU}};
      REQUIRE(isapprox(H_correct, H_sym));
    }
  }

  // Test linear chains
  for (int64_t nsites = 2; nsites < 7; ++nsites) {
    test_hubbard_symmetric_spectrum_chains(nsites);
    test_hubbard_symmetric_spectrum_chains(nsites);
    test_hubbard_symmetric_spectrum_chains(nsites);
  }

  // // test a 3x3 triangular lattice
  // Log("electron_symmetric_matrix: Hubbard 3x3 triangular");
  // std::string lfile =
  //     XDIAG_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.toml";

  // auto fl = FileToml(lfile);
  // opsum = fl["Interactions"].as<OpSum>();
  // opsum["T"] = 1.0;
  // opsum += "U" * Op("HubbardU");
  // opsum["U"] = 5.0;

  // std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
  //     {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
  //     {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
  //     {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
  //     {"Y.C1.A", 6}};
  // irreps.clear();
  // multiplicities.clear();
  // for (auto [name, mult] : rep_name_mult) {
  //   irreps.push_back(read_representation(fl, name));
  //   multiplicities.push_back(mult);
  // }
  // test_electron_symmetric_spectra(opsum, 9, irreps, multiplicities);

  // // test a 3x3 triangular lattice with Heisenberg terms
  // Log("electron_symmetric_matrix: Hubbard 3x3 triangular(+ Heisenberg terms)");
  // auto opsum_hb = opsum;
  // for (auto [cpl, op] : opsum) {
  //   if (op.type() == "Hop") {
  //     opsum_hb += "J" * Op("SdotS", {op[0], op[1]});
  //   }
  // }
  // opsum_hb["J"] = 0.4;
  // test_electron_symmetric_spectra(opsum_hb, 9, irreps, multiplicities);

  // // test a 3x3 triangular lattice with complex hoppings
  // {
  //   Log.out("electron_symmetric_matrix: Hubbard 3x3 triangular (complex)");
  //   std::string lfile = XDIAG_DIRECTORY
  //       "/misc/data/triangular.9.tup.phi.tdn.nphi.sublattices.tsl.toml";

  //   auto fl = FileToml(lfile);
  //   auto opsum = fl["Interactions"].as<OpSum>();
  //   opsum += "U" * Op("HubbardU");
  //   opsum["TPHI"] = complex(0.5, 0.5);
  //   opsum["JPHI"] = 0.;
  //   opsum["U"] = 5.0;

  //   std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
  //       {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
  //       {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
  //       {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
  //       {"Y.C1.A", 6}};
  //   irreps.clear();
  //   multiplicities.clear();
  //   for (auto [name, mult] : rep_name_mult) {
  //     irreps.push_back(read_representation(fl, name));
  //     multiplicities.push_back(mult);
  //   }
  //   test_electron_symmetric_spectra(opsum, 9, irreps, multiplicities);
  // }
} catch (Error const &e) {
  error_trace(e);
}
