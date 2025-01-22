#include "../../catch.hpp"

#include <iostream>

#include "../electron/testcases_electron.hpp"
#include "../spinhalf/testcases_spinhalf.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/apply.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/utils/close.hpp>

using namespace xdiag;

void test_spinhalf_symmetric_spectra(OpSum ops, int64_t nsites,
                                     std::vector<Representation> irreps,
                                     std::vector<int64_t> multiplicities) {
  assert(irreps.size() == multiplicities.size());

  for (int64_t nup = 3; nup <= nsites; ++nup) {
    // Log("Spinhalf Symmetric N: {}, nup: {}", nsites, nup);
    // Compute the full spectrum from non-symmetrized block
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
        // Log.out("nup: {}, k: {}, mult: {}, dim_nosym: {}, dim_sym: "
        //         "{} ",
        //         nup, k, multiplicity, spinhalf_nosym.size(),
        //         spinhalf.size());
        if (spinhalf.size() > 0) {

          // Compute partial spectrum from symmetrized block
          auto H_sym = matrixC(ops, spinhalf, spinhalf);
          REQUIRE(arma::norm(H_sym - H_sym.t()) < 1e-12);
          arma::vec eigs_sym_k;
          arma::eig_sym(eigs_sym_k, H_sym);

          // Check whether results are the same for real blocks
          if (isreal(spinhalf) && isreal(ops)) {
            auto H_sym_real = matrix(ops, spinhalf, spinhalf);
            arma::vec eigs_sym_k_real;
            arma::eig_sym(eigs_sym_k_real, H_sym_real);

            REQUIRE(close(eigs_sym_k, eigs_sym_k_real));
          }
          // append all the eigenvalues with multiplicity
          for (auto eig : eigs_sym_k)
            for (int64_t i = 0; i < multiplicity; ++i)
              eigs_sym.push_back(eig);
        }
      }
      std::sort(eigs_sym.begin(), eigs_sym.end());

      REQUIRE(close(arma::vec(eigs_sym), eigs_nosym));
    }
  }
}

void test_spinhalf_symmetric_spectra_no_sz(
    OpSum ops, int64_t nsites, std::vector<Representation> irreps,
    std::vector<int64_t> multiplicities) {
  assert(irreps.size() == multiplicities.size());

  // Log("Spinhalf Symmetric N: {}, nup: {}", nsites, nup);
  // Compute the full spectrum from non-symmetrized block
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

        // Compute partial spectrum from symmetrized block
        auto H_sym = matrixC(ops, spinhalf, spinhalf);
        REQUIRE(arma::norm(H_sym - H_sym.t()) < 1e-12);

        arma::vec eigs_sym_k;
        arma::eig_sym(eigs_sym_k, H_sym);

        auto eigs_sym_k_sz = std::vector<double>();
        for (int64_t nup = 0; nup <= nsites; ++nup) {
          auto spinhalf_sz = Spinhalf(nsites, nup, irrep);
          auto H_sym_sz = matrixC(ops, spinhalf_sz, spinhalf_sz);
          arma::vec es;
          arma::eig_sym(es, H_sym_sz);
          for (auto e : es)
            eigs_sym_k_sz.push_back(e);
        }
        std::sort(eigs_sym_k_sz.begin(), eigs_sym_k_sz.end());

        REQUIRE(close(eigs_sym_k, arma::vec(eigs_sym_k_sz)));

        // Check whether results are the same for real blocks
        if (isreal(spinhalf) && isreal(ops)) {
          auto H_sym_real = matrix(ops, spinhalf, spinhalf);

          arma::vec eigs_sym_k_real;
          arma::eig_sym(eigs_sym_k_real, H_sym);
          REQUIRE(close(eigs_sym_k, eigs_sym_k_real));
        }
        // append all the eigenvalues with multiplicity
        for (auto eig : eigs_sym_k)
          for (int64_t i = 0; i < multiplicity; ++i)
            eigs_sym.push_back(eig);
      }
    }
    std::sort(eigs_sym.begin(), eigs_sym.end());
    REQUIRE(close(arma::vec(eigs_sym), eigs_nosym));
  }
}

void test_spinhalf_symmetric_spectrum_chains(int64_t nsites) {
  using namespace xdiag::testcases::spinhalf;
  using xdiag::testcases::electron::get_cyclic_group_irreps_mult;

  // Without Heisenberg term
  Log.out("spinhalf_symmetric_matrix: HB chain, N: {}", nsites);
  auto [irreps, multiplicities] = get_cyclic_group_irreps_mult(nsites);
  auto ops = HBchain(nsites, 1.0, 1.0);
  test_spinhalf_symmetric_spectra(ops, nsites, irreps, multiplicities);
  test_spinhalf_symmetric_spectra_no_sz(ops, nsites, irreps, multiplicities);
}

TEST_CASE("spinhalf_symmetric_matrix", "[spinhalf]") try {

  // Test linear Heisenberg chains
  for (int64_t nsites = 3; nsites < 7; ++nsites) {
    test_spinhalf_symmetric_spectrum_chains(nsites);
  }

  // test a 3x3 triangular lattice
  {
    Log("spinhalf_symmetric_matrix: Triangular 3x3");
    std::string lfile = XDIAG_DIRECTORY
        "/misc/data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.toml";
    auto fl = FileToml(lfile);
    auto ops = fl["Interactions"].as<OpSum>();
    ops["Jz1"] = 1.00;
    ops["Jz2"] = 0.23;
    ops["Jx1"] = 0.76;
    ops["Jx2"] = 0.46;
    std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
        {"Gamma.D6.A1", 1}, {"Gamma.D6.A2", 1}, {"Gamma.D6.B1", 1},
        {"Gamma.D6.B2", 1}, {"Gamma.D6.E1", 2}, {"Gamma.D6.E2", 2},
        {"K.D3.A1", 2},     {"K.D3.A2", 2},     {"K.D3.E", 4},
        {"Y.D1.A", 6},      {"Y.D1.B", 6}};
    std::vector<Representation> irreps;
    std::vector<int64_t> multiplicities;
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(read_representation(fl, name));
      multiplicities.push_back(mult);
    }
    test_spinhalf_symmetric_spectra(ops, 9, irreps, multiplicities);
    test_spinhalf_symmetric_spectra_no_sz(ops, 9, irreps, multiplicities);
  }

  // test J1-J2-Jchi triangular lattice
  {
    Log("spinhalf_symmetric_matrix: Triangular J1J2Jchi N=12");
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/triangular.j1j2jch/"
                        "triangular.12.j1j2jch.sublattices.fsl.toml";

    auto fl = FileToml(lfile);
    auto ops = fl["Interactions"].as<OpSum>();
    ops["J1"] = 1.00;
    ops["J2"] = 0.15;
    ops["Jchi"] = -0.09;

    std::vector<std::pair<std::string, double>> rep_name_mult = {
        {"Gamma.C6.A", -6.9456000700824329641},
        {"Gamma.C6.B", -5.8410912437873072633},
        {"Gamma.C6.E1a", -3.8556417248355927541},
        {"Gamma.C6.E2a", -6.4157243358059030669},
        {"K.C3.A", -5.9197511811622431921},
        {"K.C3.Ea", -5.0281703836000861685},
        {"K.C3.Eb", -5.2045133473640809996},
        {"M.C2.A", -5.756684675081964464},
        {"M.C2.B", -5.7723510325561688816},
        {"X.C1.A", -5.9030627660522529965}};

    auto space_group = fl["Symmetries"].as<PermutationGroup>();
    int64_t nsites = 12;
    int64_t nup = 6;
    for (auto [name, energy] : rep_name_mult) {
      auto irrep = read_representation(fl, name);
      auto spinhalf = Spinhalf(nsites, nup, irrep);
      auto H = matrixC(ops, spinhalf, spinhalf);
      REQUIRE(arma::norm(H - H.t()) < 1e-12);

      arma::vec eigs;
      arma::eig_sym(eigs, H);

      // Log("{} {:.18f} {:.18f}", name, eigs(0), energy);

      REQUIRE(std::abs(eigs(0) - energy) < 1e-10);
    }
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
}
