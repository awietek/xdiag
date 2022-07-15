#include "../../catch.hpp"

#include "../electron/testcases_electron.h"
#include "../spinhalf/testcases_spinhalf.h"

#include <iostream>

#include <hydra/all.h>
using namespace hydra;

template <class bit_t>
void test_spinhalf_symmetric_spectra(BondList bondlist, Couplings couplings,
                                     PermutationGroup space_group,
                                     std::vector<Representation> irreps,
                                     std::vector<int> multiplicities) {
  int n_sites = space_group.n_sites();
  assert(irreps.size() == multiplicities.size());

  for (int nup = 0; nup <= n_sites; ++nup) {
    // Log("Spinhalf Symmetric N: {}, nup: {}", n_sites, nup);
    // Compute the full spectrum from non-symmetrized block
    auto spinhalf_nosym = Spinhalf<bit_t>(n_sites, nup);

    if (spinhalf_nosym.size() < 1000) {
      lila::Vector<double> eigs_sym;

      auto H_nosym =
          MatrixCplx(bondlist, couplings, spinhalf_nosym, spinhalf_nosym);
      REQUIRE(lila::close(H_nosym, lila::Herm(H_nosym)));
      auto eigs_nosym = lila::EigenvaluesSym(H_nosym);
      for (int k = 0; k < (int)irreps.size(); ++k) {
        auto irrep = irreps[k];
        int multiplicity = multiplicities[k];
        auto spinhalf = Spinhalf<bit_t>(n_sites, nup, space_group, irrep);
        // Log.out(
        //     "nup: {}, k: {}, mult: {}, dim_nosym: {}, dim_sym: "
        //     "{} ", nup, k, multiplicity, spinhalf_nosym.size(),
        //     spinhalf.size());
        if (spinhalf.size() > 0) {

          // Compute partial spectrum from symmetrized block
          auto H_sym = MatrixCplx(bondlist, couplings, spinhalf, spinhalf);
          REQUIRE(lila::close(H_sym, lila::Herm(H_sym)));
          auto eigs_sym_k = lila::EigenvaluesSym(H_sym);
          // Check whether results are the same for real blocks
          if (!is_complex(spinhalf.irrep()) && !(is_complex(couplings))) {
            auto H_sym_real =
                MatrixReal(bondlist, couplings, spinhalf, spinhalf);
            auto eigs_sym_k_real = lila::EigenvaluesSym(H_sym);
            REQUIRE(lila::close(eigs_sym_k, eigs_sym_k_real));
          }
          // append all the eigenvalues with multiplicity
          for (auto eig : eigs_sym_k)
            for (int i = 0; i < multiplicity; ++i)
              eigs_sym.push_back(eig);
        }
      }
      std::sort(eigs_sym.begin(), eigs_sym.end());

      // if (!lila::close(eigs_sym, eigs_nosym)) {
      // LilaPrint(eigs_sym);
      // LilaPrint(eigs_nosym);
      // }
      REQUIRE(lila::close(eigs_sym, eigs_nosym));
    }
  }
}

template <class bit_t>
void test_spinhalf_symmetric_spectra_no_sz(BondList bondlist,
                                           Couplings couplings,
                                           PermutationGroup space_group,
                                           std::vector<Representation> irreps,
                                           std::vector<int> multiplicities) {
  int n_sites = space_group.n_sites();
  assert(irreps.size() == multiplicities.size());

  // Log("Spinhalf Symmetric N: {}, nup: {}", n_sites, nup);
  // Compute the full spectrum from non-symmetrized block
  auto spinhalf_nosym = Spinhalf<bit_t>(n_sites);

  if (spinhalf_nosym.size() < 1000) {
    lila::Vector<double> eigs_sym;

    auto H_nosym =
        MatrixCplx(bondlist, couplings, spinhalf_nosym, spinhalf_nosym);
    REQUIRE(lila::close(H_nosym, lila::Herm(H_nosym)));
    auto eigs_nosym = lila::EigenvaluesSym(H_nosym);
    for (int k = 0; k < (int)irreps.size(); ++k) {
      auto irrep = irreps[k];
      int multiplicity = multiplicities[k];
      auto spinhalf = Spinhalf<bit_t>(n_sites, space_group, irrep);
      if (spinhalf.size() > 0) {

        // Compute partial spectrum from symmetrized block
        auto H_sym = MatrixCplx(bondlist, couplings, spinhalf, spinhalf);
        REQUIRE(lila::close(H_sym, lila::Herm(H_sym)));
        auto eigs_sym_k = lila::EigenvaluesSym(H_sym);
        // LilaPrint(eigs_sym_k);

        auto eigs_sym_k_sz = lila::Vector<double>();
        for (int nup = 0; nup <= n_sites; ++nup) {
          auto spinhalf_sz = Spinhalf<bit_t>(n_sites, nup, space_group, irrep);
          auto H_sym_sz =
              MatrixCplx(bondlist, couplings, spinhalf_sz, spinhalf_sz);
          auto es = lila::EigenvaluesSym(H_sym_sz);
          // LilaPrint(es);
          for (auto e : es)
            eigs_sym_k_sz.push_back(e);
        }
        std::sort(eigs_sym_k_sz.begin(), eigs_sym_k_sz.end());
        // LilaPrint(eigs_sym_k);
        // LilaPrint(eigs_sym_k_sz);

        REQUIRE(lila::close(eigs_sym_k, eigs_sym_k_sz));

        // Check whether results are the same for real blocks
        if (!is_complex(spinhalf.irrep()) && !(is_complex(couplings))) {
          auto H_sym_real = MatrixReal(bondlist, couplings, spinhalf, spinhalf);
          auto eigs_sym_k_real = lila::EigenvaluesSym(H_sym);
          REQUIRE(lila::close(eigs_sym_k, eigs_sym_k_real));
        }
        // append all the eigenvalues with multiplicity
        for (auto eig : eigs_sym_k)
          for (int i = 0; i < multiplicity; ++i)
            eigs_sym.push_back(eig);
      }
    }
    std::sort(eigs_sym.begin(), eigs_sym.end());
    // LilaPrint(eigs_sym);
    // LilaPrint(eigs_nosym);
    REQUIRE(lila::close(eigs_sym, eigs_nosym));
  }
}

template <class bit_t>
void test_spinhalf_symmetric_spectrum_chains(int n_sites) {
  using namespace hydra::testcases::spinhalf;
  using hydra::testcases::electron::get_cyclic_group_irreps_mult;

  // Without Heisenberg term
  Log.out("spinhalf_symmetric_matrix: HB chain, N: {}", n_sites);
  auto [space_group, irreps, multiplicities] =
      get_cyclic_group_irreps_mult(n_sites);
  auto [bondlist, couplings] = HBchain(n_sites, 1.0, 1.0);
  test_spinhalf_symmetric_spectra<bit_t>(bondlist, couplings, space_group,
                                         irreps, multiplicities);
  test_spinhalf_symmetric_spectra_no_sz<bit_t>(bondlist, couplings, space_group,
                                               irreps, multiplicities);
}

TEST_CASE("spinhalf_symmetric_matrix", "[blocks][spinhalf_symmetric]") {

  // Test linear Heisenberg chains
  for (int n_sites = 3; n_sites < 7; ++n_sites) {
    test_spinhalf_symmetric_spectrum_chains<uint16_t>(n_sites);
    test_spinhalf_symmetric_spectrum_chains<uint32_t>(n_sites);
    test_spinhalf_symmetric_spectrum_chains<uint64_t>(n_sites);
  }

  // test a 3x3 triangular lattice
  {
    Log("spinhalf_symmetric_matrix: Triangular 3x3");
    std::string lfile = "data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.lat";

    auto bondlist = read_bondlist(lfile);
    Couplings couplings;
    couplings["Jz1"] = 1.00;
    couplings["Jz2"] = 0.23;
    couplings["Jx1"] = 0.76;
    couplings["Jx2"] = 0.46;

    std::vector<std::pair<std::string, int>> rep_name_mult = {
        {"Gamma.D6.A1", 1}, {"Gamma.D6.A2", 1}, {"Gamma.D6.B1", 1},
        {"Gamma.D6.B2", 1}, {"Gamma.D6.E1", 2}, {"Gamma.D6.E2", 2},
        {"K.D3.A1", 2},     {"K.D3.A2", 2},     {"K.D3.E", 4},
        {"Y.D1.A", 6},      {"Y.D1.B", 6}};

    auto permutations = hydra::read_permutations(lfile);
    auto space_group = PermutationGroup(permutations);

    std::vector<Representation> irreps;
    std::vector<int> multiplicities;
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(read_represenation(lfile, name));
      multiplicities.push_back(mult);
    }
    test_spinhalf_symmetric_spectra<uint16_t>(bondlist, couplings, space_group,
                                              irreps, multiplicities);
    test_spinhalf_symmetric_spectra_no_sz<uint32_t>(
        bondlist, couplings, space_group, irreps, multiplicities);
  }

  // test J1-J2-Jchi triangular lattice
  {
    Log("spinhalf_symmetric_matrix: Triangular J1J2Jchi N=12");
    std::string lfile = "data/triangular.j1j2jch/"
                        "triangular.12.j1j2jch.sublattices.fsl.lat";

    auto bondlist = read_bondlist(lfile);
    Couplings couplings;
    couplings["J1"] = 1.00;
    couplings["J2"] = 0.15;
    couplings["Jchi"] = 0.09;

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

    auto permutations = hydra::read_permutations(lfile);
    auto space_group = PermutationGroup(permutations);

    int n_sites = 12;
    int n_up = 6;
    for (auto [name, energy] : rep_name_mult) {
      auto irrep = read_represenation(lfile, name);
      auto spinhalf = Spinhalf<uint16_t>(n_sites, n_up, space_group, irrep);
      auto H = MatrixCplx(bondlist, couplings, spinhalf, spinhalf);
      REQUIRE(lila::close(H, lila::Herm(H)));

      auto eigs = lila::EigenvaluesSym(H);
      // Log("{} {:.18f} {:.18f}", name, eigs(0), energy);

      REQUIRE(std::abs(eigs(0) - energy) < 1e-10);
    }
  }
}
