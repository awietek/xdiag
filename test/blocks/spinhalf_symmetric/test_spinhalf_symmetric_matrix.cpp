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
    // lila::Log("Spinhalf Symmetric N: {}, nup: {}", n_sites, nup);
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
        auto spinhalf =
            SpinhalfSymmetric<bit_t>(n_sites, nup, space_group, irrep);
        // lila::Log.out(
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
      //   LilaPrint(eigs_sym);
      //   LilaPrint(eigs_nosym);
      // }
      REQUIRE(lila::close(eigs_sym, eigs_nosym));

    }
  }
}

template <class bit_t>
void test_spinhalf_symmetric_spectrum_chains(int n_sites) {
  using namespace hydra::testcases::spinhalf;
  using hydra::testcases::electron::get_cyclic_group_irreps_mult;

  // Without Heisenberg term
  lila::Log.out("HB chain, symmetric spectra test, n_sites: {}", n_sites);
  auto [space_group, irreps, multiplicities] =
      get_cyclic_group_irreps_mult(n_sites);
  auto [bondlist, couplings] = HBchain(n_sites, 1.0, 1.0);
  test_spinhalf_symmetric_spectra<bit_t>(bondlist, couplings, space_group,
                                         irreps, multiplicities);
}

TEST_CASE("spinhalf_symmetric_matrix", "[blocks][spinhalf]") {

  // Test linear Heisenberg chains
  for (int n_sites = 3; n_sites < 7; ++n_sites) {
    test_spinhalf_symmetric_spectrum_chains<uint16_t>(n_sites);
    test_spinhalf_symmetric_spectrum_chains<uint32_t>(n_sites);
    test_spinhalf_symmetric_spectrum_chains<uint64_t>(n_sites);
  }

  // test a 3x3 triangular lattice
  lila::Log("SpinhalfSymmetric spectra test: Triangular 3x3");
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
}
