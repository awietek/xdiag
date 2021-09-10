#include "../../catch.hpp"

#include <iostream>

#include "../electron/testcases_electron.h"
#include <hydra/all.h>

using namespace hydra;

template <class bit_t>
void test_symmetric_spectra(BondList bondlist, Couplings couplings,
                            PermutationGroup space_group,
                            std::vector<Representation> irreps,
                            std::vector<int> multiplicities) {
  int n_sites = space_group.n_sites();
  assert(irreps.size() == multiplicities.size());

  for (int nup = 0; nup <= n_sites; ++nup) {
    for (int ndn = 0; ndn <= n_sites; ++ndn) {

      // Compute the full spectrum from non-symmetrized block
      auto electron_nosym = Electron<bit_t>(n_sites, nup, ndn);
      if (electron_nosym.size() < 1000) {

        auto H_nosym =
            MatrixCplx(bondlist, couplings, electron_nosym, electron_nosym);
        REQUIRE(lila::close(H_nosym, lila::Herm(H_nosym)));
        auto eigs_nosym = lila::EigenvaluesSym(H_nosym);

        lila::Vector<double> eigs_sym;
        for (int k = 0; k < (int)irreps.size(); ++k) {
          auto irrep = irreps[k];
          int multiplicity = multiplicities[k];

          auto electron = ElectronSymmetricSimple<bit_t>(n_sites, nup, ndn,
                                                         space_group, irrep);
          // lila::Log.out(
          //     "nup: {}, ndn: {}, k: {}, mult: {}, dim_nosym: {}, dim_sym:
          //     {}", nup, ndn, k, multiplicity, electron_nosym.size(),
          //     electron.size());

          if (electron.size() > 0) {

            // Compute partial spectrum from symmetrized block
            auto H_sym = MatrixCplx(bondlist, couplings, electron, electron);
            REQUIRE(lila::close(H_sym, lila::Herm(H_sym)));
            auto eigs_sym_k = lila::EigenvaluesSym(H_sym);

            // Check whether results are the same for real blocks
            if (!is_complex(electron.irrep()) && !(is_complex(couplings))) {
              auto H_sym_real =
                  MatrixReal(bondlist, couplings, electron, electron);
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

        // Check if all eigenvalues agree
        // lila::Log.out("{} {} {}", nup, ndn, eigs_sym(0));
        REQUIRE(lila::close(eigs_sym, eigs_nosym));
      }
    }
  }
}

template <class bit_t>
void test_hubbard_symmetric_spectrum_chains(int n_sites) {
  using namespace hydra::testcases::electron;

  // Without Heisenberg term
  lila::Log.out("Hubbard chain, symmetric spectra test, n_sites: {}", n_sites);
  auto [space_group, irreps, multiplicities] =
      get_cyclic_group_irreps_mult(n_sites);
  auto [bondlist, couplings] = get_linear_chain(n_sites, 1.0, 5.0);
  test_symmetric_spectra<bit_t>(bondlist, couplings, space_group, irreps,
                                multiplicities);

  // With Heisenberg term
  lila::Log.out(
      "Hubbard chain, symmetric spectra test, n_sites: {} (+ Heisenberg terms)",
      n_sites);
  auto [bondlist_hb, couplings_hb] =
      get_linear_chain_hb(n_sites, 1.0, 5.0, 0.4);
  test_symmetric_spectra<bit_t>(bondlist_hb, couplings_hb, space_group, irreps,
                                multiplicities);
}

TEST_CASE("ElectronSymmetric_Matrix", "[models][ElectronSymmetric]") {
  using namespace hydra::testcases::electron;

  lila::Log.out("ElectronSymmetric_Matrix <-> ElectronSymmetricSimple_Matrix "
                "cross-check");
  for (int n_sites = 0; n_sites < 7; ++n_sites) {
    lila::Log.out("N: {}", n_sites);
    for (int nup = 0; nup < n_sites; ++nup) {
      for (int ndn = 0; ndn < n_sites; ++ndn) {
        auto [space_group, irreps, multiplicities] =
            get_cyclic_group_irreps_mult(n_sites);

        BondList bonds;
        for (int i = 0; i < n_sites; ++i) {
          bonds << Bond("HOP", "T", {i, (i + 1) % n_sites});
        }
        Couplings cpls;
        cpls["T"] = 1.0;
        cpls["U"] = 5.0;

        for (auto irrep : irreps) {

          auto block1 =
              ElectronSymmetric(n_sites, nup, ndn, space_group, irrep);
          auto block2 =
              ElectronSymmetricSimple(n_sites, nup, ndn, space_group, irrep);

          auto Hnew = MatrixCplx(bonds, cpls, block1, block1);
          auto Hsimple = MatrixCplx(bonds, cpls, block2, block2);

          // LilaPrint(Hnew);
          // LilaPrint(Hsimple);
          REQUIRE(lila::close(Hnew, Hsimple));
        }
      }
    }
  }
}
