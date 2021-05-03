#include "../catch.hpp"

#include <iostream>

#include "testcases_electron.h"
#include <hydra/all.h>

using namespace hydra;


template <class bit_t>
void test_symmetric_spectra(BondList bondlist, Couplings couplings,
                            SpaceGroup<bit_t> space_group,
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
            matrix_cplx(bondlist, couplings, electron_nosym, electron_nosym);
        REQUIRE(lila::close(H_nosym, lila::Herm(H_nosym)));
        auto eigs_nosym = lila::EigenvaluesSym(H_nosym);

        lila::Vector<double> eigs_sym;
        for (int k = 0; k < (int)irreps.size(); ++k) {
          auto irrep = irreps[k];
          int multiplicity = multiplicities[k];

          auto electron =
              ElectronSymmetric<bit_t>(n_sites, nup, ndn, space_group, irrep);
          // HydraLog.out(
          //     "nup: {}, ndn: {}, k: {}, mult: {}, dim_nosym: {}, dim_sym: {}",
          //     nup, ndn, k, multiplicity, electron_nosym.size(),
          //     electron.size());

          if (electron.size() > 0) {

            // Compute partial spectrum from symmetrized block
            auto H_sym = matrix_cplx(bondlist, couplings, electron, electron);
            REQUIRE(lila::close(H_sym, lila::Herm(H_sym)));
            auto eigs_sym_k = lila::EigenvaluesSym(H_sym);

            // append all the eigenvalues with multiplicity
            for (auto eig : eigs_sym_k)
              for (int i = 0; i < multiplicity; ++i)
                eigs_sym.push_back(eig);
          }
        }
        std::sort(eigs_sym.begin(), eigs_sym.end());

        // Check if all eigenvalues agree
        REQUIRE(lila::close(eigs_sym, eigs_nosym));
      }
    }
  }
}

template <class bit_t>
void test_hubbard_symmetric_spectrum_chains(int n_sites) {
  using namespace hydra::testcases::electron;
  HydraLog.out("Hubbard chain, symmetric spectra test, n_sites: {}", n_sites);
  auto [bondlist, couplings] = get_linear_chain(n_sites, 1.0, 5.0);
  auto [space_group, irreps, multiplicities] =
      get_cyclic_group_irreps_mult<bit_t>(n_sites);
  test_symmetric_spectra(bondlist, couplings, space_group, irreps,
                         multiplicities);
}

TEST_CASE("electron_symmetric_matrix", "[models]") {
  using namespace hydra::testcases::electron;
  
  // Check matrices agains Weisse & Fehske
  int n_sites = 4;
  int nup = 3;
  int ndn = 2;
  double t = 1.0;
  double U = 5.0;
  auto [bondlist, couplings] = get_linear_chain(n_sites, t, U);
  auto [space_group, irreps, multiplicities] =
      get_cyclic_group_irreps_mult<uint16>(n_sites);
  for (int k = 0; k < (int)irreps.size(); ++k) {
    auto irrep = irreps[k];
    auto electron =
        ElectronSymmetric<uint16>(n_sites, nup, ndn, space_group, irrep);
    auto H_sym = matrix_cplx(bondlist, couplings, electron, electron);
    complex U2 = 2 * U;
    complex UU = U;
    complex tp = t;
    complex tm = -t;
    complex it = complex(0, t);
    if (k == 0) {
      lila::Matrix<complex> H_correct = {
          {U2, tm, tm, tp, tp, 0.}, {tm, U2, tm, tm, 0., tp},
          {tm, tm, U2, 0., tm, tm}, {tp, tm, 0., UU, tm, tp},
          {tp, 0., tm, tm, UU, tm}, {0., tp, tm, tp, tm, UU}};
      REQUIRE(lila::close(H_correct, H_sym));
    }
    if (k == 1) {
      lila::Matrix<complex> H_correct = {
          {U2, tm, it, it, tp, 0.},       {tm, U2, tm, tm, 2. * it, tp},
          {-it, tm, U2, 0., tm, it},      {-it, tm, 0., UU, tm, it},
          {tp, -2. * it, tm, tm, UU, tm}, {0., tp, -it, -it, tm, UU}};
      REQUIRE(lila::close(H_correct, H_sym));
    }
    if (k == 2) {
      lila::Matrix<complex> H_correct = {
          {U2, tm, tp, tm, tp, 0.}, {tm, U2, tm, tm, 0., tp},
          {tp, tm, U2, 0., tm, tp}, {tm, tm, 0., UU, tm, tm},
          {tp, 0., tm, tm, UU, tm}, {0., tp, tp, tm, tm, UU}};
      REQUIRE(lila::close(H_correct, H_sym));
    }
    if (k == 3) {
      lila::Matrix<complex> H_correct = {
          {U2, tm, -it, -it, tp, 0.},    {tm, U2, tm, tm, -2. * it, tp},
          {it, tm, U2, 0., tm, -it},     {it, tm, 0., UU, tm, -it},
          {tp, 2. * it, tm, tm, UU, tm}, {0., tp, it, it, tm, UU}};
      REQUIRE(lila::close(H_correct, H_sym));
    }
  }

  // Test linear chains
  for (int n_sites = 2; n_sites < 7; ++n_sites) {
    test_hubbard_symmetric_spectrum_chains<hydra::uint16>(n_sites);
    test_hubbard_symmetric_spectrum_chains<hydra::uint32>(n_sites);
    test_hubbard_symmetric_spectrum_chains<hydra::uint64>(n_sites);
  }

  // test a 3x3 triangular lattice
  HydraLog.out("Hubbard 3x3 triangular, symmetric spectra test");
  using bit_t = uint16;
  std::string lfile = "data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.lat";

  bondlist = read_bondlist(lfile);
  couplings.clear();
  couplings["T"] = 1.0;
  couplings["U"] = 5.0;
  auto permutations = read_permutations(lfile);
  space_group = SpaceGroup<bit_t>(permutations);

  std::vector<std::pair<std::string, int>> rep_name_mult = {
      {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
      {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
      {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
      {"Y.C1.A", 6}};
  irreps.clear();
  multiplicities.clear();
  for (auto [name, mult] : rep_name_mult) {
    irreps.push_back(read_represenation(lfile, name));
    multiplicities.push_back(mult);
  }
  test_symmetric_spectra(bondlist, couplings, space_group, irreps,
                         multiplicities);
}
