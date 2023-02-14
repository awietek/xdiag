#include "../../catch.hpp"

#include <iostream>

#include "../electron/testcases_electron.h"
#include "../tj/testcases_tj.h"

#include <hydra/all.h>

using namespace hydra;

void test_spectra_tj_symmetric(BondList bondlist, PermutationGroup space_group,
                               std::vector<Representation> irreps,
                               std::vector<int> multiplicities) {
  int n_sites = space_group.n_sites();
  assert(irreps.size() == multiplicities.size());

  for (int nup = 0; nup <= n_sites; ++nup) {
    for (int ndn = 0; ndn <= n_sites; ++ndn) {

      if (nup + ndn > n_sites)
        continue;

      // Compute the full spectrum from non-symmetrized block
      auto tj_nosym = tJ(n_sites, nup, ndn);
      if (tj_nosym.size() < 1000) {

        auto H_nosym = matrix_cplx(bondlist, tj_nosym, tj_nosym);
        REQUIRE(arma::norm(H_nosym - H_nosym.t()) < 1e-12);
        arma::vec eigs_nosym;
        arma::eig_sym(eigs_nosym, H_nosym);

        std::vector<double> eigs_sym;
        for (int k = 0; k < (int)irreps.size(); ++k) {
          auto irrep = irreps[k];
          int multiplicity = multiplicities[k];

          auto tj = tJ(n_sites, nup, ndn, space_group, irrep);

          if (tj.size() > 0) {

            // Compute partial spectrum from symmetrized block
            auto H_sym = matrix_cplx(bondlist, tj, tj);
            REQUIRE(arma::norm(H_sym - H_sym.t()) < 1e-12);
            arma::vec eigs_sym_k;
            arma::eig_sym(eigs_sym_k, H_sym);

            // Check whether results are the same for real blocks
            if (is_real(tj.irrep()) && bondlist.is_real()) {
              auto H_sym_real = matrix_real(bondlist, tj, tj);
              arma::vec eigs_sym_k_real;
              arma::eig_sym(eigs_sym_k_real, H_sym_real);
              REQUIRE(close(eigs_sym_k, eigs_sym_k_real));
            }
            // append all the eigenvalues with multiplicity
            for (auto eig : eigs_sym_k)
              for (int i = 0; i < multiplicity; ++i)
                eigs_sym.push_back(eig);
          }
        }
        std::sort(eigs_sym.begin(), eigs_sym.end());

        // Check if all eigenvalues agree
        REQUIRE(close(arma::vec(eigs_sym), eigs_nosym));
      }
    }
  }
}

void test_tj_symmetric_spectrum_chains(int n_sites) {
  using namespace hydra::testcases::tj;
  using namespace hydra::testcases::electron;

  Log.out("tj_symmetric_matrix: tJ chain, symmetric spectra test, n_sites: {}",
          n_sites);
  auto bondlist = tJchain(n_sites, 1.0, 0.4);
  auto [space_group, irreps, multiplicities] =
      get_cyclic_group_irreps_mult(n_sites);
  test_spectra_tj_symmetric(bondlist, space_group, irreps, multiplicities);
}

TEST_CASE("tj_symmetric_matrix", "[blocks][tj]") {
  using namespace hydra::testcases::tj;
  using namespace hydra::testcases::electron;

  // Test linear chains
  for (int n_sites = 2; n_sites < 7; ++n_sites) {
    test_tj_symmetric_spectrum_chains(n_sites);
  }

  {
    // test a 3x3 triangular lattice
    Log("tj_symmetric_matrix: tJ 3x3 triangular, symmetric spectra test");
    std::string lfile =
        HYDRA_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.lat";

    auto bondlist = read_bondlist(lfile);
    bondlist["T"] = 1.0;
    bondlist["J"] = 0.4;
    auto permutations = hydra::read_permutations(lfile);
    auto space_group = PermutationGroup(permutations);

    std::vector<std::pair<std::string, int>> rep_name_mult = {
        {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
        {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
        {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
        {"Y.C1.A", 6}};

    std::vector<Representation> irreps;
    std::vector<int> multiplicities;
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(read_representation(lfile, name));
      multiplicities.push_back(mult);
    }
    test_spectra_tj_symmetric(bondlist, space_group, irreps, multiplicities);
  }

  {
    // test a 3x3 triangular lattice with complex flux
    Log("tj_symmetric_matrix: tJ 3x3 triangular staggered flux, "
        "symmetric spectra test, complex");
    std::string lfile = HYDRA_DIRECTORY
        "/misc/data/triangular.9.tup.phi.tdn.nphi.sublattices.tsl.lat";

    auto bondlist = read_bondlist(lfile);
    std::vector<double> etas{0.0, 0.1, 0.2, 0.3};
    auto permutations = hydra::read_permutations(lfile);
    auto space_group = PermutationGroup(permutations);

    std::vector<std::pair<std::string, int>> rep_name_mult = {
        {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
        {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
        {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
        {"Y.C1.A", 6}};

    std::vector<Representation> irreps;
    std::vector<int> multiplicities;
    for (auto [name, mult] : rep_name_mult) {
      irreps.push_back(read_representation(lfile, name));
      multiplicities.push_back(mult);
    }

    for (auto eta : etas) {
      // Log("eta: {:.2f}", eta);
      bondlist["TPHI"] = complex(cos(eta * M_PI), sin(eta * M_PI));
      bondlist["JPHI"] = complex(cos(2 * eta * M_PI), sin(2 * eta * M_PI));
      test_spectra_tj_symmetric(bondlist, space_group, irreps, multiplicities);
    }
  }
}
