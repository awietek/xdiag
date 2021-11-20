#include "../../catch.hpp"

#include <iostream>

#include "../electron/testcases_electron.h"
#include "../tj/testcases_tj.h"

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

      if (nup + ndn > n_sites)
        continue;

      // Compute the full spectrum from non-symmetrized block
      auto tj_nosym = tJ<bit_t>(n_sites, nup, ndn);
      if (tj_nosym.size() < 1000) {

        auto H_nosym = MatrixCplx(bondlist, couplings, tj_nosym, tj_nosym);
        REQUIRE(lila::close(H_nosym, lila::Herm(H_nosym)));
        auto eigs_nosym = lila::EigenvaluesSym(H_nosym);

        lila::Vector<double> eigs_sym;
        for (int k = 0; k < (int)irreps.size(); ++k) {
          auto irrep = irreps[k];
          int multiplicity = multiplicities[k];

          auto tj = tJSymmetric<bit_t>(n_sites, nup, ndn, space_group, irrep);

          if (tj.size() > 0) {

            // Compute partial spectrum from symmetrized block
            auto H_sym = MatrixCplx(bondlist, couplings, tj, tj);
            REQUIRE(lila::close(H_sym, lila::Herm(H_sym)));
            auto eigs_sym_k = lila::EigenvaluesSym(H_sym);

	    // Check whether results are the same for real blocks
	    if (!is_complex(tj.irrep())) {
	      auto H_sym_real = MatrixReal(bondlist, couplings, tj, tj);
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
        REQUIRE(lila::close(eigs_sym, eigs_nosym));
      }
    }
  }
}

template <class bit_t> void test_tj_symmetric_spectrum_chains(int n_sites) {
  using namespace hydra::testcases::tj;
  using namespace hydra::testcases::electron;

  lila::Log.out("Tj chain, symmetric spectra test, n_sites: {}", n_sites);
  auto [bondlist, couplings] = tJchain(n_sites, 1.0, 5.0);
  auto [space_group, irreps, multiplicities] =
      get_cyclic_group_irreps_mult(n_sites);
  test_symmetric_spectra<uint32>(bondlist, couplings, space_group, irreps,
                                 multiplicities);
}

TEST_CASE("tj_symmetric_matrix", "[blocks][tj_symmetric]") {
  using namespace hydra::testcases::tj;
  using namespace hydra::testcases::electron;

  // Test linear chains
  for (int n_sites = 2; n_sites < 7; ++n_sites) {
    test_tj_symmetric_spectrum_chains<uint16_t>(n_sites);
    test_tj_symmetric_spectrum_chains<uint32_t>(n_sites);
    test_tj_symmetric_spectrum_chains<uint64_t>(n_sites);
  }

  // test a 3x3 triangular lattice
  lila::Log.out("Tj 3x3 triangular, symmetric spectra test");
  std::string lfile = "data/triangular.9.hop.sublattices.tsl.lat";

  auto bondlist = read_bondlist(lfile);
  Couplings couplings;
  couplings["T"] = 1.0;
  couplings["U"] = 5.0;
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
    irreps.push_back(read_represenation(lfile, name));
    multiplicities.push_back(mult);
  }
  test_symmetric_spectra<uint32>(bondlist, couplings, space_group, irreps,
                                 multiplicities);
}
