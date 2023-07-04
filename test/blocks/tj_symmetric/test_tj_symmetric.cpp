#include "../../catch.hpp"

#include <iostream>

#include <hydra/blocks/tj/tj.h>

using namespace hydra;

void check_dimensions_sum_up_tj_symmetric(int n_sites, PermutationGroup group,
                                          std::vector<Representation> irreps) {
  using combinatorics::binomial;

  Log.out("tj_symmetric: dimension sum test. N: {}", n_sites);

  idx_t sum_of_dims = 0;

  for (int nup = 0; nup <= n_sites; ++nup) {
    for (int ndn = 0; ndn <= n_sites - nup; ++ndn) {
      idx_t sum_of_dims_updn = 0;
      for (auto irrep : irreps) {

        auto block = tJ(n_sites, nup, ndn, group, irrep);
        sum_of_dims += block.size();
        sum_of_dims_updn += block.size();
      }
      REQUIRE(sum_of_dims_updn ==
              binomial(n_sites, nup) * binomial(n_sites - nup, ndn));
    }
  }
  REQUIRE(sum_of_dims == idx_t(pow(3, n_sites)));
}

TEST_CASE("tj_symmetric", "[tj]") {

  // Test linear chains
  Log("tj_symmetric: chain test");
  for (int n_sites = 1; n_sites < 7; ++n_sites) {

    // test cyclic group
    std::vector<Permutation> permutation_array;
    for (int sym = 0; sym < n_sites; ++sym) {

      std::vector<int> pv;
      for (int site = 0; site < n_sites; ++site) {
        int newsite = (site + sym) % n_sites;
        pv.push_back(newsite);
      }
      permutation_array.push_back(Permutation(pv));
    }
    auto group = PermutationGroup(permutation_array);

    // Create irreps
    std::vector<Representation> irreps;
    for (int k = 0; k < n_sites; ++k) {
      std::vector<complex> chis;
      for (int l = 0; l < n_sites; ++l)
        chis.push_back({std::cos(2 * M_PI * l * k / n_sites),
                        std::sin(2 * M_PI * l * k / n_sites)});
      irreps.push_back(Representation(chis));
    }

    check_dimensions_sum_up_tj_symmetric(n_sites, group, irreps);
  }

  // test a 3x3 triangular lattice
  int n_sites = 9;
  Log("tj_symmetric: triangular 3x3 test");

  std::string lfile =
      HYDRA_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.lat";
  auto permutations = hydra::read_permutations(lfile);
  auto group = PermutationGroup(permutations);

  std::vector<std::pair<std::string, int>> rep_name_mult = {
      {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
      {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
      {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
      {"Y.C1.A", 6}};

  std::vector<Representation> irreps;
  for (auto [name, mult] : rep_name_mult) {
    for (int i = 0; i < mult; ++i) {
      irreps.push_back(read_representation(lfile, name));
    }
  }

  check_dimensions_sum_up_tj_symmetric(n_sites, group, irreps);
}
