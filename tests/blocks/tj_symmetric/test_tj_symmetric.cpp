#include "../../catch.hpp"

#include <iostream>

#include <xdiag/blocks/tj/tj.hpp>
#include <xdiag/utils/print_macro.hpp>

using namespace xdiag;

void check_dimensions_sum_up_tj_symmetric(int64_t n_sites,
                                          PermutationGroup group,
                                          std::vector<Representation> irreps) {
  using combinatorics::binomial;

  Log.out("tj_symmetric: dimension sum test. N: {}", n_sites);

  int64_t sum_of_dims = 0;

  for (int64_t nup = 0; nup <= n_sites; ++nup) {
    for (int64_t ndn = 0; ndn <= n_sites - nup; ++ndn) {
      int64_t sum_of_dims_updn = 0;
      for (auto irrep : irreps) {

        auto block = tJ(n_sites, nup, ndn, group, irrep);
        sum_of_dims += block.size();
        sum_of_dims_updn += block.size();
      }
      REQUIRE(sum_of_dims_updn ==
              binomial(n_sites, nup) * binomial(n_sites - nup, ndn));
    }
  }
  int64_t p = 1;
  for (int64_t i = 0; i < n_sites; ++i) {
    p *= 3;
  }
  REQUIRE(sum_of_dims == p);
}

TEST_CASE("tj_symmetric", "[tj]") {

  // Test linear chains
  Log("tj_symmetric: chain test");
  for (int64_t n_sites = 1; n_sites < 7; ++n_sites) {

    // test cyclic group
    std::vector<Permutation> permutation_array;
    for (int64_t sym = 0; sym < n_sites; ++sym) {

      std::vector<int64_t> pv;
      for (int64_t site = 0; site < n_sites; ++site) {
        int64_t newsite = (site + sym) % n_sites;
        pv.push_back(newsite);
      }
      permutation_array.push_back(Permutation(pv));
    }
    auto group = PermutationGroup(permutation_array);

    // Create irreps
    std::vector<Representation> irreps;
    for (int64_t k = 0; k < n_sites; ++k) {
      std::vector<complex> chis;
      for (int64_t l = 0; l < n_sites; ++l)
        chis.push_back({std::cos(2 * M_PI * l * k / n_sites),
                        std::sin(2 * M_PI * l * k / n_sites)});
      irreps.push_back(Representation(chis));
    }

    check_dimensions_sum_up_tj_symmetric(n_sites, group, irreps);
  }

  // test a 3x3 triangular lattice
  int64_t n_sites = 9;
  Log("tj_symmetric: triangular 3x3 test");

  std::string lfile =
      XDIAG_DIRECTORY "/misc/data/triangular.9.hop.sublattices.tsl.lat";
  auto permutations = xdiag::read_permutations(lfile);
  auto group = PermutationGroup(permutations);

  std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
      {"Gamma.D3.A1", 1}, {"Gamma.D3.A2", 1}, {"Gamma.D3.E", 2},
      {"K0.D3.A1", 1},    {"K0.D3.A2", 1},    {"K0.D3.E", 2},
      {"K1.D3.A1", 1},    {"K1.D3.A2", 1},    {"K1.D3.E", 2},
      {"Y.C1.A", 6}};

  std::vector<Representation> irreps;
  for (auto [name, mult] : rep_name_mult) {
    for (int64_t i = 0; i < mult; ++i) {
      irreps.push_back(read_representation(lfile, name));
    }
  }

  check_dimensions_sum_up_tj_symmetric(n_sites, group, irreps);
}
