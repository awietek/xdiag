#include "../../catch.hpp"

#include <iostream>

#include <xdiag/blocks/spinhalf/spinhalf_apply.hpp>
#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/matrix.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/utils/close.hpp>
#include <xdiag/symmetries/group_action/group_action.hpp>

using namespace xdiag;
using namespace xdiag::combinatorics;

void test_indices_spinhalf_symmetric(Spinhalf const &block) {

  GroupAction group_action(block.permutation_group());
  // auto indexing = block.indexing_sym_sz_conserved();
  // for (int64_t idx = 0; idx < block.size(); ++idx) {
  //   bit_t state = indexing.state(idx);
  //   bit_t rep = symmetries::representative(state, group_action);
  //   REQUIRE(rep == state);
  //   int64_t idx2 = indexing.index(state);
  //   // auto norm = block.indexing_.norm(state);
  //   REQUIRE(idx == idx2);
  // }
}

void test_spinchain_blocks(int64_t n_sites) {

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
  auto space_group = PermutationGroup(permutation_array);

  // Create irrep with K momentum
  int64_t sum_of_dims = 0;
  for (int64_t nup = 0; nup <= n_sites; ++nup) {
    int64_t sum_of_dims_up = 0;

    for (int64_t k = 0; k < n_sites; ++k) {

      // Create characters
      std::vector<complex> chis;
      for (int64_t l = 0; l < n_sites; ++l)
        chis.push_back({std::cos(2 * M_PI * l * k / n_sites),
                        std::sin(2 * M_PI * l * k / n_sites)});
      auto irrep = Representation(chis);

      auto block = Spinhalf(n_sites, nup, space_group, irrep);

      sum_of_dims += block.size();
      sum_of_dims_up += block.size();

      if (block.size() > 0) {
        test_indices_spinhalf_symmetric(block);
      }
    }
    REQUIRE(sum_of_dims_up == binomial(n_sites, nup));
  }
  REQUIRE(sum_of_dims == pow(2, n_sites));
}

TEST_CASE("spinhalf_symmetric", "[spinhalf]") {

  // Test the tJ chain
  for (int64_t n_sites = 1; n_sites < 8; ++n_sites) {
    Log.out("spinhalf_symmetric: block test: Spinhalf Chain {}", n_sites);
    test_spinchain_blocks(n_sites);
  }

  // test a 3x3 triangular lattice
  Log.out("spinhalf_symmetric: block test: Triangular 3x3");
  int64_t n_sites = 9;

  std::vector<std::pair<std::string, int64_t>> rep_name_mult = {
      {"Gamma.D6.A1", 1}, {"Gamma.D6.A2", 1}, {"Gamma.D6.B1", 1},
      {"Gamma.D6.B2", 1}, {"Gamma.D6.E1", 2}, {"Gamma.D6.E2", 2},
      {"K.D3.A1", 2},     {"K.D3.A2", 2},     {"K.D3.E", 4},
      {"Y.D1.A", 6},      {"Y.D1.B", 6}};

  std::string lfile = XDIAG_DIRECTORY
      "/misc/data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.lat";
  auto permutations = xdiag::read_permutations(lfile);
  auto space_group = PermutationGroup(permutations);

  int64_t sum_dim = 0;
  for (int64_t nup = 0; nup <= n_sites; ++nup) {
    int64_t sum_dim_up = 0;

    for (auto [name, mult] : rep_name_mult) {

      auto irrep = read_representation(lfile, name);
      auto block = Spinhalf(n_sites, nup, space_group, irrep);

      int64_t dim = block.size() * mult;
      test_indices_spinhalf_symmetric(block);

      sum_dim_up += dim;
      sum_dim += dim;
    }
    REQUIRE(sum_dim_up == binomial(n_sites, nup));
  }
  REQUIRE(sum_dim == (int64_t)pow(2, n_sites));
}
