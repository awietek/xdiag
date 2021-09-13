#include "../../catch.hpp"

#include <iostream>

#include <hydra/all.h>

using namespace hydra;
using namespace hydra::combinatorics;

template <class bit_t>
void test_indices_spinhalf_symmetric(SpinhalfSymmetric<bit_t> const &block) {
  
  PermutationGroupAction group_action(block.permutation_group());
  for (idx_t idx = 0; idx < block.size(); ++idx) {
    bit_t state = block.indexing_.state(idx);
    bit_t rep = utils::representative(state, group_action);
    REQUIRE(rep == state);
    idx_t idx2 = block.indexing_.index(state);
    // auto norm = block.indexing_.norm(state);
    // lila::Log.out("{} {} {} ({},{})", idx, state, idx2, lila::real(norm),
    //               lila::imag(norm));
    REQUIRE(idx == idx2);
  }
}

template <class bit_t> void test_spinchain_blocks(int n_sites) {

  // cyclic group as space group
  std::vector<std::vector<int>> permutations;
  for (int sym = 0; sym < n_sites; ++sym) {
    std::vector<int> permutation;
    for (int site = 0; site < n_sites; ++site) {
      int newsite = (site + sym) % n_sites;
      permutation.push_back(newsite);
    }
    permutations.push_back(permutation);
  }
  auto space_group = PermutationGroup(permutations);

  // Create irrep with K momentum
  idx_t sum_of_dims = 0;
  for (int nup = 0; nup <= n_sites; ++nup) {
    idx_t sum_of_dims_up = 0;

    for (int k = 0; k < n_sites; ++k) {

      // Create characters
      std::vector<complex> chis;
      for (int l = 0; l < n_sites; ++l)
        chis.push_back({std::cos(2 * M_PI * l * k / n_sites),
                        std::sin(2 * M_PI * l * k / n_sites)});
      auto irrep = Representation(chis);

      auto block = SpinhalfSymmetric<bit_t>(n_sites, nup, space_group, irrep);

      sum_of_dims += block.size();
      sum_of_dims_up += block.size();

      if (block.size() > 0) {
        test_indices_spinhalf_symmetric(block);
      }
      lila::Log.out("{} {} {} {}", n_sites, nup, k, block.size());

    }
    REQUIRE(sum_of_dims_up == binomial(n_sites, nup));
  }
  REQUIRE(sum_of_dims == pow(2, n_sites));
}

TEST_CASE("SpinhalfSymmetric", "[models]") {

  // Test the tJ chain
  for (int n_sites = 1; n_sites < 8; ++n_sites) {
    lila::Log.out("SpinhalfSymmetric block test: Spinhalf Chain {}", n_sites);

    test_spinchain_blocks<hydra::uint16>(n_sites);
    test_spinchain_blocks<hydra::uint32>(n_sites);
    test_spinchain_blocks<hydra::uint64>(n_sites);
  }

  // test a 3x3 triangular lattice
  lila::Log.out("SpinhalfSymmetric block test: Triangular 3x3");
  int n_sites = 9;

  std::vector<std::pair<std::string, int>> rep_name_mult = {
      {"Gamma.D6.A1", 1}, {"Gamma.D6.A2", 1}, {"Gamma.D6.B1", 1},
      {"Gamma.D6.B2", 1}, {"Gamma.D6.E1", 2}, {"Gamma.D6.E2", 2},
      {"K.D3.A1", 2},     {"K.D3.A2", 2},     {"K.D3.E", 4},
      {"Y.D1.A", 6},      {"Y.D1.B", 6}};

  std::string lfile = "data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.lat";
  auto permutations = hydra::read_permutations(lfile);
  auto space_group = PermutationGroup(permutations);

  idx_t sum_dim = 0;
  for (int nup = 0; nup <= n_sites; ++nup) {
    idx_t sum_dim_up = 0;

    for (auto [name, mult] : rep_name_mult) {
      auto irrep = read_represenation(lfile, name);
      auto block = SpinhalfSymmetric<uint16>(n_sites, nup, space_group, irrep);

      idx_t dim = block.size() * mult;
      test_indices_spinhalf_symmetric(block);

      sum_dim_up += dim;
      sum_dim += dim;
    }
    REQUIRE(sum_dim_up == binomial(n_sites, nup));
  }
  REQUIRE(sum_dim == (idx_t)pow(2, n_sites));
}
