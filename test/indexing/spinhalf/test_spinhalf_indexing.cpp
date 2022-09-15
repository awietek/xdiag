#include "../../catch.hpp"

#include <iostream>

#include <hydra/all.h>

using namespace hydra;

template <typename bit_t, class Indexing>
void test_spinhalf_indexing_state(Indexing const &indexing, bit_t state) {

  auto const &group_action = indexing.group_action();
  auto const &irrep = indexing.irrep();

  idx_t idx_rep = indexing.index(state);
  if (idx_rep == invalid_index) {
    double norm = symmetries::norm(state, group_action, irrep);
    REQUIRE(std::abs(norm) < 1e-6);
  } else {

    bit_t rep = symmetries::representative(state, group_action);
    bit_t rep_lookup = indexing.representative(state);
    REQUIRE(rep == rep_lookup);

    REQUIRE(rep == indexing.state(idx_rep));

    // int n_sites = indexing.n_sites();
    // std::cout << BSTR(state) << "\n";

    double norm = symmetries::norm(state, group_action, irrep);
    REQUIRE(norm == indexing.norm(idx_rep));

    auto [idxr1, sym] = indexing.index_sym(state);
    (void)idxr1;
    REQUIRE(group_action.apply(sym, state) == rep);

    auto [idxr2, syms] = indexing.index_syms(state);
    (void)idxr2;
    for (int sym : syms) {
      REQUIRE(group_action.apply(sym, state) == rep);
    }

  }
}

template <typename bit_t>
void check_indexing_symmetric_sz(
    indexing::spinhalf::IndexingSymmetricSz<bit_t> const &indexing) {

  int n_sites = indexing.n_sites();
  int n_up = indexing.n_up();

  if (indexing.size() > 0) {
    for (auto state : combinatorics::Combinations(n_sites, n_up)) {
      test_spinhalf_indexing_state(indexing, state);
    }
  }
}

template <typename bit_t>
void check_indexing_symmetric_no_sz(
    indexing::spinhalf::IndexingSymmetricNoSz<bit_t> const &indexing) {

  int n_sites = indexing.n_sites();
  if (indexing.size() > 0) {

    for (auto state : combinatorics::Subsets(n_sites)) {
      test_spinhalf_indexing_state(indexing, state);
    }
  }
}

template <class bit_t> void test_indexing_symmetric() {

  using namespace hydra::indexing::spinhalf;

  Log("IndexingSymmetric: triangular 3x3");
  int n_sites = 9;
  std::string lfile = "data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.lat";
  auto permutations = hydra::read_permutations(lfile);
  auto perm_group = PermutationGroup(permutations);
  // std::vector<std::string> irrep_names = {
  //     "Gamma.D6.A1", "Gamma.D6.A2", "Gamma.D6.B1", "Gamma.D6.B2",
  //     "Gamma.D6.E1", "Gamma.D6.E2", "K.D3.A1",     "K.D3.A2",
  //     "K.D3.E",      "Y.D1.A",      "Y.D1.B"};

  std::vector<std::string> irrep_names = {"Gamma.D6.B2"};
  
  for (auto irrep_name : irrep_names) {
    Log("irrep: {}", irrep_name);
    auto irrep = read_represenation(lfile, irrep_name);
    auto idxng = IndexingSymmetricNoSz<bit_t>(n_sites, perm_group, irrep);
    // check_indexing_symmetric_no_sz<bit_t>(idxng);

    for (int nup = 3; nup <= 3; ++nup) {
      auto idxng = IndexingSymmetricSz<bit_t>(n_sites, nup, perm_group, irrep);
      check_indexing_symmetric_sz<bit_t>(idxng);
    }
  }
}

TEST_CASE("spinhalf_indexing", "[indexing]") {
  Log("Test SpinhalfIndexingSublattice");
  // Log("uint16_t");
  // test_indexing_symmetric<uint16_t>();
  // Log("uint32_t");
  // test_indexing_symmetric<uint32_t>();
  // Log("uint64_t");
  test_indexing_symmetric<uint64_t>();
  Log("Done");
}
