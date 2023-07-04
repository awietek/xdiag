#include "../../catch.hpp"

#include <iostream>

#include <hydra/common.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/subsets.h>
#include <hydra/indexing/spinhalf/indexing_sublattice.h>
#include <hydra/indexing/spinhalf/indexing_symmetric_sz.h>
#include <hydra/indexing/spinhalf/indexing_symmetric_no_sz.h>
#include <hydra/symmetries/operations/symmetry_operations.h>

using namespace hydra;

template <typename bit_t, class Indexing1, class Indexing2>
void compare_state(bit_t state, Indexing1 const &indexing1,
                   Indexing2 const &indexing2) {
  auto ga1 = indexing1.group_action();
  auto ga2 = indexing2.group_action();
  idx_t idx1 = indexing1.index(state);
  idx_t idx2 = indexing2.index(state);
  if (idx1 != invalid_index) {
    REQUIRE(idx1 == idx2);
    bit_t rep1 = indexing1.representative(state);
    bit_t rep2 = indexing2.representative(state);
    REQUIRE(rep1 == rep2);
    bit_t norm1 = indexing1.norm(idx1);
    bit_t norm2 = indexing2.norm(idx2);
    REQUIRE(norm1 == norm2);
    auto [idx1a, sym1] = indexing1.index_sym(state);
    auto [idx2a, sym2] = indexing2.index_sym(state);
    REQUIRE(idx1a == idx1);
    REQUIRE(idx2a == idx2);
    REQUIRE(ga1.apply(sym1, state) == rep1);
    REQUIRE(ga2.apply(sym2, state) == rep2);
    auto [idx1b, syms1] = indexing1.index_syms(state);
    auto [idx2b, syms2] = indexing2.index_syms(state);
    REQUIRE(idx1b == idx1);
    REQUIRE(idx2b == idx2);
    REQUIRE(syms1.size() > 0);
    REQUIRE(syms2.size() > 0);
    for (auto s : syms1) {
      REQUIRE(ga1.apply(s, state) == rep1);
    }
    for (auto s : syms2) {
      REQUIRE(ga2.apply(s, state) == rep2);
    }
  }
}

template <typename bit_t, class Indexing1, class Indexing2>
void compare_indices_sz(Indexing1 const &indexing1, Indexing2 const &indexing2,
                        int nup) {
  REQUIRE(indexing1.size() == indexing2.size());
  for (idx_t idx = 0; idx < indexing1.size(); ++idx) {
    REQUIRE(indexing1.state(idx) == indexing2.state(idx));
  }

  int n_sites = indexing1.n_sites();
  REQUIRE(n_sites == indexing2.n_sites());
  for (auto state : combinatorics::Combinations<bit_t>(n_sites, nup)) {
    compare_state(state, indexing1, indexing2);
  }
}

template <typename bit_t, class Indexing1, class Indexing2>
void compare_indices_no_sz(Indexing1 const &indexing1,
                           Indexing2 const &indexing2) {
  REQUIRE(indexing1.size() == indexing2.size());
  for (idx_t idx = 0; idx < indexing1.size(); ++idx) {
    REQUIRE(indexing1.state(idx) == indexing2.state(idx));
  }

  int n_sites = indexing1.n_sites();
  REQUIRE(n_sites == indexing2.n_sites());
  for (auto state : combinatorics::Subsets<bit_t>(n_sites)) {
    compare_state(state, indexing1, indexing2);
  }
}

template <class bit_t> void test_spinhalf_indexing_sublattice() {
  using indexing_no_sz_t = indexing::spinhalf::IndexingSymmetricNoSz<bit_t>;
  using indexing_sz_t = indexing::spinhalf::IndexingSymmetricSz<bit_t>;
  using indexing_sl1_t = indexing::spinhalf::IndexingSublattice<bit_t, 1>;
  using indexing_sl2_t = indexing::spinhalf::IndexingSublattice<bit_t, 2>;
  using indexing_sl3_t = indexing::spinhalf::IndexingSublattice<bit_t, 3>;
  using indexing_sl4_t = indexing::spinhalf::IndexingSublattice<bit_t, 4>;
  using indexing_sl5_t = indexing::spinhalf::IndexingSublattice<bit_t, 5>;

  {
    Log("SpinhalfIndexingSublattice: 1 sublattice");
    int n_sites = 8;
    std::string lfile =
        HYDRA_DIRECTORY "/misc/data/square.8.heisenberg.2sl.lat";
    auto permutations = hydra::read_permutations(lfile);
    auto perm_group = PermutationGroup(permutations);
    std::vector<std::string> irrep_names = {
        "Gamma.D4.A1", "Gamma.D4.A2", "Gamma.D4.B1", "Gamma.D4.B2",
        "Gamma.D4.E",  "M.D4.A1",     "M.D4.A2",     "M.D4.B1",
        "M.D4.B2",     "M.D4.E",      "Sigma.D1.A",  "Sigma.D1.B",
        "X.D2.A1",     "X.D2.A2",     "X.D2.B1",     "X.D2.B2"};

    for (auto irrep_name : irrep_names) {
      auto irrep = read_representation(lfile, irrep_name);
      auto idxng = indexing_no_sz_t(n_sites, perm_group, irrep);
      auto idxng_sl = indexing_sl1_t(n_sites, perm_group, irrep);
      compare_indices_no_sz<bit_t, indexing_no_sz_t, indexing_sl1_t>(idxng,
                                                                     idxng_sl);
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto idxng = indexing_sz_t(n_sites, nup, perm_group, irrep);
        auto idxng_sl = indexing_sl1_t(n_sites, nup, perm_group, irrep);
        compare_indices_sz<bit_t, indexing_sz_t, indexing_sl1_t>(idxng,
                                                                 idxng_sl, nup);
      }
    }
  }

  {
    Log("SpinhalfIndexingSublattice: 2 sublattice");
    int n_sites = 8;
    std::string lfile =
        HYDRA_DIRECTORY "/misc/data/square.8.heisenberg.2sl.lat";
    auto permutations = hydra::read_permutations(lfile);
    auto perm_group = PermutationGroup(permutations);
    std::vector<std::string> irrep_names = {
        "Gamma.D4.A1", "Gamma.D4.A2", "Gamma.D4.B1", "Gamma.D4.B2",
        "Gamma.D4.E",  "M.D4.A1",     "M.D4.A2",     "M.D4.B1",
        "M.D4.B2",     "M.D4.E",      "Sigma.D1.A",  "Sigma.D1.B",
        "X.D2.A1",     "X.D2.A2",     "X.D2.B1",     "X.D2.B2"};

    for (auto irrep_name : irrep_names) {
      auto irrep = read_representation(lfile, irrep_name);
      auto idxng = indexing_no_sz_t(n_sites, perm_group, irrep);
      auto idxng_sl = indexing_sl2_t(n_sites, perm_group, irrep);
      compare_indices_no_sz<bit_t, indexing_no_sz_t, indexing_sl2_t>(idxng,
                                                                     idxng_sl);
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto idxng = indexing_sz_t(n_sites, nup, perm_group, irrep);
        auto idxng_sl = indexing_sl2_t(n_sites, nup, perm_group, irrep);
        compare_indices_sz<bit_t, indexing_sz_t, indexing_sl2_t>(idxng,
                                                                 idxng_sl, nup);
      }
    }
  }

  {
    Log("SpinhalfIndexingSublattice: 3 sublattice");
    int n_sites = 9;
    std::string lfile =
        HYDRA_DIRECTORY "/misc/data/square.9.heisenberg.3sl.lat";
    auto permutations = hydra::read_permutations(lfile);
    auto perm_group = PermutationGroup(permutations);
    std::vector<std::string> irrep_names = {
        "Gamma.D2.A1", "Gamma.D2.A2", "Gamma.D2.B1",
        "Gamma.D2.B2", "Delta.C1.A",  "Sigma0.D1.A",
        "Sigma0.D1.B", "Sigma1.D1.A", "Sigma1.D1.B"};

    for (auto irrep_name : irrep_names) {
      auto irrep = read_representation(lfile, irrep_name);
      auto idxng = indexing_no_sz_t(n_sites, perm_group, irrep);
      auto idxng_sl = indexing_sl3_t(n_sites, perm_group, irrep);
      compare_indices_no_sz<bit_t, indexing_no_sz_t, indexing_sl3_t>(idxng,
                                                                     idxng_sl);
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto idxng = indexing_sz_t(n_sites, nup, perm_group, irrep);
        auto idxng_sl = indexing_sl3_t(n_sites, nup, perm_group, irrep);
        compare_indices_sz<bit_t, indexing_sz_t, indexing_sl3_t>(idxng,
                                                                 idxng_sl, nup);
      }
    }
  }

  {
    Log("SpinhalfIndexingSublattice: 3 sublattice (triangular)");
    int n_sites = 9;
    std::string lfile = HYDRA_DIRECTORY
        "/misc/data/triangular.9.Jz1Jz2Jx1Jx2D1.sublattices.tsl.lat";
    auto permutations = hydra::read_permutations(lfile);
    auto perm_group = PermutationGroup(permutations);
    std::vector<std::string> irrep_names = {
        "Gamma.D6.A1", "Gamma.D6.A2", "Gamma.D6.B1", "Gamma.D6.B2",
        "Gamma.D6.E1", "Gamma.D6.E2", "K.D3.A1",     "K.D3.A2",
        "K.D3.E",      "Y.D1.A",      "Y.D1.B"};

    for (auto irrep_name : irrep_names) {
      auto irrep = read_representation(lfile, irrep_name);
      auto idxng = indexing_no_sz_t(n_sites, perm_group, irrep);
      auto idxng_sl = indexing_sl3_t(n_sites, perm_group, irrep);
      compare_indices_no_sz<bit_t, indexing_no_sz_t, indexing_sl3_t>(idxng,
                                                                     idxng_sl);
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto idxng = indexing_sz_t(n_sites, nup, perm_group, irrep);
        auto idxng_sl = indexing_sl3_t(n_sites, nup, perm_group, irrep);
        compare_indices_sz<bit_t, indexing_sz_t, indexing_sl3_t>(idxng,
                                                                 idxng_sl, nup);
      }
    }
  }

  {
    Log("SpinhalfIndexingSublattice: 4 sublattice");
    int n_sites = 8;
    std::string lfile =
        HYDRA_DIRECTORY "/misc/data/square.8.heisenberg.4sl.lat";
    auto permutations = hydra::read_permutations(lfile);
    auto perm_group = PermutationGroup(permutations);
    std::vector<std::string> irrep_names = {
        "Gamma.D4.A1", "Gamma.D4.A2", "Gamma.D4.B1", "Gamma.D4.B2",
        "Gamma.D4.E",  "M.D4.A1",     "M.D4.A2",     "M.D4.B1",
        "M.D4.B2",     "M.D4.E",      "Sigma.D1.A",  "Sigma.D1.B",
        "X.D2.A1",     "X.D2.A2",     "X.D2.B1",     "X.D2.B2"};

    for (auto irrep_name : irrep_names) {
      auto irrep = read_representation(lfile, irrep_name);
      auto idxng = indexing_no_sz_t(n_sites, perm_group, irrep);
      auto idxng_sl = indexing_sl4_t(n_sites, perm_group, irrep);
      compare_indices_no_sz<bit_t, indexing_no_sz_t, indexing_sl4_t>(idxng,
                                                                     idxng_sl);
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto idxng = indexing_sz_t(n_sites, nup, perm_group, irrep);
        auto idxng_sl = indexing_sl4_t(n_sites, nup, perm_group, irrep);
        compare_indices_sz<bit_t, indexing_sz_t, indexing_sl4_t>(idxng,
                                                                 idxng_sl, nup);
      }
    }
  }

  {
    Log("SpinhalfIndexingSublattice: 5 sublattice");
    int n_sites = 10;
    std::string lfile =
        HYDRA_DIRECTORY "/misc/data/square.10.heisenberg.5sl.lat";
    auto permutations = hydra::read_permutations(lfile);
    auto perm_group = PermutationGroup(permutations);
    std::vector<std::string> irrep_names = {
        "Gamma.C2.A", "Gamma.C2.B", "Delta0.C1.A", "Delta1.C1.A", "X.C2.A",
        "X.C2.B",     "Z0.C1.A",    "Z1.C1.A",     "Z2.C1.A",     "Z3.C1.A"};

    for (auto irrep_name : irrep_names) {
      auto irrep = read_representation(lfile, irrep_name);
      auto idxng = indexing_no_sz_t(n_sites, perm_group, irrep);
      auto idxng_sl = indexing_sl5_t(n_sites, perm_group, irrep);
      compare_indices_no_sz<bit_t, indexing_no_sz_t, indexing_sl5_t>(idxng,
                                                                     idxng_sl);
      for (int nup = 0; nup <= n_sites; ++nup) {
        auto idxng = indexing_sz_t(n_sites, nup, perm_group, irrep);
        auto idxng_sl = indexing_sl5_t(n_sites, nup, perm_group, irrep);
        compare_indices_sz<bit_t, indexing_sz_t, indexing_sl5_t>(idxng,
                                                                 idxng_sl, nup);
      }
    }
  }
}

TEST_CASE("SpinhalfIndexingSublattice", "[symmetries]") {
  Log("Test SpinhalfIndexingSublattice");
  Log("uint16_t");
  test_spinhalf_indexing_sublattice<uint16_t>();
  Log("uint32_t");
  test_spinhalf_indexing_sublattice<uint32_t>();
  Log("uint64_t");
  test_spinhalf_indexing_sublattice<uint64_t>();
  Log("Done");
}
