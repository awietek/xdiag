#include "../catch.hpp"

#include <hydra/all.h>

using namespace hydra;
using namespace hydra::combinatorics;
using namespace hydra::symmetries;

static PermutationGroup cyclic_group(int n_sites) {
  // test cyclic group
  std::vector<int> permutation_array;
  for (int sym = 0; sym < n_sites; ++sym) {
    for (int site = 0; site < n_sites; ++site) {
      int newsite = (site + sym) % n_sites;
      permutation_array.push_back(newsite);
    }
  }
  return PermutationGroup(n_sites, n_sites, permutation_array);
}

template <typename bit_t> void test_stabilizer_symmetries(int n_sites) {
  auto group_action = PermutationGroupLookup<bit_t>(cyclic_group(n_sites));

  for (bit_t bits : Subsets(n_sites)) {
    auto stab_syms = stabilizer_symmetries(bits, group_action);
    for (int sym : stab_syms) {
      REQUIRE(group_action.apply(sym, bits) == bits);
    }
  }
}

template <typename bit_t> void test_representative(int n_sites) {
  auto group_action = PermutationGroupLookup<bit_t>(cyclic_group(n_sites));

  for (bit_t bits : Subsets(n_sites)) {
    bit_t rep = representative(bits, group_action);
    int hit_ctr = 0;
    int n_sym = group_action.n_symmetries();
    for (int sym = 0; sym < n_sym; ++sym) {
      REQUIRE(group_action.apply(sym, bits) >= rep);
      if (group_action.apply(sym, bits) == rep)
        hit_ctr++;
    }
    if (n_sym)
      REQUIRE(hit_ctr > 0);
  }
}

template <typename bit_t> void test_representative_subset(int n_sites) {
  auto group_action = PermutationGroupLookup<bit_t>(cyclic_group(n_sites));
  std::vector<int> subset;
  for (int i = 0; i < n_sites; i += 2)
    subset.push_back(i);

  for (bit_t bits : Subsets(n_sites)) {
    bit_t rep = representative_subset(bits, group_action, subset);
    int hit_ctr = 0;
    int n_sym = subset.size();
    for (int sym : subset) {
      REQUIRE(group_action.apply(sym, bits) >= rep);
      if (group_action.apply(sym, bits) == rep)
        hit_ctr++;
    }
    if (n_sym > 0)
      REQUIRE(hit_ctr > 0);
  }
}

template <typename bit_t> void test_representative_sym(int n_sites) {
  auto group_action = PermutationGroupLookup<bit_t>(cyclic_group(n_sites));

  for (bit_t bits : Subsets(n_sites)) {
    auto [rep, rsym] = representative_sym(bits, group_action);
    int hit_ctr = 0;
    int n_sym = group_action.n_symmetries();
    for (int sym = 0; sym < n_sym; ++sym) {
      REQUIRE(group_action.apply(sym, bits) >= rep);
      if (group_action.apply(sym, bits) == rep)
        hit_ctr++;
    }
    if (n_sym > 0) {
      REQUIRE(hit_ctr > 0);
      REQUIRE(rsym < n_sym);
      REQUIRE(group_action.apply(rsym, bits) == rep);
    }
  }
}

template <typename bit_t> void test_representative_sym_subset(int n_sites) {
  auto group_action = PermutationGroupLookup<bit_t>(cyclic_group(n_sites));
  std::vector<int> subset;
  for (int i = 0; i < n_sites; i += 2)
    subset.push_back(i);

  for (bit_t bits : Subsets(n_sites)) {
    auto [rep, rsym] = representative_sym_subset(bits, group_action, subset);
    int hit_ctr = 0;
    int n_sym = group_action.n_symmetries();
    for (int sym : subset) {
      REQUIRE(group_action.apply(sym, bits) >= rep);
      if (group_action.apply(sym, bits) == rep)
        hit_ctr++;
    }
    if (n_sym > 0) {
      REQUIRE(hit_ctr > 0);
      REQUIRE(group_action.apply(rsym, bits) == rep);
    }
  }
}

TEST_CASE("symmetry_operations", "[symmetries]") {

  for (int n_sites = 0; n_sites <= 5; ++n_sites) {
    test_stabilizer_symmetries<uint16_t>(n_sites);
    test_stabilizer_symmetries<uint32_t>(n_sites);
    test_stabilizer_symmetries<uint64_t>(n_sites);
  }

  for (int n_sites = 0; n_sites <= 5; ++n_sites) {
    test_representative<uint16_t>(n_sites);
    test_representative<uint32_t>(n_sites);
    test_representative<uint64_t>(n_sites);
  }

  for (int n_sites = 0; n_sites <= 5; ++n_sites) {
    test_representative_subset<uint16_t>(n_sites);
    test_representative_subset<uint32_t>(n_sites);
    test_representative_subset<uint64_t>(n_sites);
  }

  for (int n_sites = 0; n_sites <= 5; ++n_sites) {
    test_representative_sym<uint16_t>(n_sites);
    test_representative_sym<uint32_t>(n_sites);
    test_representative_sym<uint64_t>(n_sites);
  }

  for (int n_sites = 0; n_sites <= 5; ++n_sites) {
    test_representative_sym_subset<uint16_t>(n_sites);
    test_representative_sym_subset<uint32_t>(n_sites);
    test_representative_sym_subset<uint64_t>(n_sites);
  }
}
