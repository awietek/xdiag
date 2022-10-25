#include "../../catch.hpp"

#include <hydra/all.h>

using namespace hydra;
using namespace hydra::combinatorics;
using namespace hydra::symmetries;

static PermutationGroup cyclic_group(int n_sites) {
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
  return PermutationGroup(permutation_array);
}

template <typename bit_t> void test_stabilizer_symmetries(int n_sites) {
  auto group_action = GroupActionLookup<bit_t>(cyclic_group(n_sites));

  for (bit_t bits : Subsets(n_sites)) {
    auto stab_syms = stabilizer_symmetries(bits, group_action);
    for (int sym : stab_syms) {
      REQUIRE(group_action.apply(sym, bits) == bits);
    }
  }
}

template <typename bit_t> void test_representative(int n_sites) {
  auto group_action = GroupActionLookup<bit_t>(cyclic_group(n_sites));

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
  auto group_action = GroupActionLookup<bit_t>(cyclic_group(n_sites));
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
  auto group_action = GroupActionLookup<bit_t>(cyclic_group(n_sites));

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
  auto group_action = GroupActionLookup<bit_t>(cyclic_group(n_sites));
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

template <typename bit_t> void test_norm(int n_sites) {
  auto group_action = GroupActionLookup<bit_t>(cyclic_group(n_sites));
  std::vector<complex> chis(n_sites, 1.0);
  Representation irrep(chis);
  for (bit_t bits : Subsets(n_sites)) {
    double nrm = norm(bits, group_action, irrep);
    auto stabilizer = stabilizer_symmetries(bits, group_action);
    REQUIRE(close(nrm * nrm, (double)stabilizer.size()));
  }
}

template <typename bit_t>
void test_representatives_indices_symmetries_limits(int n_sites) {

  auto group_action = GroupActionLookup<bit_t>(cyclic_group(n_sites));

  for (int npar = 0; npar <= n_sites; ++npar) {
    auto lintable = indexing::LinTable<bit_t>(n_sites, npar);
    auto [reps, idces, syms, limits] =
        representatives_indices_symmetries_limits<bit_t>(
            indexing::CombinationsIndexing<bit_t>(n_sites, npar), group_action);
    idx_t n_reps = reps.size();
    for (idx_t i = 0; i < n_reps; ++i) {
      REQUIRE(reps[i] == representative(reps[i], group_action));
    }

    idx_t idx = 0;
    for (bit_t state : Combinations<bit_t>(n_sites, npar)) {
      idx_t k = idces[idx];
      bit_t rep = reps[k];
      REQUIRE(rep == representative(state, group_action));

      auto [l, u] = limits[idx];
      for (span_size_t i = l; i < u; ++i) {
        REQUIRE(rep == group_action.apply(syms[i], state));
      }

      ++idx;
    }
  }
}



TEST_CASE("symmetry_operations", "[symmetries]") {
  using namespace hydra;

  Log("Testing symmetry_operations");
  int max_N = 6;

  for (int n_sites = 0; n_sites <= max_N; ++n_sites) {
    test_stabilizer_symmetries<uint16_t>(n_sites);
    test_stabilizer_symmetries<uint32_t>(n_sites);
    test_stabilizer_symmetries<uint64_t>(n_sites);
  }

  for (int n_sites = 0; n_sites <= max_N; ++n_sites) {
    test_representative<uint16_t>(n_sites);
    test_representative<uint32_t>(n_sites);
    test_representative<uint64_t>(n_sites);
  }

  for (int n_sites = 0; n_sites <= max_N; ++n_sites) {
    test_representative_subset<uint16_t>(n_sites);
    test_representative_subset<uint32_t>(n_sites);
    test_representative_subset<uint64_t>(n_sites);
  }

  for (int n_sites = 0; n_sites <= max_N; ++n_sites) {
    test_representative_sym<uint16_t>(n_sites);
    test_representative_sym<uint32_t>(n_sites);
    test_representative_sym<uint64_t>(n_sites);
  }

  for (int n_sites = 0; n_sites <= max_N; ++n_sites) {
    test_representative_sym_subset<uint16_t>(n_sites);
    test_representative_sym_subset<uint32_t>(n_sites);
    test_representative_sym_subset<uint64_t>(n_sites);
  }

  for (int n_sites = 0; n_sites <= max_N; ++n_sites) {
    test_norm<uint16_t>(n_sites);
    test_norm<uint32_t>(n_sites);
    test_norm<uint64_t>(n_sites);
  }

  for (int n_sites = 1; n_sites <= max_N; ++n_sites) {
    test_representatives_indices_symmetries_limits<uint16_t>(n_sites);
    test_representatives_indices_symmetries_limits<uint32_t>(n_sites);
    test_representatives_indices_symmetries_limits<uint64_t>(n_sites);
  }


}
