#include "../../catch.hpp"

#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/combinations_indexing.hpp>
#include <xdiag/combinatorics/lin_table.hpp>
#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/combinatorics/subsets_indexing.hpp>
#include <xdiag/symmetries/group_action/group_action_lookup.hpp>
#include <xdiag/symmetries/operations/group_action_operations.hpp>
#include <xdiag/symmetries/operations/representative_list.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/utils/close.hpp>

using namespace xdiag;
using namespace xdiag::combinatorics;
using namespace xdiag::symmetries;

static PermutationGroup cyclic_group(int64_t nsites) {
  // test cyclic group
  std::vector<Permutation> permutation_array;
  for (int64_t sym = 0; sym < nsites; ++sym) {

    std::vector<int64_t> pv;
    for (int64_t site = 0; site < nsites; ++site) {
      int64_t newsite = (site + sym) % nsites;
      pv.push_back(newsite);
    }
    permutation_array.push_back(Permutation(pv));
  }
  return PermutationGroup(permutation_array);
}

template <typename bit_t> void test_stabilizer_symmetries(int64_t nsites) {
  auto group_action = GroupActionLookup<bit_t>(cyclic_group(nsites));

  for (bit_t bits : Subsets(nsites)) {
    auto stab_syms = stabilizer_symmetries(bits, group_action);
    for (int64_t sym : stab_syms) {
      REQUIRE(group_action.apply(sym, bits) == bits);
    }
  }
}

template <typename bit_t> void test_representative(int64_t nsites) {
  auto group_action = GroupActionLookup<bit_t>(cyclic_group(nsites));

  for (bit_t bits : Subsets(nsites)) {
    bit_t rep = representative(bits, group_action);
    int64_t hit_ctr = 0;
    int64_t n_sym = group_action.n_symmetries();
    for (int64_t sym = 0; sym < n_sym; ++sym) {
      REQUIRE(group_action.apply(sym, bits) >= rep);
      if (group_action.apply(sym, bits) == rep)
        hit_ctr++;
    }
    if (n_sym)
      REQUIRE(hit_ctr > 0);
  }
}

template <typename bit_t> void test_representative_subset(int64_t nsites) {
  auto group_action = GroupActionLookup<bit_t>(cyclic_group(nsites));
  std::vector<int64_t> subset;
  for (int64_t i = 0; i < nsites; i += 2)
    subset.push_back(i);

  for (bit_t bits : Subsets(nsites)) {
    bit_t rep = representative_subset(bits, group_action, subset);
    int64_t hit_ctr = 0;
    int64_t n_sym = subset.size();
    for (int64_t sym : subset) {
      REQUIRE(group_action.apply(sym, bits) >= rep);
      if (group_action.apply(sym, bits) == rep)
        hit_ctr++;
    }
    if (n_sym > 0)
      REQUIRE(hit_ctr > 0);
  }
}

template <typename bit_t> void test_representative_sym(int64_t nsites) {
  auto group_action = GroupActionLookup<bit_t>(cyclic_group(nsites));

  for (bit_t bits : Subsets(nsites)) {
    auto [rep, rsym] = representative_sym(bits, group_action);
    int64_t hit_ctr = 0;
    int64_t n_sym = group_action.n_symmetries();
    for (int64_t sym = 0; sym < n_sym; ++sym) {
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

template <typename bit_t> void test_representative_sym_subset(int64_t nsites) {
  auto group_action = GroupActionLookup<bit_t>(cyclic_group(nsites));
  std::vector<int64_t> subset;
  for (int64_t i = 0; i < nsites; i += 2)
    subset.push_back(i);

  for (bit_t bits : Subsets(nsites)) {
    auto [rep, rsym] = representative_sym_subset(bits, group_action, subset);
    int64_t hit_ctr = 0;
    int64_t n_sym = group_action.n_symmetries();
    for (int64_t sym : subset) {
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

template <typename bit_t> void test_norm(int64_t nsites) {
  auto group_action = GroupActionLookup<bit_t>(cyclic_group(nsites));
  arma::cx_vec chis(nsites, arma::fill::ones);
  Representation irrep(group_action.permutation_group(), chis);
  for (bit_t bits : Subsets(nsites)) {
    double nrm = norm(bits, group_action, chis);
    auto stabilizer = stabilizer_symmetries(bits, group_action);
    REQUIRE(close(nrm * nrm, (double)stabilizer.size()));
  }
}

template <typename bit_t>
void test_representatives_indices_symmetries_limits(int64_t nsites) {

  auto group_action = GroupActionLookup<bit_t>(cyclic_group(nsites));

  for (int64_t npar = 0; npar <= nsites; ++npar) {
    auto lintable = combinatorics::LinTable<bit_t>(nsites, npar);
    auto [reps, idces, syms, limits] =
        representatives_indices_symmetries_limits<bit_t>(
            combinatorics::CombinationsIndexing<bit_t>(nsites, npar),
            group_action);
    int64_t n_reps = reps.size();
    for (int64_t i = 0; i < n_reps; ++i) {
      REQUIRE(reps[i] == representative(reps[i], group_action));
    }

    int64_t idx = 0;
    for (bit_t state : Combinations<bit_t>(nsites, npar)) {
      int64_t k = idces[idx];
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
  using namespace xdiag;

  Log("Testing symmetry_operations");
  int64_t max_N = 6;

  for (int64_t nsites = 1; nsites <= max_N; ++nsites) {
    test_stabilizer_symmetries<uint16_t>(nsites);
    test_stabilizer_symmetries<uint32_t>(nsites);
    test_stabilizer_symmetries<uint64_t>(nsites);
  }

  for (int64_t nsites = 1; nsites <= max_N; ++nsites) {
    test_representative<uint16_t>(nsites);
    test_representative<uint32_t>(nsites);
    test_representative<uint64_t>(nsites);
  }

  for (int64_t nsites = 1; nsites <= max_N; ++nsites) {
    test_representative_subset<uint16_t>(nsites);
    test_representative_subset<uint32_t>(nsites);
    test_representative_subset<uint64_t>(nsites);
  }

  for (int64_t nsites = 1; nsites <= max_N; ++nsites) {
    test_representative_sym<uint16_t>(nsites);
    test_representative_sym<uint32_t>(nsites);
    test_representative_sym<uint64_t>(nsites);
  }

  for (int64_t nsites = 1; nsites <= max_N; ++nsites) {
    test_representative_sym_subset<uint16_t>(nsites);
    test_representative_sym_subset<uint32_t>(nsites);
    test_representative_sym_subset<uint64_t>(nsites);
  }

  for (int64_t nsites = 1; nsites <= max_N; ++nsites) {
    test_norm<uint16_t>(nsites);
    test_norm<uint32_t>(nsites);
    test_norm<uint64_t>(nsites);
  }

  for (int64_t nsites = 1; nsites <= max_N; ++nsites) {
    test_representatives_indices_symmetries_limits<uint16_t>(nsites);
    test_representatives_indices_symmetries_limits<uint32_t>(nsites);
    test_representatives_indices_symmetries_limits<uint64_t>(nsites);
  }
}
