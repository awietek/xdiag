#include "../../catch.hpp"

using namespace xdiag;
using namespace xdiag::combinatorics;
using namespace xdiag::symmetries;

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

template <typename bit_t>
void test_representatives_indices_symmetries_limits_omp(int n_sites) {

  auto group_action = GroupActionLookup<bit_t>(cyclic_group(n_sites));

  for (int npar = 0; npar <= n_sites; ++npar) {
    auto lintable = indexing::LinTable<bit_t>(n_sites, npar);
    auto [reps1, idces1, syms1, limits1] =
        representatives_indices_symmetries_limits<bit_t>(
            indexing::CombinationsIndexing<bit_t>(n_sites, npar), group_action);
    auto [reps2, idces2, syms2, limits2] =
        representatives_indices_symmetries_limits_omp<bit_t>(
            indexing::CombinationsIndexing<bit_t>(n_sites, npar), group_action);

    REQUIRE(reps1 == reps2);
    REQUIRE(idces1 == idces2);
    REQUIRE(syms1 == syms2);
    REQUIRE(limits1 == limits2);
  }
}

TEST_CASE("representative_list_omp", "[symmetries]") {
  using namespace xdiag;

#ifdef _OPENMP
  Log("testing representative_list_omp");
  int max_N = 6;

  for (int n_sites = 0; n_sites <= max_N; ++n_sites) {
    test_representatives_indices_symmetries_limits_omp<uint16_t>(n_sites);
    test_representatives_indices_symmetries_limits_omp<uint32_t>(n_sites);
    test_representatives_indices_symmetries_limits_omp<uint64_t>(n_sites);
  }
#endif
}
