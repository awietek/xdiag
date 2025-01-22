#include "../../catch.hpp"

using namespace xdiag;
using namespace xdiag::combinatorics;
using namespace xdiag::symmetries;

static PermutationGroup cyclic_group(int nsites) {
  // test cyclic group
  std::vector<Permutation> permutation_array;
  for (int sym = 0; sym < nsites; ++sym) {

    std::vector<int> pv;
    for (int site = 0; site < nsites; ++site) {
      int newsite = (site + sym) % nsites;
      pv.push_back(newsite);
    }
    permutation_array.push_back(Permutation(pv));
  }
  return PermutationGroup(permutation_array);
}

template <typename bit_t>
void test_representatives_indices_symmetries_limits_omp(int nsites) {

  auto group_action = GroupActionLookup<bit_t>(cyclic_group(nsites));

  for (int npar = 0; npar <= nsites; ++npar) {
    auto lintable = indexing::LinTable<bit_t>(nsites, npar);
    auto [reps1, idces1, syms1, limits1] =
        representatives_indices_symmetries_limits<bit_t>(
            indexing::CombinationsIndexing<bit_t>(nsites, npar), group_action);
    auto [reps2, idces2, syms2, limits2] =
        representatives_indices_symmetries_limits_omp<bit_t>(
            indexing::CombinationsIndexing<bit_t>(nsites, npar), group_action);

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

  for (int nsites = 0; nsites <= max_N; ++nsites) {
    test_representatives_indices_symmetries_limits_omp<uint16_t>(nsites);
    test_representatives_indices_symmetries_limits_omp<uint32_t>(nsites);
    test_representatives_indices_symmetries_limits_omp<uint64_t>(nsites);
  }
#endif
}
