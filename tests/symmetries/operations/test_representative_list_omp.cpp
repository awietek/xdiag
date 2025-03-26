#ifdef _OPENMP
#include "../../catch.hpp"

#include <xdiag/combinatorics/combinations.hpp>
#include <xdiag/combinatorics/combinations_indexing.hpp>
#include <xdiag/combinatorics/lin_table.hpp>
#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/combinatorics/subsets_indexing.hpp>
#include <xdiag/io/read.hpp>
#include <xdiag/symmetries/group_action/group_action_lookup.hpp>
#include <xdiag/symmetries/operations/representative_list.hpp>
#include <xdiag/symmetries/operations/representative_list_omp.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;
using namespace xdiag::combinatorics;
using namespace xdiag::symmetries;

template <typename bit_t>
void test_representatives_indices_symmetries_limits_omp() {
  std::string lfile =
      XDIAG_DIRECTORY "/misc/data/triangular.j1j2jch/"
                      "triangular.12.j1j2jch.sublattices.fsl.toml";

  int nsites = 12;
  auto fl = FileToml(lfile);
  // auto ops = fl["Interactions"].as<OpSum>();

  // auto irrep = read_representation(fl, "Gamma.C6.A");
  // auto characters = irrep.characters().as<arma::vec>();

  auto irrep = read_representation(fl, "K.C3.Eb");
  auto characters = irrep.characters().as<arma::cx_vec>();

  auto group_action = GroupActionLookup<bit_t>(irrep.group());

  for (int nup = 0; nup <= nsites; ++nup) {
    // auto rawindex = CombinationsIndexing<bit_t>(nsites, nup);

    auto rawindex = SubsetsIndexing<bit_t>(nsites);
    auto [reps1, idces1, syms1, limits1, norms1] =
        representatives_indices_symmetries_limits_norms<bit_t>(
            rawindex, group_action, characters);
    auto [reps2, idces2, syms2, limits2, norms2] =
        representatives_indices_symmetries_limits_norms_omp<bit_t>(
            rawindex, group_action, characters);

    REQUIRE(reps1 == reps2);
    REQUIRE(idces1 == idces2);
    REQUIRE(syms1 == syms2);
    REQUIRE(limits1 == limits2);
    REQUIRE(norms1 == norms2);
  }
}

TEST_CASE("representative_list_omp", "[symmetries]") {
  using namespace xdiag;

  Log("testing representative_list_omp");
  test_representatives_indices_symmetries_limits_omp<uint32_t>();
  test_representatives_indices_symmetries_limits_omp<uint64_t>();
}
#endif
