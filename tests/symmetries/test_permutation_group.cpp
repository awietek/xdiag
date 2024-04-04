#include "../catch.hpp"

#include <iostream>

#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>

TEST_CASE("permutation_group", "[symmetries]") {
  using namespace xdiag;
  Log("Test PermutationGroup");

  std::string lfile =
      XDIAG_DIRECTORY "/misc/data/triangular.j1j2jch/"
                      "triangular.12.j1j2jch.sublattices.fsl.lat";
  auto group = PermutationGroup(read_permutations(lfile));

  for (int sym = 0; sym < group.size(); ++sym) {
    auto p = group[sym];
    auto pinv = group[group.inverse(sym)];
    // XDiagPrint(p * pinv);
    auto id = identity_permutation(group.n_sites());
    REQUIRE(p * pinv == id);
  }

  Log("done");
}
