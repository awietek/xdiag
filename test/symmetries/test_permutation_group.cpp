#include "../catch.hpp"

#include <iostream>

#include <hydra/all.h>

TEST_CASE("permutation_group", "[symmetries]") {
  using namespace hydra;
  Log("Test PermutationGroup");

  std::string lfile = "data/triangular.j1j2jch/"
                      "triangular.12.j1j2jch.sublattices.fsl.lat";
  auto group = PermutationGroup(read_permutations(lfile));

  for (int sym = 0; sym < group.size(); ++sym) {
    auto p = group[sym];
    auto pinv = group[group.inverse(sym)];
    // HydraPrint(p * pinv);
    auto id = identity_permutation(group.n_sites());
    REQUIRE(p * pinv == id);
  }

  Log("done");
}
