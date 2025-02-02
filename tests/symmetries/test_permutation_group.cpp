#include "../catch.hpp"

#include <iostream>

#include <xdiag/combinatorics/subsets.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/symmetries/operations/symmetry_operations.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/utils/logger.hpp>

TEST_CASE("permutation_group", "[symmetries]") try {
  using namespace xdiag;
  Log("Test PermutationGroup");

  std::string lfile =
      XDIAG_DIRECTORY "/misc/data/triangular.j1j2jch/"
                      "triangular.12.j1j2jch.sublattices.fsl.toml";
  auto fl = FileToml(lfile);
  auto group = fl["Symmetries"].as<PermutationGroup>();
  for (int sym = 0; sym < group.size(); ++sym) {
    auto p = group[sym];
    auto pinv = group[group.inverse(sym)];
    // XDIAG_SHOW(p * pinv);
    auto id = Permutation(group.nsites());
    REQUIRE(p * pinv == id);
  }
  Log("done");
} catch (xdiag::Error e) {
  error_trace(e);
}
