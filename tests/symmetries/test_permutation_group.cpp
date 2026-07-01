// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <iostream>

#include <xdiag/config.hpp>
#include <xdiag/io/file_toml.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/utils/logger.hpp>

TEST_CASE("permutation_group", "[symmetries]") try {
  using namespace xdiag;
  Log("Test PermutationGroup");

  // Test inverse lookup on a file-loaded group
  {
    std::string lfile =
        XDIAG_DIRECTORY "/misc/data/triangular.j1j2jch/"
                        "triangular.12.j1j2jch.sublattices.fsl.toml";
    auto fl = FileToml(lfile);
    auto group = fl["Symmetries"].as<PermutationGroup>();
    for (int sym = 0; sym < group.size(); ++sym) {
      auto p = group[sym];
      auto pinv = group[group.inv(sym)];
      auto id = Permutation(group.nsites());
      REQUIRE(p * pinv == id);
    }
  }

  for (int64_t n = 2; n < 7; ++n) {
    auto group = cyclic_group(n);

    // Construction from arma::Mat (columns = permutations) matches vector ctor
    {
      arma::Mat<int64_t> mat(n, n);
      for (int64_t sym = 0; sym < n; ++sym)
        for (int64_t i = 0; i < n; ++i)
          mat(i, sym) = group[sym][i];
      REQUIRE(PermutationGroup(mat) == group);
    }

    // ptr(sym) points to the correct permutation data
    for (int64_t sym = 0; sym < group.size(); ++sym) {
      auto p = group[sym];
      auto const *data = group.ptr(sym);
      for (int64_t i = 0; i < group.nsites(); ++i)
        REQUIRE(data[i] == p[i]);
    }

    // multiply(s1, s2) is consistent with permutation product
    for (int64_t s1 = 0; s1 < group.size(); ++s1)
      for (int64_t s2 = 0; s2 < group.size(); ++s2)
        REQUIRE(group[group.multiply(s1, s2)] == group[s1] * group[s2]);

    // operator== / !=
    REQUIRE(group == cyclic_group(n));
    REQUIRE(group != cyclic_group(n + 1));
  }

  // subgroup: even-index elements of Z_6 form a valid subgroup
  {
    auto group = cyclic_group(6);
    auto sub = subgroup(group, {0, 2, 4});
    REQUIRE(sub.size() == 3);
    REQUIRE(sub.nsites() == 6);
    for (int64_t sym = 0; sym < sub.size(); ++sym)
      REQUIRE(sub[sym] * sub[sub.inv(sym)] == Permutation(sub.nsites()));
  }

  // Invalid group construction should throw
  // REQUIRE_THROWS(PermutationGroup(std::vector<Permutation>{}));
  // REQUIRE_THROWS(PermutationGroup({Permutation({1, 0})})); // missing identity

  Log("done");
} catch (xdiag::Error e) {
  error_trace(e);
}
