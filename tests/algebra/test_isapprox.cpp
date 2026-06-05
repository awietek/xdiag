// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <optional>

#include <tests/catch.hpp>

#include <xdiag/algebra/algebras/electron_algebra.hpp>
#include <xdiag/algebra/algebras/spin_algebra.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;
using namespace xdiag::algebra;

TEST_CASE("isapprox", "[operators]") try {
  using namespace arma;

  Algebra spin = spin_algebra();
  Algebra elec = electron_algebra();

  Log("Testing isapprox of operators");

  // SdotS is symmetric under exchange of its two sites
  REQUIRE(
      isapprox(OpSum(Op("SdotS", {0, 1})), OpSum(Op("SdotS", {1, 0})), spin));

  // ScalarChirality is invariant under cyclic permutations of its sites
  REQUIRE(isapprox(OpSum(Op("ScalarChirality", {0, 1, 2})),
                   OpSum(Op("ScalarChirality", {1, 2, 0})), spin));

  // ... and flips sign under an odd permutation
  REQUIRE(!isapprox(OpSum(Op("ScalarChirality", {0, 1, 2})),
                    OpSum(Op("ScalarChirality", {2, 1, 0})), spin));

  // ScalarChirality{2,1,0} = -ScalarChirality{0,1,2}
  OpSum os1 = 1.0 * Op("ScalarChirality", {0, 1, 2});
  OpSum os2 = -1.0 * Op("ScalarChirality", {2, 1, 0});
  REQUIRE(isapprox(os1, os2, spin));

  os1 = 1.0 * Op("ScalarChirality", {0, 1, 2});
  os2 = 1.0 * Op("ScalarChirality", {2, 1, 0});
  std::optional<Scalar> factor = isapprox_multiple(os1, os2, spin);
  REQUIRE(factor);
  REQUIRE(isapprox(*factor, -1.0));

  // The bond operators Exchange (spin) and Hop / Hopup / Hopdn (electron) are
  // hermitian and symmetric under exchange of their two sites: O{0,1} == O{1,0}.
  std::vector<std::string> types = {"Exchange", "Hop", "Hopup", "Hopdn"};
  for (std::string const &type : types) {
    Algebra alg = (type == "Exchange") ? spin : elec;

    // symmetric under site exchange (same complex coefficient)
    os1 = complex(0, 1) * Op(type, {0, 1});
    os2 = complex(0, 1) * Op(type, {1, 0});
    REQUIRE(isapprox(os1, os2, alg));

    // a relative sign between the two orderings makes them differ
    os2 = complex(0, -1) * Op(type, {1, 0});
    REQUIRE(!isapprox(os1, os2, alg));
  }

  Log("done");
} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
