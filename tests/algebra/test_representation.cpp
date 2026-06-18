// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <cmath>
#include <optional>

#include "../catch.hpp"

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/representation.hpp>
#include <xdiag/algebra/symmetrize.hpp>
#include <xdiag/matrices/matrix.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/symmetries/representation_set.hpp>
#include <xdiag/utils/error.hpp>

using namespace xdiag;
using namespace xdiag::algebra;
using namespace arma;

TEST_CASE("opsumrepresentation", "[algebra]") try {

  // ===========================================================================
  // representation() — U(1) / charge sectors
  // ===========================================================================
  // For a charge Representation the type names the action ("nup"/"ndn"/"np").
  // The charge passed in is ignored and recomputed from the OpSum. The algebra
  // only matters for "Matrix" Ops (via its local dimension d); d == 2 suffices.
  auto nup = [](OpSum const &ops) -> std::optional<int64_t> {
    std::optional<Representation> rep =
        representation(ops, Representation("nup", 0), matrix_algebra(2, 2));
    return rep ? std::optional<int64_t>(rep->charge()) : std::nullopt;
  };
  auto ndn = [](OpSum const &ops) -> std::optional<int64_t> {
    std::optional<Representation> rep =
        representation(ops, Representation("ndn", 0), matrix_algebra(2, 2));
    return rep ? std::optional<int64_t>(rep->charge()) : std::nullopt;
  };

  // The result is itself a charge Representation labelled by the action.
  {
    std::optional<Representation> rep = representation(
        OpSum(Op("S+", 0)), Representation("nup", 999), matrix_algebra(2, 2));
    REQUIRE(rep);
    REQUIRE(rep->is_charge());
    REQUIRE(rep->type() == "nup");
    REQUIRE(rep->charge() == 1); // input charge 999 was ignored
  }

  REQUIRE(*nup(OpSum(Op("S+", 0))) == 1);
  REQUIRE(*ndn(OpSum(Op("S+", 0))) == -1);
  REQUIRE(*nup(OpSum(Op("S-", 0))) == -1);
  REQUIRE(*ndn(OpSum(Op("S-", 0))) == 1);

  REQUIRE(*nup(OpSum(Op("Cdagup", 0))) == 1);
  REQUIRE(*ndn(OpSum(Op("Cdagup", 0))) == 0);
  REQUIRE(*nup(OpSum(Op("Cup", 0))) == -1);
  REQUIRE(*ndn(OpSum(Op("Cdagdn", 0))) == 1);
  REQUIRE(*ndn(OpSum(Op("Cdn", 0))) == -1);

  REQUIRE(*nup(OpSum(Op("SdotS", {0, 1}))) == 0);
  REQUIRE(*nup(OpSum(Op("Sz", 0))) == 0);

  // Different monomials carrying different charges -> no well-defined sector
  REQUIRE(nup(OpSum(Op("S+", 0)) + OpSum(Op("Sz", 0))) == std::nullopt);

  // "Matrix" Ops: charge from the particle-number difference of nonzero entries
  mat sz({{0.5, 0.0}, {0.0, -0.5}});
  mat sm({{0., 1.}, {0., 0.}});
  mat sp({{0., 0.}, {1., 0.}});
  cx_mat sy(mat({{0., 0.}, {0., 0.}}), mat({{0., -0.5}, {0.5, 0.}}));

  REQUIRE(*nup(OpSum(Op("Matrix", 0, sz))) == 0);
  REQUIRE(*nup(OpSum(Op("Matrix", 0, sp))) == 1);
  REQUIRE(*nup(OpSum(Op("Matrix", 0, sm))) == -1);
  REQUIRE(nup(OpSum(Op("Matrix", 0, sy))) == std::nullopt);
  REQUIRE(*nup(OpSum(Op("Matrix", {0, 1}, mat(kron(sp, sp))))) == 2);
  REQUIRE(*nup(OpSum(Op("Matrix", {0, 1}, mat(kron(sp, sm))))) == 0);

  // ===========================================================================
  // representation() — PermutationGroup irreps
  // ===========================================================================
  // Symmetrizing an Op with an irrep yields an OpSum that transforms under that
  // irrep; representation() must recover it from the group action alone.
  for (int64_t nsites = 3; nsites < 7; ++nsites) {
    Representation trivial = cyclic_group_irrep(nsites, 0);
    for (int64_t k = 0; k < nsites; ++k) {
      Representation irrep = cyclic_group_irrep(nsites, k);
      OpSum ops = symmetrize(Op("Sz", 0), irrep);

      std::optional<Representation> rep =
          representation(ops, irrep, spin_algebra(nsites));
      REQUIRE(rep);
      REQUIRE(rep->is_permutation());
      REQUIRE(isapprox(*rep, irrep));

      // Only the GROUP of the input is used; its characters are ignored.
      // Passing the trivial irrep (same group, different characters) recovers
      // the same result.
      std::optional<Representation> rep2 =
          representation(ops, trivial, spin_algebra(nsites));
      REQUIRE(rep2);
      REQUIRE(isapprox(*rep2, irrep));
    }
  }

  // An OpSum that is not covariant under the group has no well-defined sector:
  // a single-site Sz is not mapped to a multiple of itself by translations.
  {
    std::optional<Representation> none = representation(
        OpSum(Op("Sz", 0)), cyclic_group_irrep(4, 0), spin_algebra(4));
    REQUIRE_FALSE(none);
  }

  // ===========================================================================
  // representations() — a whole RepresentationSet at once
  // ===========================================================================
  {
    Representation perm = cyclic_group_irrep(4, 1);
    OpSum ops = symmetrize(Op("Sz", 0), perm);

    // Bogus input charges (ignored); the three types are
    // SitePermutation/nup/ndn
    RepresentationSet irreps(
        {perm, Representation("nup", 77), Representation("ndn", 99)});
    RepresentationSet result = representations(ops, irreps, spin_algebra(4));

    REQUIRE(result.size() == 3);
    REQUIRE(result.charge("nup") == 0); // Sz terms conserve nup and ndn
    REQUIRE(result.charge("ndn") == 0);
    REQUIRE(result.has_type("SitePermutation"));
    REQUIRE(isapprox(Representation(*result.group("SitePermutation"),
                                    *result.characters("SitePermutation")),
                     perm));
  }

  // Symmetries under which the OpSum has no well-defined sector are dropped.
  {
    OpSum mixed = OpSum(Op("S+", 0)) + OpSum(Op("Sz", 0)); // nup ill-defined
    RepresentationSet irreps({Representation("nup", 0)});
    RepresentationSet result =
        representations(mixed, irreps, matrix_algebra(1, 2));
    REQUIRE(result.size() == 0);
  }

} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
