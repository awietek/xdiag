#include <tests/catch.hpp>

#include <vector>

#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/symmetries/representation_set.hpp>
#include <xdiag/utils/error.hpp>

using namespace xdiag;

TEST_CASE("representation_set", "[symmetries]") try {

  // ---------------------------------------------------------------------------
  // Construction and per-type lookup
  // ---------------------------------------------------------------------------
  {
    RepresentationSet set(
        {Representation("Nup", 2), Representation("Ndn", 3)});
    REQUIRE(set.size() == 2);
    REQUIRE(set.has_type("Nup"));
    REQUIRE(set.has_type("Ndn"));
    REQUIRE_FALSE(set.has_type("Ntot"));
    REQUIRE(set.charge("Nup") == 2);
    REQUIRE(set.charge("Ndn") == 3);
    REQUIRE(set.isreal());
    REQUIRE(isreal(set));

    // unknown type -> std::nullopt
    REQUIRE(set.charge("Ntot") == std::nullopt);
    // charge representation has no group / characters -> std::nullopt
    REQUIRE(set.group("Nup") == std::nullopt);
    REQUIRE(set.characters("Nup") == std::nullopt);
  }

  // mixing a permutation representation with charge representations
  {
    PermutationGroup group = cyclic_group(4);
    RepresentationSet mixed(
        {Representation(group), Representation("Nup", 1)});
    REQUIRE(mixed.size() == 2);
    REQUIRE(mixed.group("SitePermutation") == group);
    REQUIRE(mixed.characters("SitePermutation").has_value());
    REQUIRE(mixed.characters("SitePermutation").value().size() == 4);
    REQUIRE(mixed.charge("Nup") == 1);
    // wrong-kind / unknown lookups yield std::nullopt
    REQUIRE(mixed.charge("SitePermutation") == std::nullopt);
    REQUIRE(mixed.group("Nup") == std::nullopt);
  }

  // default-constructed set is empty and real
  {
    RepresentationSet empty;
    REQUIRE(empty.size() == 0);
    REQUIRE(empty.isreal());
  }

  // duplicate types are rejected (the set is keyed by type)
  REQUIRE_THROWS(RepresentationSet(
      {Representation("Nup", 1), Representation("Nup", 2)}));

  // ---------------------------------------------------------------------------
  // Equality is independent of insertion order
  // ---------------------------------------------------------------------------
  {
    RepresentationSet a(
        {Representation("Nup", 2), Representation("Ndn", 3)});
    RepresentationSet b(
        {Representation("Ndn", 3), Representation("Nup", 2)});
    REQUIRE(a == b);
    REQUIRE_FALSE(a != b);

    // different charge -> unequal
    RepresentationSet c(
        {Representation("Nup", 2), Representation("Ndn", 4)});
    REQUIRE(a != c);

    // different size -> unequal
    RepresentationSet d({Representation("Nup", 2)});
    REQUIRE(a != d);

    // same size, different type -> unequal
    RepresentationSet e(
        {Representation("Nup", 2), Representation("Sz", 3)});
    REQUIRE(a != e);
  }

  // isreal is false as soon as one member is complex
  {
    RepresentationSet cplx({cyclic_group_irrep(5, 1)});
    REQUIRE_FALSE(cplx.isreal());
    REQUIRE_FALSE(isreal(cplx));
  }

  // ---------------------------------------------------------------------------
  // isapprox
  // ---------------------------------------------------------------------------
  {
    RepresentationSet a(
        {Representation("Nup", 2), Representation("Ndn", 3)});
    RepresentationSet b(
        {Representation("Ndn", 3), Representation("Nup", 2)});
    RepresentationSet c(
        {Representation("Nup", 2), Representation("Ndn", 4)});
    RepresentationSet d({Representation("Nup", 2)});
    RepresentationSet e(
        {Representation("Nup", 2), Representation("Sz", 3)});

    REQUIRE(isapprox(a, b));        // order-independent
    REQUIRE_FALSE(isapprox(a, c));  // different charge
    REQUIRE_FALSE(isapprox(a, d));  // different size
    REQUIRE_FALSE(isapprox(a, e));  // different type set

    // member-function form agrees with the free function
    REQUIRE(a.isapprox(b));
    REQUIRE_FALSE(a.isapprox(c));
  }

  // ---------------------------------------------------------------------------
  // multiply / operator*
  // ---------------------------------------------------------------------------
  {
    // charge representations: charges add per type
    RepresentationSet s1(
        {Representation("Nup", 2), Representation("Ndn", 1)});
    RepresentationSet s2(
        {Representation("Nup", 5), Representation("Ndn", 4)});

    RepresentationSet prod = s1 * s2;
    REQUIRE(prod.charge("Nup") == 7);
    REQUIRE(prod.charge("Ndn") == 5);
    // multiply() and operator* agree, also independent of order
    REQUIRE(multiply(s1, s2) == prod);
    REQUIRE(s2 * s1 == prod);
    REQUIRE(s1.multiply(s2) == prod);

    // permutation representations: characters multiply per type
    RepresentationSet p1({cyclic_group_irrep(6, 2)});
    RepresentationSet p2({cyclic_group_irrep(6, 3)});
    RepresentationSet pprod = p1 * p2;
    REQUIRE(pprod.characters("SitePermutation").has_value());
    REQUIRE(isapprox(pprod.characters("SitePermutation").value(),
                     cyclic_group_irrep(6, 5).characters()));

    // multiply throws on incompatible sets
    // size mismatch
    REQUIRE_THROWS(multiply(s1, RepresentationSet({Representation("Nup", 1)})));
    // same size, mismatched types
    REQUIRE_THROWS(multiply(
        s1, RepresentationSet(
                {Representation("Nup", 1), Representation("Sz", 1)})));
  }

} catch (xdiag::Error e) {
  xdiag::error_trace(e);
  throw;
}
