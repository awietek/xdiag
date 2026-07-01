#include <tests/catch.hpp>

#include <sstream>
#include <vector>

#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/representation.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;
using namespace arma;

TEST_CASE("representation", "[symmetries]") try {

  // ---------------------------------------------------------------------------
  // Permutation-type representations (group + characters)
  // ---------------------------------------------------------------------------

  // Trivial representation: all-ones real characters
  for (int64_t n = 2; n < 7; ++n) {
    PermutationGroup group = cyclic_group(n);
    Representation trivial = Representation(group);
    REQUIRE(trivial.type() == "SitePermutation");
    REQUIRE(trivial.is_permutation());
    REQUIRE_FALSE(trivial.is_charge());
    REQUIRE(trivial.group() == group);
    REQUIRE(trivial.isreal());
    REQUIRE(trivial.characters().size() == n);
    for (double c : trivial.characters().as<arma::vec>()) {
      REQUIRE(c == Approx(1.0));
    }
    // group()/characters() are defined, charge() is not
    REQUIRE_THROWS(trivial.charge());
  }

  // isreal distinguishes real (k=0, k=n/2) from complex (other k) irreps
  for (int64_t n = 5; n < 8; ++n) {
    REQUIRE(cyclic_group_irrep(n, 0).isreal());
    REQUIRE(isreal(cyclic_group_irrep(n, 0)));
    REQUIRE_FALSE(cyclic_group_irrep(n, 1).isreal());
    REQUIRE_FALSE(isreal(cyclic_group_irrep(n, 1)));
  }
  // pi-momentum irrep of an even group is real
  REQUIRE(cyclic_group_irrep(6, 3).isreal());

  // ---------------------------------------------------------------------------
  // Charge-type representations (type label + integer charge)
  // ---------------------------------------------------------------------------
  {
    Representation r = Representation("Nup", 3);
    REQUIRE(r.type() == "Nup");
    REQUIRE(r.is_charge());
    REQUIRE_FALSE(r.is_permutation());
    REQUIRE(r.charge() == 3);
    REQUIRE(r.isreal()); // charge reps are always real
    // charge() is defined, group()/characters() are not
    REQUIRE_THROWS(r.group());
    REQUIRE_THROWS(r.characters());
  }

  // Default-constructed Representation has an (empty) type
  {
    Representation r;
    REQUIRE(r.isreal());
    REQUIRE(r.type().empty());
  }

  // ---------------------------------------------------------------------------
  // operator== / operator!=
  // ---------------------------------------------------------------------------
  {
    REQUIRE(cyclic_group_irrep(6, 2) == cyclic_group_irrep(6, 2));
    REQUIRE(cyclic_group_irrep(6, 2) != cyclic_group_irrep(6, 3));

    REQUIRE(Representation("Nup", 2) == Representation("Nup", 2));
    REQUIRE(Representation("Nup", 2) != Representation("Nup", 3));
    REQUIRE(Representation("Nup", 2) != Representation("Ndn", 2));

    // different kinds (charge vs permutation) are never equal
    REQUIRE(Representation("Nup", 0) != cyclic_group_irrep(4, 0));
  }

  // ---------------------------------------------------------------------------
  // isapprox
  // ---------------------------------------------------------------------------
  for (int64_t n = 3; n < 7; ++n) {
    Representation r0 = cyclic_group_irrep(n, 0);
    REQUIRE(isapprox(r0, cyclic_group_irrep(n, 0)));
    REQUIRE_FALSE(isapprox(r0, cyclic_group_irrep(n, 1)));
  }
  // charge reps
  REQUIRE(isapprox(Representation("Nup", 5), Representation("Nup", 5)));
  REQUIRE_FALSE(isapprox(Representation("Nup", 5), Representation("Nup", 6)));
  // different type strings are not approximately equal
  REQUIRE_FALSE(isapprox(Representation("Nup", 5), Representation("Ndn", 5)));
  REQUIRE_FALSE(isapprox(Representation("Nup", 0), cyclic_group_irrep(4, 0)));

  // ---------------------------------------------------------------------------
  // multiply / operator*
  // ---------------------------------------------------------------------------

  // irrep * conj(irrep) is always real
  for (int64_t n = 3; n < 8; ++n) {
    for (int64_t k = 0; k < n; ++k) {
      Representation irrep = cyclic_group_irrep(n, k);
      arma::cx_vec chars_hc = arma::conj(irrep.characters().as<arma::cx_vec>());
      Representation irrep_hc = Representation(irrep.group(), Vector(chars_hc));
      REQUIRE((irrep * irrep_hc).isreal());
    }
  }

  // C_n characters multiply: irrep(k1) * irrep(k2) == irrep((k1+k2) mod n)
  for (int64_t n = 3; n < 7; ++n) {
    for (int64_t k1 = 0; k1 < n; ++k1) {
      for (int64_t k2 = 0; k2 < n; ++k2) {
        Representation prod =
            multiply(cyclic_group_irrep(n, k1), cyclic_group_irrep(n, k2));
        REQUIRE(prod.is_permutation());
        REQUIRE(isapprox(prod, cyclic_group_irrep(n, (k1 + k2) % n)));
        // operator* agrees with multiply
        REQUIRE(isapprox(cyclic_group_irrep(n, k1) * cyclic_group_irrep(n, k2),
                         prod));
      }
    }
  }

  // charge reps: tensor product adds charges
  {
    Representation prod = Representation("Nup", 2) * Representation("Nup", 5);
    REQUIRE(prod.type() == "Nup");
    REQUIRE(prod.is_charge());
    REQUIRE(prod.charge() == 7);
  }

  // multiply throws when the two representations are incompatible
  {
    // different type strings
    REQUIRE_THROWS(
        multiply(Representation("Nup", 1), Representation("Ndn", 1)));
    // different kinds
    REQUIRE_THROWS(
        multiply(cyclic_group_irrep(4, 0), Representation("Nup", 1)));
    // same kind, different PermutationGroup
    REQUIRE_THROWS(
        multiply(cyclic_group_irrep(4, 0), cyclic_group_irrep(6, 0)));
  }

  // ---------------------------------------------------------------------------
  // Construction error handling
  // ---------------------------------------------------------------------------
  {
    PermutationGroup group = cyclic_group(4);
    // number of characters must match the group size
    REQUIRE_THROWS(Representation(group, std::vector<double>{1.0, 1.0}));
    // characters must obey c(g) * c(h) = c(gh)
    REQUIRE_THROWS(
        Representation(group, std::vector<double>{1.0, 2.0, 1.0, 2.0}));
  }

  // ---------------------------------------------------------------------------
  // Streaming / to_string
  // ---------------------------------------------------------------------------
  {
    std::ostringstream ss_real, ss_cplx, ss_charge;
    ss_real << cyclic_group_irrep(4, 0);
    ss_cplx << cyclic_group_irrep(5, 1);
    ss_charge << Representation("Nup", 2);
    REQUIRE_FALSE(ss_real.str().empty());
    REQUIRE_FALSE(ss_cplx.str().empty());
    REQUIRE_FALSE(ss_charge.str().empty());

    REQUIRE_FALSE(to_string(cyclic_group_irrep(4, 0)).empty());
    REQUIRE_FALSE(to_string(Representation("Nup", 2)).empty());
  }

} catch (xdiag::Error e) {
  xdiag::error_trace(e);
}
