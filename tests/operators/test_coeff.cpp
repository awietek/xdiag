// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/math/scalar.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/utils/error.hpp>

TEST_CASE("Coeff", "[operators]") try {
  using namespace xdiag;

  // --- Construction and isscalar/isstring ---

  Coeff cd(1.5);
  REQUIRE(cd.isscalar());
  REQUIRE_FALSE(cd.isstring());

  Coeff cc(complex(1.0, 2.0));
  REQUIRE(cc.isscalar());
  REQUIRE_FALSE(cc.isstring());

  Coeff cs(std::string("J"));
  REQUIRE(cs.isstring());
  REQUIRE_FALSE(cs.isscalar());

  Coeff ccs(Scalar(3.0));
  REQUIRE(ccs.isscalar());
  REQUIRE_FALSE(ccs.isstring());

  Coeff cchar("K");
  REQUIRE(cchar.isstring());
  REQUIRE_FALSE(cchar.isscalar());

  // --- scalar() accessor ---

  REQUIRE(cd.scalar() == Scalar(1.5));
  REQUIRE(cc.scalar() == Scalar(complex(1.0, 2.0)));
  REQUIRE(ccs.scalar() == Scalar(3.0));
  // REQUIRE_THROWS(cs.scalar()); // string Coeff: must throw

  // --- string() accessor ---

  REQUIRE(cs.string() == "J");
  REQUIRE(cchar.string() == "K");
  // REQUIRE_THROWS(cd.string()); // scalar Coeff: must throw

  // --- free-function isscalar/isstring ---

  REQUIRE(isscalar(cd));
  REQUIRE(isstring(cs));
  REQUIRE(scalar(cd) == Scalar(1.5));
  REQUIRE(string(cs) == "J");

  // --- operator== / operator!= ---

  REQUIRE(Coeff(1.0) == Coeff(1.0));
  // double 1.0 and complex(1.0,0.0) are different variant alternatives inside
  // Scalar
  REQUIRE(Coeff(1.0) != Coeff(complex(1.0, 0.0)));
  REQUIRE(Coeff("J") == Coeff("J"));
  REQUIRE(Coeff("J") != Coeff("K"));
  REQUIRE(Coeff(1.0) != Coeff("J"));

  // --- operator*(Coeff, Coeff) ---

  Coeff r1(2.0);
  Coeff r2(3.0);
  Coeff cx1(complex(1.0, 1.0));
  Coeff cx2(complex(0.0, 1.0));

  // double * double
  Coeff prod_rr = r1 * r2;
  REQUIRE(prod_rr.isscalar());
  REQUIRE(prod_rr.scalar() == Scalar(6.0));

  // double * complex
  Coeff prod_rc = r1 * cx1;
  REQUIRE(prod_rc.isscalar());
  REQUIRE(prod_rc.scalar() == Scalar(complex(2.0, 2.0)));

  // complex * complex
  Coeff prod_cc = cx1 * cx2; // (1+i)*i = i - 1 = (-1, 1)
  REQUIRE(prod_cc.isscalar());
  REQUIRE(prod_cc.scalar() == Scalar(complex(-1.0, 1.0)));

  // --- operator*(Coeff, Scalar) and operator*(Scalar, Coeff) ---

  Coeff prod_cs = r1 * Scalar(4.0);
  REQUIRE(prod_cs.scalar() == Scalar(8.0));

  Coeff prod_sc = Scalar(5.0) * r2;
  REQUIRE(prod_sc.scalar() == Scalar(15.0));

  // --- throws when either side is a string Coeff ---

  // REQUIRE_THROWS(cs * r1);
  // REQUIRE_THROWS(r1 * cs);
  // REQUIRE_THROWS(cs * cs);

} catch (xdiag::Error e) {
  error_trace(e);
}
