// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <vector>
#include <xdiag/operators/hc.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/error.hpp>

TEST_CASE("Monomial", "[operators]") try {
  using namespace xdiag;

  Op A("A", 0);
  Op B("B", 1);
  Op C("C", 2);

  // --- Default construction ---
  {
    Monomial m;
    REQUIRE(m.empty());
    REQUIRE(m.size() == 0);
  }

  // --- Construction from single Op ---
  {
    Monomial m(A);
    REQUIRE(m.size() == 1);
    REQUIRE(m[0] == A);
  }

  // --- Initializer-list construction ---
  {
    Monomial m{A, B, C};
    REQUIRE(m.size() == 3);
    REQUIRE(m[0] == A);
    REQUIRE(m[1] == B);
    REQUIRE(m[2] == C);
  }

  // --- Explicit vector construction ---
  {
    std::vector<Op> v = {A, B, C};
    Monomial m(v);
    REQUIRE(m.size() == 3);
    REQUIRE(m[0] == A);
    REQUIRE(m[1] == B);
    REQUIRE(m[2] == C);
  }

  // --- ops() accessor ---
  {
    Monomial m{A, B};
    std::vector<Op> ops = m.ops();
    REQUIRE(ops.size() == 2);
    REQUIRE(ops[0] == A);
    REQUIRE(ops[1] == B);
  }

  // --- begin()/end() range-for ---
  {
    Monomial m{A, B, C};
    std::vector<Op> collected;
    for (auto const &op : m) {
      collected.push_back(op);
    }
    REQUIRE(collected.size() == 3);
    REQUIRE(collected[0] == A);
    REQUIRE(collected[1] == B);
    REQUIRE(collected[2] == C);
  }

  // --- operator*=(Op): appends one op ---
  {
    Monomial m{A, B};
    m *= C;
    REQUIRE(m.size() == 3);
    REQUIRE(m[2] == C);
  }

  // --- operator*=(Monomial): appends all ops ---
  {
    Monomial m1{A, B};
    Monomial m2{C};
    m1 *= m2;
    REQUIRE(m1.size() == 3);
    REQUIRE(m1[2] == C);
  }

  // --- operator*(Op) and operator*(Monomial): non-mutating ---
  {
    Monomial m{A, B};
    Monomial r1 = m * C;
    REQUIRE(r1.size() == 3);
    REQUIRE(r1[2] == C);
    REQUIRE(m.size() == 2); // original unchanged

    Monomial m2{C};
    Monomial r2 = m * m2;
    REQUIRE(r2.size() == 3);
    REQUIRE(r2[0] == A);
    REQUIRE(r2[1] == B);
    REQUIRE(r2[2] == C);
  }

  // --- Free Op*Op -> Monomial of length 2 ---
  {
    Monomial m = A * B;
    REQUIRE(m.size() == 2);
    REQUIRE(m[0] == A);
    REQUIRE(m[1] == B);
  }

  // --- Free Op*Monomial -> length n+1 ---
  {
    Monomial m{B, C};
    Monomial r = A * m;
    REQUIRE(r.size() == 3);
    REQUIRE(r[0] == A);
    REQUIRE(r[1] == B);
    REQUIRE(r[2] == C);
  }

  // --- operator== / operator!= ---
  {
    Monomial m1{A, B, C};
    Monomial m2{A, B, C};
    Monomial m3{A, C, B};
    REQUIRE(m1 == m2);
    REQUIRE(m1 != m3);
  }

  // --- operator< (lexicographic) ---
  {
    // "A" < "B" lexicographically by type string
    Monomial mA{Op("A", 0)};
    Monomial mB{Op("B", 0)};
    REQUIRE(mA < mB);
    REQUIRE_FALSE(mB < mA);

    // Same type, lower site first
    Monomial mA0{Op("X", 0)};
    Monomial mA1{Op("X", 1)};
    REQUIRE(mA0 < mA1);
    REQUIRE_FALSE(mA1 < mA0);

    // Shorter prefix is less when prefix matches
    Monomial mShort{Op("X", 0)};
    Monomial mLong{Op("X", 0), Op("X", 1)};
    REQUIRE(mShort < mLong);
    REQUIRE_FALSE(mLong < mShort);
  }

  // --- hc(Monomial): hermitian conjugate ---
  {
    Op sp = Op("S+", 0);
    Op sm = Op("S-", 0);

    // Single-op: hc of S+ is S-
    Monomial m_sp{sp};
    Monomial m_sm{sm};
    REQUIRE(hc(m_sp) == OpSum(m_sm));
    REQUIRE(hc(m_sm) == OpSum(m_sp));

    // Multi-op: hc reverses order and conjugates each
    // hc({S+(0), S-(1)}) == {S+(1), S-(0)}
    Monomial m2{Op("S+", 0), Op("S-", 1)};
    Monomial m2_hc{Op("S+", 1), Op("S-", 0)};
    REQUIRE(hc(m2) == OpSum(m2_hc));
  }

} catch (xdiag::Error e) {
  error_trace(e);
}

TEST_CASE("Term", "[operators]") try {
  using namespace xdiag;

  Op A("A", 0);
  Op B("B", 1);

  Monomial m1{A, B};
  Monomial m2{A, B};
  Monomial m3{B, A};

  Coeff c1(1.0);
  Coeff c2(1.0);
  Coeff c3(2.0);

  Term t1{c1, m1};
  Term t2{c2, m2};
  Term t3{c3, m1};
  Term t4{c1, m3};

  // Same coeff and monomial
  REQUIRE(t1 == t2);

  // Different coeff
  REQUIRE(t1 != t3);

  // Different monomial
  REQUIRE(t1 != t4);

} catch (xdiag::Error e) {
  error_trace(e);
}
