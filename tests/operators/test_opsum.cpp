// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <iostream>
#include <sstream>
#include <vector>
#include <xdiag/operators/collect.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>

TEST_CASE("OpSum construction", "[operators]") try {
  using namespace xdiag;

  Op A("A", int64_t(0));
  Op B("B", int64_t(0));
  Monomial mAB{A, B};

  // Empty
  {
    OpSum ops;
    REQUIRE(ops.size() == 0);
  }

  // From Op
  {
    OpSum ops(A);
    REQUIRE(ops.size() == 1);
  }

  // From Monomial
  {
    OpSum ops(mAB);
    REQUIRE(ops.size() == 1);
    REQUIRE(ops.terms()[0].monomial == mAB);
  }

  // (double, Op)
  {
    OpSum ops(2.0, A);
    REQUIRE(ops.size() == 1);
    REQUIRE(ops.terms()[0].coeff == Coeff(2.0));
  }

  // (complex, Op)
  {
    OpSum ops(complex(1.0, 2.0), A);
    REQUIRE(ops.size() == 1);
    REQUIRE(ops.terms()[0].coeff == Coeff(complex(1.0, 2.0)));
  }

  // (string, Op)
  {
    OpSum ops(std::string("J"), A);
    REQUIRE(ops.size() == 1);
    REQUIRE(ops.terms()[0].coeff == Coeff("J"));
  }

  // (Coeff, Monomial)
  {
    OpSum ops(Coeff(3.0), mAB);
    REQUIRE(ops.size() == 1);
    REQUIRE(ops.terms()[0].coeff == Coeff(3.0));
    REQUIRE(ops.terms()[0].monomial == mAB);
  }

} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("OpSum vector space", "[operators]") try {
  using namespace xdiag;

  Op A("A", int64_t(0));
  Op B("B", int64_t(1));
  Op C("C", int64_t(2));

  // += Op
  {
    OpSum ops;
    ops += A;
    ops += B;
    REQUIRE(ops.size() == 2);
  }

  // += OpSum
  {
    OpSum ops1(1.0, A);
    OpSum ops2(1.0, B);
    ops1 += ops2;
    REQUIRE(ops1.size() == 2);
  }

  // + operator (non-mutating)
  {
    OpSum ops1(1.0, A);
    OpSum ops2(1.0, B);
    OpSum ops3 = ops1 + ops2;
    REQUIRE(ops3.size() == 2);
    REQUIRE(ops1.size() == 1); // unchanged
  }

  // unary -: all coefficients negated
  {
    OpSum ops(2.0, A);
    OpSum neg = -ops;
    REQUIRE(neg.size() == 1);
    REQUIRE(neg.terms()[0].coeff == Coeff(-2.0));
  }

  // - operator
  {
    OpSum ops1(3.0, A);
    OpSum ops2(1.0, A);
    OpSum diff = ops1 - ops2;
    REQUIRE(diff.size() == 2);
  }

  // *=(double)
  {
    OpSum ops(2.0, A);
    ops *= 3.0;
    REQUIRE(ops.terms()[0].coeff == Coeff(6.0));
  }

  // /=(double)
  {
    OpSum ops(6.0, A);
    ops /= 2.0;
    REQUIRE(ops.terms()[0].coeff == Coeff(3.0));
  }

  // *=(complex): widens coefficient to complex
  {
    OpSum ops(2.0, A);
    ops *= complex(0.0, 1.0); // 2 * i = 2i
    REQUIRE(ops.terms()[0].coeff == Coeff(complex(0.0, 2.0)));
  }

} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("OpSum algebra product", "[operators]") try {
  using namespace xdiag;

  Op A("A", int64_t(0));
  Op B("B", int64_t(0));
  Op C("C", int64_t(0));
  Op D("D", int64_t(0));

  // Op * Op -> 2-op monomial via free algebra product on OpSum
  {
    OpSum lhs(1.0, A);
    OpSum rhs(1.0, B);
    OpSum prod = lhs * rhs;
    REQUIRE(prod.size() == 1);
    Monomial expected{A, B};
    REQUIRE(prod.terms()[0].monomial == expected);
  }

  // Distributive product: (A + B) * (C + D) -> 4 terms: {A,C},{A,D},{B,C},{B,D}
  {
    OpSum lhs = 1.0 * A + 1.0 * B;
    OpSum rhs = 1.0 * C + 1.0 * D;
    OpSum prod = lhs * rhs;
    REQUIRE(prod.size() == 4);

    Monomial mAC{A, C};
    Monomial mAD{A, D};
    Monomial mBC{B, C};
    Monomial mBD{B, D};

    // Collect all monomials from product terms
    std::vector<Monomial> monos;
    for (auto const &t : prod) {
      monos.push_back(t.monomial);
    }
    REQUIRE(std::find(monos.begin(), monos.end(), mAC) != monos.end());
    REQUIRE(std::find(monos.begin(), monos.end(), mAD) != monos.end());
    REQUIRE(std::find(monos.begin(), monos.end(), mBC) != monos.end());
    REQUIRE(std::find(monos.begin(), monos.end(), mBD) != monos.end());
  }

} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("OpSum named params", "[operators]") try {
  using namespace xdiag;

  Op X("X", int64_t(0));

  // Build with named param, set value, plain() resolves
  {
    OpSum ops = std::string("J") * X;
    ops["J"] = Scalar(2.0);
    OpSum resolved = ops.plain();
    REQUIRE(resolved.size() == 1);
    REQUIRE(resolved.terms()[0].coeff == Coeff(2.0));
  }

  // Accessing missing param throws (const accessor uses .at())
  {
    OpSum const ops(X);
    // REQUIRE_THROWS(ops["missing"]);
  }

  // J1/J2 Heisenberg chain (original test)
  {
    OpSum ops;
    ops += "J1" * Op("SdotS", {0, 1});
    ops += "J1" * Op("SdotS", {1, 2});
    ops += "J1" * Op("SdotS", {2, 3});
    ops += "J1" * Op("SdotS", {3, 4});
    ops += "J1" * Op("SdotS", {4, 5});
    ops += "J1" * Op("SdotS", {5, 0});
    ops += "J2" * Op("SdotS", {0, 2});
    ops += "J2" * Op("SdotS", {1, 3});
    ops += "J2" * Op("SdotS", {2, 4});
    ops += "J2" * Op("SdotS", {3, 5});
    ops += "J2" * Op("SdotS", {4, 0});
    ops += "J2" * Op("SdotS", {5, 1});

    double a = 1.0;
    complex b(2.0, 3.0);

    ops["J1"] = a;
    ops["J2"] = b;

    OpSum ops2;
    ops2 += a * Op("SdotS", {0, 1});
    ops2 += a * Op("SdotS", {1, 2});
    ops2 += a * Op("SdotS", {2, 3});
    ops2 += a * Op("SdotS", {3, 4});
    ops2 += a * Op("SdotS", {4, 5});
    ops2 += a * Op("SdotS", {5, 0});
    ops2 += b * Op("SdotS", {0, 2});
    ops2 += b * Op("SdotS", {1, 3});
    ops2 += b * Op("SdotS", {2, 4});
    ops2 += b * Op("SdotS", {3, 5});
    ops2 += b * Op("SdotS", {4, 0});
    ops2 += b * Op("SdotS", {5, 1});

    REQUIRE(ops.plain() == ops2);
  }

} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("OpSum access", "[operators]") try {
  using namespace xdiag;

  Op A("A", int64_t(0));
  Op B("B", int64_t(1));

  OpSum ops = 1.0 * A + 2.0 * B;

  // terms()
  REQUIRE(ops.terms().size() == 2);

  // size()
  REQUIRE(ops.size() == 2);

  // params() starts empty when no named coefficients used
  REQUIRE(ops.params().empty());

  // begin()/end() range-for yields Terms
  {
    int count = 0;
    for (auto const &t : ops) {
      REQUIRE(t.monomial.size() == 1);
      ++count;
    }
    REQUIRE(count == 2);
  }

} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("OpSum collect", "[operators]") try {
  using namespace xdiag;

  Op A("A", int64_t(0));

  // A + A -> 2*A (same monomial combined)
  {
    OpSum ops = 1.0 * A + 1.0 * A;
    OpSum c = collect(ops);
    REQUIRE(c.size() == 1);
    REQUIRE(c.terms()[0].coeff == Coeff(2.0));
  }

  // A + (-A) -> empty (zero coefficient removed)
  {
    OpSum ops = 1.0 * A + (-1.0) * A;
    OpSum c = collect(ops);
    REQUIRE(c.size() == 0);
  }

} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("OpSum equality", "[operators]") try {
  using namespace xdiag;

  Op A("A", int64_t(0));
  Op B("B", int64_t(1));

  OpSum ops1 = 1.0 * A + 2.0 * B;
  OpSum ops2 = 1.0 * A + 2.0 * B;
  OpSum ops3 = 1.0 * A + 3.0 * B;

  REQUIRE(ops1 == ops2);
  REQUIRE(ops1 != ops3);

} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("OpSum operand overloads", "[operators]") try {
  using namespace xdiag;

  Op A("A", int64_t(0));
  Op B("B", int64_t(1));
  Monomial mAB{A, B};

  // --- operator+ / += with Op and Monomial operands ---
  {
    OpSum ops(A);
    ops += B;    // += Op
    ops += mAB;  // += Monomial
    REQUIRE(ops.size() == 3);
  }
  {
    OpSum base(A);
    OpSum s1 = base + B;   // + Op
    OpSum s2 = base + mAB; // + Monomial
    REQUIRE(s1.size() == 2);
    REQUIRE(s2.size() == 2);
  }

  // --- operator- / -= with Op and Monomial operands ---
  {
    OpSum ops = 2.0 * A + 3.0 * B;
    ops -= A;    // -= Op
    ops -= mAB;  // -= Monomial
    REQUIRE(collect(ops).size() == 3); // A(coeff1), B, -(A*B)
  }
  {
    OpSum base = 2.0 * A;
    OpSum d1 = base - B;   // - Op
    OpSum d2 = base - mAB; // - Monomial
    REQUIRE(d1.size() == 2);
    REQUIRE(d2.size() == 2);
  }

  // --- unary minus negates every coefficient ---
  {
    OpSum ops = 2.0 * A + 3.0 * B;
    OpSum neg = -ops;
    OpSum sum = collect(ops + neg);
    REQUIRE(sum.size() == 0); // ops + (-ops) == 0
  }

  // --- algebra product with Op / Monomial / OpSum operands ---
  {
    OpSum a(A);
    OpSum pOp = a * B;     // OpSum * Op
    OpSum pMono = a * mAB; // OpSum * Monomial
    OpSum pSum = a * OpSum(B); // OpSum * OpSum
    REQUIRE(pOp.size() == 1);
    REQUIRE(pMono.size() == 1);
    REQUIRE(pSum.size() == 1);

    OpSum acc(A);
    acc *= B;   // *= Op
    REQUIRE(acc.size() == 1);
    OpSum acc2(A);
    acc2 *= mAB; // *= Monomial
    REQUIRE(acc2.size() == 1);
    OpSum acc3(A);
    acc3 *= OpSum(B); // *= OpSum
    REQUIRE(acc3.size() == 1);
  }

} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("OpSum scalar scaling paths", "[operators]") try {
  using namespace xdiag;

  Op A("A", int64_t(0));

  // --- /= complex ---
  {
    OpSum ops = complex(4.0, 0.0) * A;
    ops /= complex(2.0, 0.0);
    REQUIRE(collect(ops).terms()[0].coeff == Coeff(complex(2.0, 0.0)));
  }

  // --- scaling an OpSum carrying a named parameter resolves it via params_ ---
  {
    OpSum ops = std::string("J") * A;
    ops["J"] = Scalar(3.0);
    ops *= 2.0; // hits the params_.find branch in operator*=(Scalar)
    OpSum resolved = ops.plain();
    REQUIRE(resolved.terms()[0].coeff == Coeff(6.0));
  }

  // --- merging two OpSums that both carry named parameters ---
  {
    Op B("B", int64_t(1));
    OpSum o1 = std::string("J1") * A;
    o1["J1"] = Scalar(1.5);
    OpSum o2 = std::string("J2") * B;
    o2["J2"] = Scalar(2.5);
    o1 += o2; // merge_params combines the parameter maps
    REQUIRE(o1.params().size() == 2);
    OpSum resolved = o1.plain();
    REQUIRE(resolved.size() == 2);
  }

} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("OpSum free operators", "[operators]") try {
  using namespace xdiag;

  Op A("A", int64_t(0));
  Op B("B", int64_t(1));
  Monomial mAB{A, B};
  Monomial mBA{B, A};
  complex z(1.0, -2.0);

  auto single = [](OpSum const &ops, Coeff const &c, Monomial const &m) {
    return (ops.size() == 1) && (ops.terms()[0].coeff == c) &&
           (ops.terms()[0].monomial == m);
  };

  // --- scalar * Op and Op * scalar ---
  {
    REQUIRE(single(2.0 * A, Coeff(2.0), Monomial(A)));
    REQUIRE(single(A * 2.0, Coeff(2.0), Monomial(A)));
    REQUIRE(single(z * A, Coeff(z), Monomial(A)));
    REQUIRE(single(A * z, Coeff(z), Monomial(A)));
    REQUIRE(single(std::string("J") * A, Coeff("J"), Monomial(A)));
    REQUIRE(single(A * std::string("J"), Coeff("J"), Monomial(A)));
    REQUIRE(single(Coeff(3.0) * A, Coeff(3.0), Monomial(A)));
    REQUIRE(single(A * Coeff(3.0), Coeff(3.0), Monomial(A)));
    REQUIRE(single(Scalar(4.0) * A, Coeff(4.0), Monomial(A)));
    REQUIRE(single(A * Scalar(4.0), Coeff(4.0), Monomial(A)));
  }

  // --- scalar * Monomial and Monomial * scalar ---
  {
    REQUIRE(single(2.0 * mAB, Coeff(2.0), mAB));
    REQUIRE(single(mAB * 2.0, Coeff(2.0), mAB));
    REQUIRE(single(z * mAB, Coeff(z), mAB));
    REQUIRE(single(mAB * z, Coeff(z), mAB));
    REQUIRE(single(std::string("J") * mAB, Coeff("J"), mAB));
    REQUIRE(single(mAB * std::string("J"), Coeff("J"), mAB));
    REQUIRE(single(Coeff(3.0) * mAB, Coeff(3.0), mAB));
    REQUIRE(single(mAB * Coeff(3.0), Coeff(3.0), mAB));
    REQUIRE(single(Scalar(4.0) * mAB, Coeff(4.0), mAB));
    REQUIRE(single(mAB * Scalar(4.0), Coeff(4.0), mAB));
  }

  // --- scaling an OpSum from either side, and division ---
  {
    OpSum ops = 2.0 * A + 3.0 * B;
    OpSum expected = 4.0 * A + 6.0 * B;
    REQUIRE(2.0 * ops == expected);
    REQUIRE(ops * 2.0 == expected);
    REQUIRE(Scalar(2.0) * ops == expected);
    REQUIRE(ops * Scalar(2.0) == expected);
    REQUIRE(ops / 0.5 == expected);
    REQUIRE(ops / Scalar(0.5) == expected);

    OpSum expectedz = (2.0 * z) * A + (3.0 * z) * B;
    REQUIRE(z * ops == expectedz);
    REQUIRE(ops * z == expectedz);
    REQUIRE(ops * Scalar(z) == expectedz);

    complex zinv = 1.0 / z;
    OpSum expected_div = (2.0 * zinv) * A + (3.0 * zinv) * B;
    REQUIRE(ops / z == expected_div);
    REQUIRE(ops / Scalar(z) == expected_div);
    REQUIRE(ops == 2.0 * A + 3.0 * B); // operand untouched
  }

  // --- additive operators with a Monomial on the left ---
  {
    OpSum s = mAB + mBA;
    REQUIRE(s.size() == 2);
    REQUIRE(s.terms()[0].monomial == mAB);
    REQUIRE(s.terms()[1].monomial == mBA);
    REQUIRE(s.terms()[1].coeff == Coeff(1.0));

    OpSum d = mAB - mBA;
    REQUIRE(d.size() == 2);
    REQUIRE(d.terms()[1].coeff == Coeff(-1.0));

    OpSum n = -mAB;
    REQUIRE(single(n, Coeff(-1.0), mAB));

    OpSum rhs = 2.0 * A;
    OpSum s2 = mAB + rhs;
    REQUIRE(s2.size() == 2);
    REQUIRE(s2.terms()[0].monomial == mAB);
    REQUIRE(s2.terms()[1].coeff == Coeff(2.0));

    OpSum d2 = mAB - rhs;
    REQUIRE(d2.size() == 2);
    REQUIRE(d2.terms()[1].coeff == Coeff(-2.0));
  }

  // --- algebra products with Op / Monomial on the left ---
  {
    OpSum rhs = 2.0 * B + 3.0 * A;
    OpSum p = A * rhs;
    REQUIRE(p.size() == 2);
    REQUIRE(p.terms()[0].monomial == mAB);
    REQUIRE(p.terms()[0].coeff == Coeff(2.0));
    REQUIRE(p.terms()[1].monomial == Monomial{A, A});
    REQUIRE(p.terms()[1].coeff == Coeff(3.0));

    OpSum q = mAB * rhs;
    REQUIRE(q.size() == 2);
    REQUIRE(q.terms()[0].monomial == Monomial{A, B, B});
    REQUIRE(q.terms()[1].monomial == Monomial{A, B, A});
    REQUIRE(q.terms()[1].coeff == Coeff(3.0));

    // non-commutativity: A * ops != ops * A in general
    REQUIRE(A * rhs != rhs * A);
  }
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("OpSum constructors and assignment", "[operators]") try {
  using namespace xdiag;

  Op A("A", int64_t(0));
  Op B("B", int64_t(1));
  Monomial mAB{A, B};
  complex z(0.5, 2.0);

  // --- remaining constructor overloads ---
  {
    REQUIRE(OpSum(Scalar(2.0), A).terms()[0].coeff == Coeff(2.0));
    REQUIRE(OpSum(Scalar(z), mAB).terms()[0].coeff == Coeff(z));
    REQUIRE(OpSum(Coeff("J"), A).terms()[0].coeff == Coeff("J"));
    REQUIRE(OpSum(std::string("J"), mAB).terms()[0].coeff == Coeff("J"));
    REQUIRE(OpSum(2.0, mAB).terms()[0].monomial == mAB);
    REQUIRE(OpSum(z, mAB).terms()[0].coeff == Coeff(z));
    REQUIRE(OpSum(z, A).terms()[0].coeff == Coeff(z));
    REQUIRE(!OpSum(z, A).isreal());
    REQUIRE(OpSum(2.0, A).isreal());
    REQUIRE(OpSum().empty());
    REQUIRE(empty(OpSum()));
    REQUIRE(!empty(OpSum(A)));
  }

  // --- assignment from Op / Monomial replaces terms and clears params ---
  {
    OpSum ops = std::string("J") * A + 2.0 * B;
    ops["J"] = 1.5;
    REQUIRE(ops.params().size() == 1);

    ops = B;
    REQUIRE(ops.size() == 1);
    REQUIRE(ops.terms()[0].monomial == Monomial(B));
    REQUIRE(ops.terms()[0].coeff == Coeff(1.0));
    REQUIRE(ops.params().empty());

    ops = mAB;
    REQUIRE(ops.size() == 1);
    REQUIRE(ops.terms()[0].monomial == mAB);
  }

  // --- const parameter access ---
  {
    OpSum ops = std::string("J") * A;
    ops["J"] = Scalar(2.5);
    OpSum const &cops = ops;
    REQUIRE(cops["J"] == Scalar(2.5));
    REQUIRE_THROWS_AS(cops["missing"], Error);
  }

  // --- isreal with string coefficients depends on resolved values ---
  {
    OpSum ops = std::string("J") * A;
    ops["J"] = 1.0;
    REQUIRE(ops.isreal());
    ops["J"] = z;
    REQUIRE(!ops.isreal());
    REQUIRE(!isreal(ops));
  }
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}

TEST_CASE("OpSum printing", "[operators]") try {
  using namespace xdiag;

  Op A("SdotS", {0, 1});
  Op B("Sz", int64_t(2));
  Monomial mAB{A, B};

  {
    OpSum ops = 2.0 * A + complex(0.0, 1.5) * B + std::string("J") * mAB;
    ops["J"] = Scalar(0.25);
    std::string s = to_string(ops);
    REQUIRE(s.find("SdotS") != std::string::npos);
    REQUIRE(s.find("Sz") != std::string::npos);
    REQUIRE(s.find("J") != std::string::npos);
    REQUIRE(s.find("0.25") != std::string::npos);

    std::ostringstream oss;
    oss << ops;
    REQUIRE(oss.str() == s);
  }
  {
    // empty and parameter-free OpSums print without throwing
    REQUIRE_NOTHROW(to_string(OpSum()));
    REQUIRE_NOTHROW(to_string(OpSum(A)));
  }
} catch (xdiag::Error e) {
  error_trace(e);
  throw;
}
