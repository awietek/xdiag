// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <sstream>

#include <xdiag/math/vector.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

TEST_CASE("test_vector", "[complex]") {
  Log("Testing Vector");

  arma::vec rv = {1.0, 2.0, 3.0};
  arma::cx_vec cv = {complex(1.0, 1.0), complex(2.0, -1.0), complex(0.0, 3.0)};

  // --- Construction and type query ---
  {
    Vector vr(rv);
    REQUIRE(vr.isreal());
    REQUIRE(vr.is<arma::vec>());
    REQUIRE(!vr.is<arma::cx_vec>());
    REQUIRE(vr.size() == 3);

    Vector vc(cv);
    REQUIRE(!vc.isreal());
    REQUIRE(vc.is<arma::cx_vec>());
    REQUIRE(!vc.is<arma::vec>());
    REQUIRE(vc.size() == 3);
  }

  // --- as<T>: real always works; complex only from complex ---
  {
    Vector vr(rv);
    REQUIRE(arma::approx_equal(vr.as<arma::vec>(), rv, "absdiff", 1e-14));
    // as<arma::cx_vec> widens silently
    arma::cx_vec promoted = vr.as<arma::cx_vec>();
    REQUIRE(std::imag(promoted(0)) == Approx(0.0));
    REQUIRE(std::imag(promoted(1)) == Approx(0.0));
    REQUIRE(std::imag(promoted(2)) == Approx(0.0));

    Vector vc(cv);
    // REQUIRE_THROWS(vc.as<arma::vec>());
    REQUIRE(arma::approx_equal(vc.as<arma::cx_vec>(), cv, "absdiff", 1e-14));
  }

  // --- real(), imag() ---
  {
    Vector vc(cv);
    arma::vec r = vc.real();
    arma::vec i = vc.imag();
    REQUIRE(r(0) == Approx(1.0));
    REQUIRE(i(0) == Approx(1.0));
    REQUIRE(r(2) == Approx(0.0));
    REQUIRE(i(2) == Approx(3.0));

    // Vector vr(rv);
    // REQUIRE(arma::norm(vr.imag()) == Approx(0.0));
  }

  // --- conj ---
  {
    Vector vr(rv);
    Vector cr = vr.conj();
    REQUIRE(cr.isreal());
    REQUIRE(arma::approx_equal(cr.as<arma::vec>(), rv, "absdiff", 1e-14));

    Vector vc(cv);
    Vector cc = vc.conj();
    REQUIRE(!cc.isreal());
    arma::cx_vec expected = arma::conj(cv);
    REQUIRE(
        arma::approx_equal(cc.as<arma::cx_vec>(), expected, "absdiff", 1e-14));
    // free function
    REQUIRE(arma::approx_equal(conj(vc).as<arma::cx_vec>(), expected, "absdiff",
                               1e-14));
  }

  // --- to_real ---
  {
    Vector vr(rv);
    REQUIRE(vr.to_real().isreal());

    arma::cx_vec near_real = {complex(1.0, 1e-15), complex(2.0, 0.0),
                              complex(3.0, -1e-15)};
    Vector vnr(near_real);
    REQUIRE(vnr.to_real(1e-12).isreal());

    Vector vc(cv);
    // REQUIRE_THROWS(vc.to_real());
  }

  // --- Vector space: same-type arithmetic ---
  {
    Vector a(rv), b(rv * 2.0);
    Vector sum = a + b;
    REQUIRE(sum.isreal());
    REQUIRE(
        arma::approx_equal(sum.as<arma::vec>(), rv * 3.0, "absdiff", 1e-14));

    Vector diff = b - a;
    REQUIRE(arma::approx_equal(diff.as<arma::vec>(), rv, "absdiff", 1e-14));

    Vector scaled = a * Scalar(3.0);
    REQUIRE(
        arma::approx_equal(scaled.as<arma::vec>(), rv * 3.0, "absdiff", 1e-14));

    Vector div = b / Scalar(2.0);
    REQUIRE(arma::approx_equal(div.as<arma::vec>(), rv, "absdiff", 1e-14));

    Vector neg = -a;
    REQUIRE(arma::approx_equal(neg.as<arma::vec>(), -rv, "absdiff", 1e-14));
  }

  // --- Promotion: real + complex widens ---
  {
    Vector vr(rv);
    Vector vc(cv);
    Vector sum = vr + vc;
    REQUIRE(!sum.isreal());
    arma::cx_vec expected =
        arma::cx_vec(rv, arma::vec(3, arma::fill::zeros)) + cv;
    REQUIRE(
        arma::approx_equal(sum.as<arma::cx_vec>(), expected, "absdiff", 1e-14));
  }

  // --- Scalar promotion on multiply ---
  {
    Vector vr(rv);
    Vector vc = vr * Scalar(complex(0.0, 1.0));
    REQUIRE(!vc.isreal());
    arma::cx_vec expected(arma::vec(3, arma::fill::zeros), rv);
    REQUIRE(
        arma::approx_equal(vc.as<arma::cx_vec>(), expected, "absdiff", 1e-14));
  }

  // --- isapprox ---
  {
    Vector a(rv);
    Vector b(rv + 1e-14);
    REQUIRE(isapprox(a, b));
    REQUIRE(!isapprox(a, Vector(rv * 2.0)));
  }

  // --- compound assignment ---
  {
    Vector a(rv);
    a += Vector(rv);
    REQUIRE(arma::approx_equal(a.as<arma::vec>(), rv * 2.0, "absdiff", 1e-14));
    a -= Vector(rv);
    REQUIRE(arma::approx_equal(a.as<arma::vec>(), rv, "absdiff", 1e-14));
    a *= Scalar(2.0);
    REQUIRE(arma::approx_equal(a.as<arma::vec>(), rv * 2.0, "absdiff", 1e-14));
    a /= Scalar(2.0);
    REQUIRE(arma::approx_equal(a.as<arma::vec>(), rv, "absdiff", 1e-14));
  }
}

TEST_CASE("test_vector_extended", "[complex]") {
  Log("Testing Vector extended API");

  arma::vec rv = {1.0, 2.0, 3.0};
  arma::cx_vec cv = {complex(1.0, 1.0), complex(2.0, -1.0), complex(0.0, 3.0)};
  arma::cx_vec rvc = arma::conv_to<arma::cx_vec>::from(rv);
  complex z(0.5, -1.5);

  // --- construction from std::vector ---
  {
    Vector vr(std::vector<double>{1.0, 2.0, 3.0});
    REQUIRE(vr.isreal());
    REQUIRE(arma::approx_equal(vr.as<arma::vec>(), rv, "absdiff", 1e-14));

    Vector vc(std::vector<complex>{complex(1.0, 1.0), complex(2.0, -1.0),
                                   complex(0.0, 3.0)});
    REQUIRE(!vc.isreal());
    REQUIRE(arma::approx_equal(vc.as<arma::cx_vec>(), cv, "absdiff", 1e-14));
  }

  // --- equality and inequality ---
  {
    Vector a(rv), b(rv), c(cv);
    arma::vec rv2 = rv * 2.0;
    Vector d(rv2);
    REQUIRE(a == b);
    REQUIRE(!(a != b));
    REQUIRE(a != c);
    REQUIRE(a != d);
    REQUIRE(c == Vector(cv));
  }

  // --- free real / imag ---
  {
    Vector vc(cv);
    REQUIRE(arma::approx_equal(real(vc), arma::real(cv), "absdiff", 1e-14));
    REQUIRE(arma::approx_equal(imag(vc), arma::imag(cv), "absdiff", 1e-14));
    Vector vr(rv);
    REQUIRE(arma::approx_equal(real(vr), rv, "absdiff", 1e-14));
    REQUIRE(arma::norm(imag(vr)) == Approx(0.0));
    REQUIRE(imag(vr).n_elem == 3);
    REQUIRE(isreal(vr));
    REQUIRE(!isreal(vc));
  }

  // --- mixed compound assignment widens the left operand ---
  {
    Vector a(rv);
    a += Vector(cv);
    REQUIRE(!a.isreal());
    arma::cx_vec e1 = rvc + cv;
    REQUIRE(arma::approx_equal(a.as<arma::cx_vec>(), e1, "absdiff", 1e-14));

    Vector b(rv);
    b -= Vector(cv);
    REQUIRE(!b.isreal());
    arma::cx_vec e2 = rvc - cv;
    REQUIRE(arma::approx_equal(b.as<arma::cx_vec>(), e2, "absdiff", 1e-14));

    Vector c(cv);
    c += Vector(rv);
    arma::cx_vec e3 = cv + rvc;
    REQUIRE(arma::approx_equal(c.as<arma::cx_vec>(), e3, "absdiff", 1e-14));

    Vector d(cv);
    d -= Vector(rv);
    arma::cx_vec e4 = cv - rvc;
    REQUIRE(arma::approx_equal(d.as<arma::cx_vec>(), e4, "absdiff", 1e-14));

    Vector f(cv);
    f += Vector(cv);
    arma::cx_vec e5 = cv + cv;
    REQUIRE(arma::approx_equal(f.as<arma::cx_vec>(), e5, "absdiff", 1e-14));
    f -= Vector(cv);
    REQUIRE(arma::approx_equal(f.as<arma::cx_vec>(), cv, "absdiff", 1e-14));
  }

  // --- scalar compound assignment in every combination ---
  {
    Vector a(rv);
    a *= Scalar(z);
    REQUIRE(!a.isreal());
    arma::cx_vec e1 = rvc * z;
    REQUIRE(arma::approx_equal(a.as<arma::cx_vec>(), e1, "absdiff", 1e-14));

    Vector b(cv);
    b *= Scalar(2.0);
    arma::cx_vec e2 = cv * 2.0;
    REQUIRE(arma::approx_equal(b.as<arma::cx_vec>(), e2, "absdiff", 1e-14));
    b *= Scalar(z);
    arma::cx_vec e3 = cv * 2.0 * z;
    REQUIRE(arma::approx_equal(b.as<arma::cx_vec>(), e3, "absdiff", 1e-14));

    Vector c(rv);
    c /= Scalar(z);
    REQUIRE(!c.isreal());
    arma::cx_vec e4 = rvc / z;
    REQUIRE(arma::approx_equal(c.as<arma::cx_vec>(), e4, "absdiff", 1e-14));

    Vector d(cv);
    d /= Scalar(2.0);
    arma::cx_vec e5 = cv / 2.0;
    REQUIRE(arma::approx_equal(d.as<arma::cx_vec>(), e5, "absdiff", 1e-14));
    d /= Scalar(z);
    arma::cx_vec e6 = cv / 2.0 / z;
    REQUIRE(arma::approx_equal(d.as<arma::cx_vec>(), e6, "absdiff", 1e-14));

    // non-mutating forms
    Vector r(rv);
    REQUIRE(arma::approx_equal((r / Scalar(z)).as<arma::cx_vec>(), e4,
                               "absdiff", 1e-14));
    REQUIRE(arma::approx_equal((r * Scalar(z)).as<arma::cx_vec>(), e1,
                               "absdiff", 1e-14));
  }

  // --- unary minus and mixed binary operators ---
  {
    Vector vc(cv);
    arma::cx_vec ncv = -cv;
    REQUIRE(arma::approx_equal((-vc).as<arma::cx_vec>(), ncv, "absdiff",
                               1e-14));
    Vector vr(rv);
    arma::cx_vec e1 = rvc - cv;
    REQUIRE(arma::approx_equal((vr - vc).as<arma::cx_vec>(), e1, "absdiff",
                               1e-14));
    arma::cx_vec e2 = cv - rvc;
    REQUIRE(arma::approx_equal((vc - vr).as<arma::cx_vec>(), e2, "absdiff",
                               1e-14));
    arma::cx_vec e3 = cv + rvc;
    REQUIRE(arma::approx_equal((vc + vr).as<arma::cx_vec>(), e3, "absdiff",
                               1e-14));
  }

  // --- printing ---
  {
    Vector vr(rv);
    std::string s = to_string(vr);
    REQUIRE(!s.empty());
    std::ostringstream oss;
    oss << vr;
    REQUIRE(oss.str() == s);
    REQUIRE(!to_string(Vector(cv)).empty());
  }

  // --- as<arma::vec> and to_real on genuinely complex data throw ---
  {
    Vector vc(cv);
    REQUIRE_THROWS(vc.as<arma::vec>());
    REQUIRE_THROWS(vc.to_real());
  }
}
