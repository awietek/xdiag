// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

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
