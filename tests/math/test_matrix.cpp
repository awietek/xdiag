// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/math/matrix.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

TEST_CASE("test_matrix", "[complex]") {
  Log("Testing Matrix");

  arma::mat rm = {{1.0, 2.0}, {3.0, 4.0}};
  arma::cx_mat cm = {{complex(1.0, 1.0), complex(2.0, 0.0)},
                     {complex(0.0, -1.0), complex(4.0, 2.0)}};

  // --- Construction and type query ---
  {
    Matrix mr(rm);
    REQUIRE(mr.isreal());
    REQUIRE(mr.is<arma::mat>());
    REQUIRE(!mr.is<arma::cx_mat>());
    REQUIRE(mr.n_rows() == 2);
    REQUIRE(mr.n_cols() == 2);

    Matrix mc(cm);
    REQUIRE(!mc.isreal());
    REQUIRE(mc.is<arma::cx_mat>());
    REQUIRE(!mc.is<arma::mat>());
    REQUIRE(mc.n_rows() == 2);
    REQUIRE(mc.n_cols() == 2);
  }

  // --- as<T> ---
  {
    Matrix mr(rm);
    REQUIRE(arma::approx_equal(mr.as<arma::mat>(), rm, "absdiff", 1e-14));
    arma::cx_mat promoted = mr.as<arma::cx_mat>();
    REQUIRE(arma::norm(arma::imag(promoted), "inf") == Approx(0.0));

    Matrix mc(cm);
    // REQUIRE_THROWS(mc.as<arma::mat>());
    REQUIRE(arma::approx_equal(mc.as<arma::cx_mat>(), cm, "absdiff", 1e-14));
  }

  // --- real(), imag() ---
  {
    Matrix mc(cm);
    REQUIRE(mc.real()(0, 0) == Approx(1.0));
    REQUIRE(mc.imag()(0, 0) == Approx(1.0));
    REQUIRE(mc.real()(1, 0) == Approx(0.0));
    REQUIRE(mc.imag()(1, 0) == Approx(-1.0));

    Matrix mr(rm);
    REQUIRE(arma::norm(mr.imag(), "inf") == Approx(0.0));
  }

  // --- hc: conjugate transpose ---
  {
    Matrix mr(rm);
    arma::mat expected_r = arma::trans(rm);
    REQUIRE(arma::approx_equal(mr.hc().as<arma::mat>(), expected_r, "absdiff",
                               1e-14));
    REQUIRE(mr.hc().isreal());

    Matrix mc(cm);
    arma::cx_mat expected_c = arma::trans(cm);
    REQUIRE(arma::approx_equal(mc.hc().as<arma::cx_mat>(), expected_c,
                               "absdiff", 1e-14));
    // free function
    REQUIRE(arma::approx_equal(hc(mc).as<arma::cx_mat>(), expected_c, "absdiff",
                               1e-14));
  }

  // --- to_real ---
  {
    Matrix mr(rm);
    REQUIRE(mr.to_real().isreal());

    arma::cx_mat near_real =
        arma::cx_mat(rm, arma::mat(2, 2, arma::fill::zeros));
    Matrix mnr(near_real);
    REQUIRE(mnr.to_real(1e-12).isreal());
    REQUIRE(arma::approx_equal(mnr.to_real(1e-12).as<arma::mat>(), rm,
                               "absdiff", 1e-14));

    Matrix mc(cm);
    // REQUIRE_THROWS(mc.to_real());
  }

  // --- Linear combination (scalar ops) ---
  {
    arma::mat rm2 = rm * 2.0;
    arma::mat rm3 = rm * 3.0;
    Matrix a(rm), b(rm2);
    Matrix sum = a + b;
    REQUIRE(sum.isreal());
    REQUIRE(arma::approx_equal(sum.as<arma::mat>(), rm3, "absdiff", 1e-14));

    Matrix diff = b - a;
    REQUIRE(arma::approx_equal(diff.as<arma::mat>(), rm, "absdiff", 1e-14));

    Matrix scaled = a * Scalar(3.0);
    REQUIRE(arma::approx_equal(scaled.as<arma::mat>(), rm3, "absdiff", 1e-14));

    Matrix neg = -a;
    arma::mat nrm = -rm;
    REQUIRE(arma::approx_equal(neg.as<arma::mat>(), nrm, "absdiff", 1e-14));
  }

  // --- Matrix * Matrix ---
  {
    Matrix a(rm), b(rm);
    Matrix prod = a * b;
    REQUIRE(prod.isreal());
    REQUIRE(
        arma::approx_equal(prod.as<arma::mat>(), rm * rm, "absdiff", 1e-14));

    // real * complex widens
    Matrix c(cm);
    Matrix mixed = a * c;
    REQUIRE(!mixed.isreal());
    arma::cx_mat expected =
        arma::cx_mat(rm, arma::mat(2, 2, arma::fill::zeros)) * cm;
    REQUIRE(arma::approx_equal(mixed.as<arma::cx_mat>(), expected, "absdiff",
                               1e-14));
  }

  // --- Matrix * Vector ---
  {
    arma::vec rv = {1.0, 2.0};
    Matrix mr(rm);
    Vector vr(rv);
    Vector result = mr * vr;
    REQUIRE(result.isreal());
    REQUIRE(
        arma::approx_equal(result.as<arma::vec>(), rm * rv, "absdiff", 1e-14));

    // complex matrix * real vector widens
    Matrix mc(cm);
    Vector vc_result = mc * vr;
    REQUIRE(!vc_result.isreal());
    arma::cx_vec expected =
        cm * arma::cx_vec(rv, arma::vec(2, arma::fill::zeros));
    REQUIRE(arma::approx_equal(vc_result.as<arma::cx_vec>(), expected,
                               "absdiff", 1e-14));
  }

  // --- Promotion on scalar multiply ---
  {
    Matrix mr(rm);
    Matrix mc = mr * Scalar(complex(0.0, 1.0));
    REQUIRE(!mc.isreal());
    arma::cx_mat expected(arma::mat(2, 2, arma::fill::zeros), rm);
    REQUIRE(
        arma::approx_equal(mc.as<arma::cx_mat>(), expected, "absdiff", 1e-14));
  }

  // --- isapprox ---
  {
    Matrix a(rm);
    arma::mat rm_perturbed = rm + arma::mat(2, 2, arma::fill::ones) * 1e-14;
    Matrix b(rm_perturbed);
    REQUIRE(isapprox(a, b));
    arma::mat rm2 = rm * 2.0;
    REQUIRE(!isapprox(a, Matrix(rm2)));
  }

  // --- compound assignment ---
  {
    arma::mat rm2 = rm * 2.0;
    Matrix a(rm);
    a += Matrix(rm);
    REQUIRE(arma::approx_equal(a.as<arma::mat>(), rm2, "absdiff", 1e-14));
    a -= Matrix(rm);
    REQUIRE(arma::approx_equal(a.as<arma::mat>(), rm, "absdiff", 1e-14));
    a *= Scalar(2.0);
    REQUIRE(arma::approx_equal(a.as<arma::mat>(), rm2, "absdiff", 1e-14));
    a /= Scalar(2.0);
    REQUIRE(arma::approx_equal(a.as<arma::mat>(), rm, "absdiff", 1e-14));
  }

  // --- n_rows / n_cols with non-square matrix ---
  {
    arma::mat rect(3, 5, arma::fill::zeros);
    Matrix m(rect);
    REQUIRE(m.n_rows() == 3);
    REQUIRE(m.n_cols() == 5);
  }
}
