// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <sstream>

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

TEST_CASE("test_matrix_extended", "[complex]") {
  Log("Testing Matrix extended API");

  arma::mat rm = {{1.0, 2.0}, {3.0, 4.0}};
  arma::cx_mat cm = {{complex(1.0, 1.0), complex(2.0, 0.0)},
                     {complex(0.0, -1.0), complex(4.0, 2.0)}};
  arma::cx_mat rmc = arma::conv_to<arma::cx_mat>::from(rm);
  complex z(0.5, -1.5);

  // --- equality and inequality ---
  {
    Matrix a(rm), b(rm), c(cm);
    arma::mat rm2 = rm * 2.0;
    Matrix d(rm2);
    REQUIRE(a == b);
    REQUIRE(!(a != b));
    REQUIRE(a != c);
    REQUIRE(a != d);
    REQUIRE(c == Matrix(cm));
  }

  // --- free real / imag / conj / transpose ---
  {
    Matrix mc(cm);
    REQUIRE(arma::approx_equal(real(mc), arma::real(cm), "absdiff", 1e-14));
    REQUIRE(arma::approx_equal(imag(mc), arma::imag(cm), "absdiff", 1e-14));
    Matrix mr(rm);
    REQUIRE(arma::approx_equal(real(mr), rm, "absdiff", 1e-14));
    REQUIRE(arma::norm(imag(mr), "inf") == Approx(0.0));
    REQUIRE(isreal(mr));
    REQUIRE(!isreal(mc));

    arma::cx_mat cmconj = arma::conj(cm);
    REQUIRE(arma::approx_equal(conj(mc).as<arma::cx_mat>(), cmconj, "absdiff",
                               1e-14));
    REQUIRE(arma::approx_equal(mc.conj().as<arma::cx_mat>(), cmconj, "absdiff",
                               1e-14));
    REQUIRE(conj(mr).isreal());
    REQUIRE(arma::approx_equal(conj(mr).as<arma::mat>(), rm, "absdiff", 1e-14));

    arma::cx_mat cmt = cm.st();
    REQUIRE(arma::approx_equal(transpose(mc).as<arma::cx_mat>(), cmt, "absdiff",
                               1e-14));
    REQUIRE(arma::approx_equal(mc.transpose().as<arma::cx_mat>(), cmt,
                               "absdiff", 1e-14));
    arma::mat rmt = rm.t();
    REQUIRE(transpose(mr).isreal());
    REQUIRE(arma::approx_equal(transpose(mr).as<arma::mat>(), rmt, "absdiff",
                               1e-14));
    // transpose differs from hc for complex matrices
    REQUIRE(transpose(mc) != hc(mc));
  }

  // --- division by scalar ---
  {
    Matrix mr(rm);
    arma::mat e1 = rm / 2.0;
    REQUIRE(arma::approx_equal((mr / Scalar(2.0)).as<arma::mat>(), e1,
                               "absdiff", 1e-14));
    arma::cx_mat e2 = rmc / z;
    Matrix q = mr / Scalar(z);
    REQUIRE(!q.isreal());
    REQUIRE(arma::approx_equal(q.as<arma::cx_mat>(), e2, "absdiff", 1e-14));

    Matrix mc(cm);
    arma::cx_mat e3 = cm / 2.0;
    REQUIRE(arma::approx_equal((mc / Scalar(2.0)).as<arma::cx_mat>(), e3,
                               "absdiff", 1e-14));
    arma::cx_mat e4 = cm / z;
    REQUIRE(arma::approx_equal((mc / Scalar(z)).as<arma::cx_mat>(), e4,
                               "absdiff", 1e-14));

    Matrix a(rm);
    a /= Scalar(z);
    REQUIRE(arma::approx_equal(a.as<arma::cx_mat>(), e2, "absdiff", 1e-14));
    Matrix b(cm);
    b /= Scalar(2.0);
    REQUIRE(arma::approx_equal(b.as<arma::cx_mat>(), e3, "absdiff", 1e-14));
    b /= Scalar(z);
    arma::cx_mat e5 = cm / 2.0 / z;
    REQUIRE(arma::approx_equal(b.as<arma::cx_mat>(), e5, "absdiff", 1e-14));
  }

  // --- mixed compound assignment widens the left operand ---
  {
    Matrix a(rm);
    a += Matrix(cm);
    REQUIRE(!a.isreal());
    arma::cx_mat e1 = rmc + cm;
    REQUIRE(arma::approx_equal(a.as<arma::cx_mat>(), e1, "absdiff", 1e-14));

    Matrix b(rm);
    b -= Matrix(cm);
    arma::cx_mat e2 = rmc - cm;
    REQUIRE(arma::approx_equal(b.as<arma::cx_mat>(), e2, "absdiff", 1e-14));

    Matrix c(cm);
    c += Matrix(rm);
    arma::cx_mat e3 = cm + rmc;
    REQUIRE(arma::approx_equal(c.as<arma::cx_mat>(), e3, "absdiff", 1e-14));

    Matrix d(cm);
    d -= Matrix(rm);
    arma::cx_mat e4 = cm - rmc;
    REQUIRE(arma::approx_equal(d.as<arma::cx_mat>(), e4, "absdiff", 1e-14));

    Matrix f(cm);
    f -= Matrix(cm);
    REQUIRE(arma::norm(f.as<arma::cx_mat>(), "inf") == Approx(0.0));

    Matrix g(cm);
    g *= Scalar(z);
    arma::cx_mat e5 = cm * z;
    REQUIRE(arma::approx_equal(g.as<arma::cx_mat>(), e5, "absdiff", 1e-14));
  }

  // --- unary minus and mixed binary operators on complex ---
  {
    Matrix mc(cm);
    arma::cx_mat ncm = -cm;
    REQUIRE(arma::approx_equal((-mc).as<arma::cx_mat>(), ncm, "absdiff",
                               1e-14));
    Matrix mr(rm);
    arma::cx_mat e1 = rmc - cm;
    REQUIRE(arma::approx_equal((mr - mc).as<arma::cx_mat>(), e1, "absdiff",
                               1e-14));
    arma::cx_mat e2 = cm + rmc;
    REQUIRE(arma::approx_equal((mc + mr).as<arma::cx_mat>(), e2, "absdiff",
                               1e-14));
    arma::cx_mat e3 = cm * rmc;
    REQUIRE(arma::approx_equal((mc * mr).as<arma::cx_mat>(), e3, "absdiff",
                               1e-14));
    arma::cx_mat e4 = cm * cm;
    REQUIRE(arma::approx_equal((mc * mc).as<arma::cx_mat>(), e4, "absdiff",
                               1e-14));
  }

  // --- Matrix * Vector with complex vector ---
  {
    arma::cx_vec cv = {complex(1.0, 2.0), complex(-1.0, 0.5)};
    Matrix mr(rm);
    Matrix mc(cm);
    Vector vc(cv);
    arma::cx_vec e1 = rmc * cv;
    REQUIRE(arma::approx_equal((mr * vc).as<arma::cx_vec>(), e1, "absdiff",
                               1e-14));
    arma::cx_vec e2 = cm * cv;
    REQUIRE(arma::approx_equal((mc * vc).as<arma::cx_vec>(), e2, "absdiff",
                               1e-14));
  }

  // --- printing ---
  {
    Matrix mr(rm);
    std::string s = to_string(mr);
    REQUIRE(!s.empty());
    std::ostringstream oss;
    oss << mr;
    REQUIRE(oss.str() == s);
    REQUIRE(!to_string(Matrix(cm)).empty());
  }

  // --- as<arma::mat> and to_real on genuinely complex data throw ---
  {
    Matrix mc(cm);
    REQUIRE_THROWS(mc.as<arma::mat>());
    REQUIRE_THROWS(mc.to_real());
  }
}
