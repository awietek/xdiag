// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/kernels/apply.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

// The apply overloads taking an Op or a Monomial together with a Block and an
// armadillo vector are public API that nothing inside the library reaches,
// since states/apply.cpp promotes to an OpSum first. A non-hermitian operator
// is used on purpose: for a hermitian one an adjoint applied by mistake would
// give the same answer.
TEST_CASE("apply_op_monomial_overloads", "[kernels]") {
  Log("apply_op_monomial_overloads test ...");

  int64_t nsites = 6;
  int64_t nup = 2;
  Spinhalf block_in(nsites, nup);
  Spinhalf block_out(nsites, nup + 1);
  Op op("S+", std::vector<int64_t>{2});
  Monomial mono(op);

  arma::mat ref = matrix(op, block_in, block_out);
  REQUIRE(arma::norm(ref) > 1e-8); // the operator must not be trivial
  arma::cx_mat refc = arma::conv_to<arma::cx_mat>::from(ref);

  int64_t din = dim(block_in), dout = dim(block_out);
  arma::vec v(din, arma::fill::randn);
  arma::cx_vec vc(din, arma::fill::randn);
  arma::mat m(din, 3, arma::fill::randn);
  arma::cx_mat mc(din, 3, arma::fill::randn);

  SECTION("Op") {
    arma::vec w(dout, arma::fill::zeros);
    apply(op, block_in, v, block_out, w);
    REQUIRE(arma::norm(w - ref * v) < 1e-12);

    arma::cx_vec wc(dout, arma::fill::zeros);
    apply(op, block_in, vc, block_out, wc);
    REQUIRE(arma::norm(wc - refc * vc) < 1e-12);

    arma::mat wm(dout, 3, arma::fill::zeros);
    apply(op, block_in, m, block_out, wm);
    REQUIRE(arma::norm(wm - ref * m) < 1e-12);

    arma::cx_mat wmc(dout, 3, arma::fill::zeros);
    apply(op, block_in, mc, block_out, wmc);
    REQUIRE(arma::norm(wmc - refc * mc) < 1e-12);
  }

  SECTION("Monomial") {
    arma::vec w(dout, arma::fill::zeros);
    apply(mono, block_in, v, block_out, w);
    REQUIRE(arma::norm(w - ref * v) < 1e-12);

    arma::cx_vec wc(dout, arma::fill::zeros);
    apply(mono, block_in, vc, block_out, wc);
    REQUIRE(arma::norm(wc - refc * vc) < 1e-12);

    arma::mat wm(dout, 3, arma::fill::zeros);
    apply(mono, block_in, m, block_out, wm);
    REQUIRE(arma::norm(wm - ref * m) < 1e-12);

    arma::cx_mat wmc(dout, 3, arma::fill::zeros);
    apply(mono, block_in, mc, block_out, wmc);
    REQUIRE(arma::norm(wmc - refc * mc) < 1e-12);
  }
}
