// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <cmath>

#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/states/dot.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/error.hpp>

TEST_CASE("dot", "[states]") try {
  using namespace xdiag;

  auto block = Spinhalf(4, 2); // dim 6

  {
    arma::vec ar = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    arma::vec br = {6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
    arma::cx_vec ac(arma::vec{1., 0., 0., 0., 0., 0.},
                    arma::vec{0., 1., 0., 0., 0., 0.});
    arma::cx_vec bc(arma::vec{0., 1., 0., 0., 0., 0.},
                    arma::vec{1., 0., 0., 0., 0., 0.});
    auto vr = State(block, ar);
    auto wr = State(block, br);
    auto vc = State(block, ac);
    auto wc = State(block, bc);
    dot(vr, wr);
    REQUIRE_THROWS(dot(vr, wc));
    REQUIRE_THROWS(dot(vc, wr));
    REQUIRE_THROWS(dot(vc, wc));

    dotC(vr, wr);
    dotC(vr, wc);
    dotC(vc, wr);
    dotC(vc, wc);  
  }

  // --- real dot of single-column states ---
  {
    arma::vec a = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    arma::vec b = {6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
    auto v = State(block, a);
    auto w = State(block, b);
    REQUIRE(std::abs(dot(v, w) - arma::dot(a, b)) < 1e-12);
  }

  // --- complex dotC is conjugate-linear in the first argument ---
  {
    arma::cx_vec a(arma::vec{1., 0., 0., 0., 0., 0.},
                   arma::vec{0., 1., 0., 0., 0., 0.});
    arma::cx_vec b(arma::vec{0., 1., 0., 0., 0., 0.},
                   arma::vec{1., 0., 0., 0., 0., 0.});
    auto v = State(block, a);
    auto w = State(block, b);
    complex expected = arma::cdot(a, b); // conj(a).t() * b
    REQUIRE(std::abs(dotC(v, w) - expected) < 1e-12);
  }

  // --- matrix_dot / matrix_dotC of multi-column states (Gram matrices) ---
  {
    arma::mat M(6, 3, arma::fill::randu);
    auto v = State(block, M);
    arma::mat g = matrix_dot(v, v);
    REQUIRE(arma::norm(g - M.t() * M) < 1e-10);
  }
  {
    arma::cx_mat M(6, 2, arma::fill::randu);
    auto v = State(block, M);
    arma::cx_mat g = matrix_dotC(v, v);
    REQUIRE(arma::norm(g - M.t() * M) < 1e-10);
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("dot errors", "[states]") try {
  using namespace xdiag;

  auto block = Spinhalf(4, 2); // dim 6

  // --- states on different blocks cannot be combined ---
  {
    auto other = Spinhalf(4, 1); // dim 4, different block
    auto v = State(block, true);
    auto w = State(other, true);
    REQUIRE_THROWS_AS(dot(v, w), Error);
    REQUIRE_THROWS_AS(dotC(v, w), Error);
  }

  // --- plain dot requires single-column states ---
  {
    auto v = State(block, true, 2); // 2 columns
    auto w = State(block, true, 2);
    REQUIRE_THROWS_AS(dot(v, w), Error);
    REQUIRE_THROWS_AS(dotC(v, w), Error);
  }

  // --- real dot rejects complex states ---
  {
    auto v = State(block, false); // complex
    auto w = State(block, false);
    REQUIRE_THROWS_AS(dot(v, w), Error);
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
