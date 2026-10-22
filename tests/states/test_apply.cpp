// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/apply.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/dot.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/error.hpp>

TEST_CASE("states_apply", "[states]") try {
  using namespace xdiag;

  auto block = Spinhalf(4, 2);
  auto v = State(block);
  fill(v, RandomState(1234));

  {
    State w = apply(OpSum(), v);
    REQUIRE(!w.isvalid());

    w = apply(OpSum(), State());
    REQUIRE(!w.isvalid());
  }

  // --- apply(Op, State) ---
  {
    State w = apply(Op("SdotS", std::vector<int64_t>{0, 1}), v);
    REQUIRE(w.block() == v.block());
    REQUIRE(isvalid(w));

    Op op = Op("SdotS", std::vector<int64_t>{0, 1});
    apply(op, v, w);
    REQUIRE(isvalid(w));
    apply(Monomial(op), v, w);
    REQUIRE(isvalid(w));
  }

  // Check argument conversions
  {
    Op opr = Op("SdotS", std::vector<int64_t>{0, 1});
    OpSum opc = complex(0, 1) * Op("SdotS", std::vector<int64_t>{0, 1});

    for (int ncols = 1; ncols < 4; ++ncols) {
      auto vr = random_state(block, true, ncols);
      auto vc = random_state(block, false, ncols);
      auto wr = random_state(block, true, ncols);
      auto wc = random_state(block, false, ncols);
      apply(opr, vr, wr);
      apply(opr, vr, wc);
      apply(opr, vc, wr);
      apply(opr, vc, wc);
      apply(opc, vr, wr);
      apply(opc, vr, wc);
      apply(opc, vc, wr);
      apply(opc, vc, wc);
    }
  }
  // --- apply(Monomial, State) ---
  {
    Monomial mono{Op("Sz", int64_t(0)), Op("Sz", int64_t(1))};
    State w = apply(mono, v);
    REQUIRE(isvalid(w));
  }

  // --- apply(OpSum, State), real ---
  {
    OpSum ops;
    ops += Op("SdotS", std::vector<int64_t>{0, 1});
    ops += Op("SdotS", std::vector<int64_t>{1, 2});
    State w = apply(ops, v);
    REQUIRE(isreal(w));
    REQUIRE(isvalid(w));
  }

  // --- apply into a preallocated state ---
  {
    OpSum ops(Op("SdotS", std::vector<int64_t>{0, 1}));
    State w = apply(ops, v);
    State w2(block);
    apply(ops, v, w2);
    REQUIRE(isapprox(w, w2));
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("apply complex operator", "[states]") try {
  using namespace xdiag;

  auto block = Spinhalf(4, 2);
  auto v = State(block); // real start vector
  fill(v, RandomState(9));

  // A complex-coefficient OpSum applied to a real state must produce a complex
  // result (exercises the isreal(ops)==false branch and make_complex).
  OpSum ops;
  ops += complex(0.0, 1.0) * Op("SdotS", std::vector<int64_t>{0, 1});
  State w = apply(ops, v);
  REQUIRE(!isreal(w));
  REQUIRE(isvalid(w));

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("apply multi-column states", "[states]") try {
  using namespace xdiag;

  auto block = Spinhalf(5, 2);
  OpSum opr;
  opr += Op("SdotS", {0, 1});
  opr += 0.5 * Op("SdotS", {1, 3});
  opr += Op("Sz", 2);
  OpSum opc = opr + complex(0.0, 0.7) * Op("SdotS", {2, 4});

  // Reference: the result applied column by column must agree with the
  // multi-column application, for every real/complex combination of operator,
  // input, and output state.
  for (int64_t ncols = 2; ncols < 5; ++ncols) {
    auto vr = random_state(block, true, ncols, 7, false);
    auto vc = random_state(block, false, ncols, 11, false);

    // real operator, real input
    {
      State w = apply(opr, vr);
      REQUIRE(isreal(w));
      REQUIRE(w.ncols() == ncols);
      for (int64_t c = 0; c < ncols; ++c) {
        REQUIRE(isapprox(w.col(c), apply(opr, vr.col(c))));
      }

      // real output preallocated
      State wr(block, true, ncols);
      apply(opr, vr, wr);
      REQUIRE(isreal(wr));
      REQUIRE(isapprox(wr, w));

      // complex output preallocated is narrowed to real
      State wc(block, false, ncols);
      apply(opr, vr, wc);
      REQUIRE(isreal(wc));
      REQUIRE(isapprox(wc, w));
    }

    // real operator, complex input
    {
      State w = apply(opr, vc);
      REQUIRE(!isreal(w));
      for (int64_t c = 0; c < ncols; ++c) {
        REQUIRE(isapprox(w.col(c), apply(opr, vc.col(c))));
      }
      State wr(block, true, ncols);
      apply(opr, vc, wr);
      REQUIRE(!isreal(wr));
      REQUIRE(isapprox(wr, w));
      State wc(block, false, ncols);
      apply(opr, vc, wc);
      REQUIRE(isapprox(wc, w));
    }

    // complex operator, real and complex input
    {
      State w = apply(opc, vr);
      REQUIRE(!isreal(w));
      for (int64_t c = 0; c < ncols; ++c) {
        REQUIRE(isapprox(w.col(c), apply(opc, vr.col(c))));
      }
      State wr(block, true, ncols);
      apply(opc, vr, wr);
      REQUIRE(!isreal(wr));
      REQUIRE(isapprox(wr, w));

      State w2 = apply(opc, vc);
      for (int64_t c = 0; c < ncols; ++c) {
        REQUIRE(isapprox(w2.col(c), apply(opc, vc.col(c))));
      }
      State wc(block, false, ncols);
      apply(opc, vc, wc);
      REQUIRE(isapprox(wc, w2));
    }

    // the multi-column result equals the dense matrix product
    {
      arma::mat H = matrix(opr, block);
      arma::mat V = vr.matrix();
      REQUIRE(arma::norm(apply(opr, vr).matrix() - H * V) < 1e-10);
      arma::cx_mat Hc = matrixC(opc, block);
      arma::cx_mat Vc = vc.matrixC();
      REQUIRE(arma::norm(apply(opc, vc).matrixC() - Hc * Vc) < 1e-10);
    }
  }

  // mismatched column counts between input and output throw
  {
    auto v = random_state(block, true, 2);
    State w(block, true, 3);
    REQUIRE_THROWS_AS(apply(opr, v, w), Error);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
