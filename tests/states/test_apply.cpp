// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/blocks/spinhalf.hpp>
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
      auto vr = random_state(block, true);
      auto vc = random_state(block, false);
      auto wr = random_state(block, true);
      auto wc = random_state(block, false);
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
