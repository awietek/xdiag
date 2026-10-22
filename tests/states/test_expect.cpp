// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <tuple>

#include <xdiag/algebra/symmetrize.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/math/isapprox.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/apply.hpp>
#include <xdiag/states/create_state.hpp>
#include <xdiag/states/dot.hpp>
#include <xdiag/states/expect.hpp>
#include <xdiag/states/inner.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/utils/error.hpp>

using namespace xdiag;

// Naive expectation per site: <psi| build(i) |psi> for an OpSum builder. On a
// symmetry-adapted block, build(i) must already be symmetrized so the bare
// inner / apply stays inside the block.
template <typename Builder>
static arma::vec naive_expect(State const &psi, Builder build) {
  int64_t n = psi.nsites();
  arma::vec r(n);
  for (int64_t i = 0; i < n; ++i) {
    r(i) = inner(psi, build(i), psi);
  }
  return r;
}

TEST_CASE("expect", "[states]") try {
  // ===== Spin-1/2 Heisenberg ring: <Sz_i> =====
  {
    int64_t n = 4;
    OpSum H;
    for (int64_t i = 0; i < n; ++i) {
      H += Op("SdotS", {i, (i + 1) % n});
    }

    // No permutation symmetry (nup sector only): bare ops act in-block.
    {
      auto block = Spinhalf(n, n / 2);
      State psi = std::get<1>(eig0(H, block));
      arma::vec ex = expect(psi, "Sz");
      arma::vec naive =
          naive_expect(psi, [&](int64_t i) { return OpSum(Op("Sz", i)); });
      REQUIRE(isapprox(ex, naive, 1e-10, 1e-10));
    }

    // Translation symmetry (k = 0): expect must symmetrize internally; the
    // naive uses the same group-symmetrized operator.
    {
      auto group = cyclic_group(n);
      auto block = Spinhalf(n, n / 2, cyclic_group_irrep(n, 0));
      State psi = std::get<1>(eig0(H, block));
      arma::vec ex = expect(psi, "Sz");
      arma::vec naive = naive_expect(
          psi, [&](int64_t i) { return symmetrize(Op("Sz", i), group); });
      REQUIRE(isapprox(ex, naive, 1e-10, 1e-10));
    }
  }

  // ===== Bose-Hubbard ring: <N_i> =====
  {
    int64_t n = 4, d = 3, number = 3;
    OpSum H;
    for (int64_t i = 0; i < n; ++i) {
      H += Op("Hop", {i, (i + 1) % n});
    }
    H += 2.0 * Op("HubbardU");

    {
      auto block = Boson(n, d, number);
      State psi = std::get<1>(eig0(H, block));
      arma::vec ex = expect(psi, "N");
      arma::vec naive =
          naive_expect(psi, [&](int64_t i) { return OpSum(Op("N", i)); });
      REQUIRE(isapprox(ex, naive, 1e-10, 1e-10));
    }

    {
      auto group = cyclic_group(n);
      auto block = Boson(n, d, number, cyclic_group_irrep(n, 0));
      State psi = std::get<1>(eig0(H, block));
      arma::vec ex = expect(psi, "N");
      arma::vec naive = naive_expect(
          psi, [&](int64_t i) { return symmetrize(Op("N", i), group); });
      REQUIRE(isapprox(ex, naive, 1e-10, 1e-10));
    }
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}

TEST_CASE("inner overloads", "[states]") try {
  using namespace xdiag;

  auto block = Spinhalf(4, 2);
  auto v = random_state(block, true, 1, 3);
  auto w = random_state(block, true, 1, 5);
  auto vc = random_state(block, false, 1, 3);

  Op sz0("Sz", 0);
  Op sz1("Sz", 1);
  Monomial mono{sz0, sz1};
  OpSum ops(mono);
  OpSum opsc = complex(0.0, 1.0) * mono;

  // --- expectation values: Op, Monomial, OpSum agree ---
  {
    double eo = inner(sz0, v);
    double em = inner(mono, v);
    double es = inner(ops, v);
    REQUIRE(std::abs(em - es) < 1e-12);
    REQUIRE(std::abs(eo - dot(v, apply(sz0, v))) < 1e-12);
    REQUIRE(std::abs(em - dot(v, apply(mono, v))) < 1e-12);

    REQUIRE(std::abs(innerC(sz0, v) - complex(eo)) < 1e-12);
    REQUIRE(std::abs(innerC(mono, v) - complex(em)) < 1e-12);
    REQUIRE(std::abs(innerC(ops, v) - complex(es)) < 1e-12);
  }

  // --- matrix elements <v|A|w> ---
  {
    double mo = inner(v, sz0, w);
    double mm = inner(v, mono, w);
    double ms = inner(v, ops, w);
    REQUIRE(std::abs(mm - ms) < 1e-12);
    REQUIRE(std::abs(mo - dot(v, apply(sz0, w))) < 1e-12);
    REQUIRE(std::abs(innerC(v, sz0, w) - complex(mo)) < 1e-12);
    REQUIRE(std::abs(innerC(v, mono, w) - complex(mm)) < 1e-12);
    REQUIRE(std::abs(innerC(v, ops, w) - complex(ms)) < 1e-12);
  }

  // --- complex states and operators need innerC ---
  {
    REQUIRE_THROWS_AS(inner(sz0, vc), Error);
    REQUIRE_THROWS_AS(inner(mono, vc), Error);
    REQUIRE_THROWS_AS(inner(opsc, v), Error);
    REQUIRE_THROWS_AS(inner(v, sz0, vc), Error);
    REQUIRE_THROWS_AS(inner(vc, mono, v), Error);

    complex ec = innerC(opsc, v);
    REQUIRE(std::abs(ec - complex(0.0, 1.0) * inner(mono, v)) < 1e-12);
    complex em = innerC(mono, vc);
    REQUIRE(std::abs(em - dotC(vc, apply(mono, vc))) < 1e-12);
    complex eo = innerC(vc, sz0, v);
    REQUIRE(std::abs(eo - dotC(vc, apply(sz0, v))) < 1e-12);
    complex emm = innerC(v, mono, vc);
    REQUIRE(std::abs(emm - dotC(v, apply(mono, vc))) < 1e-12);
  }
} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
