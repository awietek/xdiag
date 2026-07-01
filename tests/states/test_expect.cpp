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
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
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
}
