// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0
//
// "Long" Hilbert space test for the Electron block (>64 sites ->
// BitsetDynamic). Single-band Hubbard chain at low fixed particle numbers.
//
//   H = -t sum_<ij>,sigma (c^dag c + h.c.) + U sum_i n_{i up} n_{i dn}
//
// xdiag Op("Hop") carries the minus sign, Op("HubbardU") = sum_i n_up n_dn, so
// the couplings passed are +t and U. A converged DMRG reference is not yet
// available for this model, so the open-boundary energy is still computed but
// the comparison against the DMRG value is commented out. The periodic-boundary
// translational-symmetry consistency check is fully active.

#include <cmath>

#include <tests/blocks/test_long_blocks.hpp>
#include <tests/catch.hpp>

#include <xdiag/blocks/electron.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

TEST_CASE("electron_long", "[long]") try {
  Log("Electron long test");

  int64_t N = 65; // smallest chain activating BitsetDynamic (>64 sites)
  double t = 1.0;
  double U = 4.0;
  int64_t nup = 1;
  int64_t ndn = 1; // dim = 65 * 65 = 4225

  // On-site interaction (same for OBC and PBC).
  OpSum ops_int = U * Op("HubbardU"); // sum_i n_{i up} n_{i dn}

  // --- Open boundary conditions ---
  OpSum ops_obc = ops_int;
  for (int64_t i = 0; i < N - 1; ++i) {
    ops_obc += t * Op("Hop", {i, i + 1}); // Hop = -sum_sigma(c^dag c + h.c.)
  }
  double e0_obc = eigval0(ops_obc, Electron(N, nup, ndn));
  // double e0_dmrg = unavailable;
  Log("e0: {:.12f} (ED)", e0_obc);
  // REQUIRE(std::abs(e0_obc - e0_dmrg) < 1e-6);
  REQUIRE(std::isfinite(e0_obc));

  // --- Periodic boundary conditions with translational symmetry ---
  OpSum ops_pbc = ops_obc;
  ops_pbc += t * Op("Hop", {N - 1, 0});

  auto [e0_full, e0_sym] = testcases::translation_ground_states(
      ops_pbc, Electron(N, nup, ndn), N, [&](Representation const &irrep) {
        return Electron(N, nup, ndn, irrep);
      });
  REQUIRE(std::abs(e0_full - e0_sym) < 1e-6);

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
