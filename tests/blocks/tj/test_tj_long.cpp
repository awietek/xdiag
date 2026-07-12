// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0
//
// "Long" Hilbert space test for the tJ block (>64 sites -> BitsetDynamic).
// Nearest-neighbour t-J chain at low fixed particle numbers.
//
//   H = -t sum_<ij>,sigma (c^dag c + h.c.) + J sum_<ij> (S_i.S_j - n_i n_j / 4)
//
// xdiag Op("Hop") carries the minus sign, Op("tJSdotS") = S.S - n_i n_j/4, so
// the couplings passed are +t and J. A converged DMRG reference is not yet
// available for this model, so the open-boundary energy is still computed but
// the comparison against the DMRG value is commented out. The periodic-boundary
// translational-symmetry consistency check is fully active.

#include <cmath>

#include <tests/blocks/test_long_blocks.hpp>
#include <tests/catch.hpp>

#include <xdiag/blocks/tj.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/timing.hpp>

using namespace xdiag;

TEST_CASE("tj_long", "[long]") try {
  Log("tJ long test");

  int64_t N = 65; // smallest chain activating BitsetDynamic (>64 sites)
  double t = 1.0;
  double J = 1.0;
  int64_t nup = 1;
  int64_t ndn = 1; // dim = 65 * 64 = 4160

  // --- Open boundary conditions ---
  OpSum ops_obc;
  for (int64_t i = 0; i < N - 1; ++i) {
    ops_obc += t * Op("Hop", {i, i + 1}); // Hop = -sum_sigma(c^dag c + h.c.)
    ops_obc += J * Op("tJSdotS", {i, i + 1}); // S.S - n_i n_j / 4
  }
  double e0_obc = eigval0(ops_obc, tJ(N, nup, ndn));

  Log("e0: {:.6f} (ED)", e0_obc);

  // No converged DMRG reference yet; enable once available:
  // double e0_dmrg = /* from dmrg_tj.jl */;
  // REQUIRE(std::abs(e0_obc - e0_dmrg) < 1e-6);
  REQUIRE(std::isfinite(e0_obc));

  // --- Periodic boundary conditions with translational symmetry ---
  std::vector<int64_t> Ns = {32, 65};
  for (int64_t N : Ns) {
    Log("N={}", N);
    tic();
    OpSum ops_pbc;
    for (int64_t i = 0; i < N; ++i) {
      ops_pbc += t * Op("Hop", {i, (i + 1) % N});
      ops_pbc += J * Op("tJSdotS", {i, (i + 1) % N});
    }

    auto [e0_full, e0_sym] = testcases::translation_ground_states(
        ops_pbc, tJ(N, nup, ndn), N,
        [&](Representation const &irrep) { return tJ(N, nup, ndn, irrep); });
    REQUIRE(std::abs(e0_full - e0_sym) < 1e-6);
    toc();
  }

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
