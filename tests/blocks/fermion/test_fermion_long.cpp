// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0
//
// "Long" Hilbert space test for the Fermion block (N=65 -> BitsetDynamic).
// Spinless fermions, nearest-neighbour hopping + interaction, low fixed filling.
//
//   H = -t sum_i (c^dag_i c_{i+1} + h.c.) + V sum_i n_i n_{i+1}
//
// xdiag Op("Hop") already carries the minus sign (Hop = -(c^dag c + h.c.)), so
// the coupling passed is +t. N=65 is the smallest chain activating the
// BitsetDynamic bit type (>64 sites); with Nf=2 the dimension is C(65,2)=2080.

#include <cmath>

#include <tests/catch.hpp>
#include <tests/blocks/test_long_blocks.hpp>

#include <xdiag/blocks/fermion.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

TEST_CASE("fermion_long", "[long]") try {
  Log("Fermion long test");

  int64_t N = 65; // smallest chain activating BitsetDynamic (>64 sites)
  double t = 1.0;
  double V = 1.0;
  int64_t Nf = 2; // C(65,2) = 2080

  // --- Open boundary conditions ---
  OpSum ops_obc;
  for (int64_t i = 0; i < N - 1; ++i) {
    ops_obc += t * Op("Hop", {i, i + 1}); // Hop = -(c^dag c + h.c.)
    ops_obc += V * Op("NN", {i, i + 1});  // NN  = n_i n_j
  }
  double e0_obc = eigval0(ops_obc, Fermion(N, Nf));
  double e0_dmrg = -3.9885633480006;
  Log("e0: {:.12f} (ED), {:.12f} (DMRG)", e0_obc, e0_dmrg);
  REQUIRE(std::abs(e0_obc - e0_dmrg) < 1e-6);
  REQUIRE(std::isfinite(e0_obc));

  // --- Periodic boundary conditions with translational symmetry ---
  OpSum ops_pbc = ops_obc;
  ops_pbc += t * Op("Hop", {N - 1, 0}); // wrap-around bond
  ops_pbc += V * Op("NN", {N - 1, 0});

  auto [e0_full, e0_sym] = testcases::translation_ground_states(
      ops_pbc, Fermion(N, Nf), N,
      [&](Representation const &irrep) { return Fermion(N, Nf, irrep); });
  REQUIRE(std::abs(e0_full - e0_sym) < 1e-6);

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
