// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0
//
// "Long" Hilbert space test for the Boson block (N=65 -> BitsetDynamic).
// Bose-Hubbard chain, local dimension d = 4, low fixed particle number.
//
//   H = -t sum_i (b^dag_i b_{i+1} + h.c.) + (U/2) sum_i n_i (n_i - 1)
//
// The Boson block expands Op("Hop") -> -(Adag A + h.c.) and Op("HubbardU") ->
// (1/2) sum_i N_i (N_i - 1) via the matrix algebra, so the couplings passed are
// +t and U. Reference from dmrg_boson.jl: OBC, N=128, t=U=1, d=4, Nb=3.

#include <cmath>

#include <tests/catch.hpp>
#include <tests/blocks/test_long_blocks.hpp>

#include <xdiag/blocks/boson.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

TEST_CASE("boson_long", "[long]") try {
  Log("Boson long test");

  int64_t N = 65; // smallest chain activating BitsetDynamic (>64 sites)
  int64_t d = 4;  // local boson dimension (maxOcc)
  double t = 1.0;
  double U = 1.0;
  int64_t Nb = 2; // dim ~ 2145

  // On-site interaction (same for OBC and PBC: it lives on single sites).
  OpSum ops_int = U * Op("HubbardU"); // (U/2) sum_i N_i (N_i - 1)

  // --- Open boundary conditions ---
  OpSum ops_obc = ops_int;
  for (int64_t i = 0; i < N - 1; ++i) {
    ops_obc += t * Op("Hop", {i, i + 1}); // Hop = -(b^dag b + h.c.)
  }
  double e0_obc = eigval0(ops_obc, Boson(N, d, Nb));
  double e0_dmrg = -3.98990907896486;
  Log("e0: {:.12f} (ED), {:.12f} (DMRG)", e0_obc, e0_dmrg);
  REQUIRE(std::abs(e0_obc - e0_dmrg) < 1e-6);
  REQUIRE(std::isfinite(e0_obc));

  // --- Periodic boundary conditions with translational symmetry ---
  OpSum ops_pbc = ops_obc;
  ops_pbc += t * Op("Hop", {N - 1, 0});

  auto [e0_full, e0_sym] = testcases::translation_ground_states(
      ops_pbc, Boson(N, d, Nb), N,
      [&](Representation const &irrep) { return Boson(N, d, Nb, irrep); });
  REQUIRE(std::abs(e0_full - e0_sym) < 1e-6);

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
