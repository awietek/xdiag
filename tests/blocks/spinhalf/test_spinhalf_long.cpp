// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0
//
// "Long" Hilbert space test for the Spinhalf block (N=65 -> BitsetDynamic).
// Heisenberg chain at a low, fixed magnetization sector.
//
//   H = J * sum_<ij> S_i . S_j
//
// N=65 is the smallest chain that activates the BitsetDynamic bit type (>64
// sites); with nup=2 the sector dimension is C(65,2)=2080.


#include <cmath>

#include <tests/catch.hpp>
#include <tests/blocks/test_long_blocks.hpp>

#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

TEST_CASE("spinhalf_long", "[long]") try {
  Log("Spinhalf long test");
  
  int64_t N = 65;
  double J = 1.0;
  int64_t nup = 2;

  // --- Open boundary conditions ---
  OpSum ops_obc;
  for (int64_t i = 0; i < N - 1; ++i) {
    ops_obc += J * Op("SdotS", {i, i + 1});
  }
  double e0_obc = eigval0(ops_obc, Spinhalf(N, nup));
  double e0_dmrg = 12.0059269353772;
  Log("e0: {:.12f} (ED), {:.12f} (DMRG)", e0_obc, e0_dmrg);
  REQUIRE(std::abs(e0_obc - e0_dmrg) < 1e-6);
  REQUIRE(std::isfinite(e0_obc));

  // --- Periodic boundary conditions with translational symmetry ---
  OpSum ops_pbc = ops_obc;
  ops_pbc += J * Op("SdotS", {N - 1, 0});

  auto [e0_full, e0_sym] = testcases::translation_ground_states(
      ops_pbc, Spinhalf(N, nup), N,
      [&](Representation const &irrep) { return Spinhalf(N, nup, irrep); });
  REQUIRE(std::abs(e0_full - e0_sym) < 1e-6);

} catch (xdiag::Error const &e) {
  error_trace(e);
  throw;
}
