// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cmath>
#include <functional>
#include <limits>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag::testcases {

// Ground-state energy of a translation-invariant (PBC) Hamiltonian computed two
// ways: (1) on the full block without symmetry, and (2) as the minimum over all
// translation-momentum sectors of the cyclic group on `nsites` sites. For a
// genuinely translation-invariant Hamiltonian the two must agree; this
// exercises the symmetric BitsetDynamic kernels for the "long" (>64 site)
// blocks. `sector(irrep)` returns the symmetry-adapted block for a given irrep.
//
// Returns {e0_full, e0_sym} so the caller can REQUIRE their agreement.
inline std::pair<double, double> translation_ground_states(
    OpSum const &ops_pbc, Block const &full, int64_t nsites,
    std::function<Block(Representation const &)> const &sector) {
  double e0_full = eigval0(ops_pbc, full);

  double e0_sym = std::numeric_limits<double>::max();
  for (int64_t k = 0; k < nsites; ++k) {
    Representation irrep = cyclic_group_irrep(nsites, k);
    Block block_k = sector(irrep);
    if (dim(block_k) > 0) {
      e0_sym = std::min(e0_sym, eigval0(ops_pbc, block_k));
    }
  }
  return {e0_full, e0_sym};
}

} // namespace xdiag::testcases
