// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/blocks/electron/testcases_electron.hpp>
#include <tests/catch.hpp>

#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/linalg/lobpcg/eigs_lobpcg.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/math/dot.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/xdiag_show.hpp>

#include <cmath>
#include <vector>

using namespace xdiag;

TEST_CASE("lobpcg_distributed", "[linalg]") {
  using namespace xdiag::testcases::electron;

  int nsites = 6;
  int nup = 3;
  int ndn = 3;
  OpSum ops = freefermion_alltoall(nsites);
  auto block = Electron(nsites, nup, ndn);
  auto blockd = ElectronDistributed(nsites, nup, ndn);

  int64_t neigs = 4;

  // High-level API: lowest neigs eigenvalues with guard vectors.
  Log("lobpcg: Distributed");
  EigsLobpcgResult r = eigs_lobpcg(ops, block, neigs, /*guard=*/3, /*tol=*/1e-9,
                                   /*max_iterations=*/500);
  double e0 = eigval0(ops, block);
  REQUIRE(isapprox(r.eigenvalues(0), e0));
}
