// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <tests/catch.hpp>

#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/kernels/sparse/coo_matrix.hpp>
#include <xdiag/kernels/sparse/csc_matrix.hpp>
#include <xdiag/kernels/sparse/csr_matrix.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/apply.hpp>
#include <xdiag/states/fill.hpp>
#include <xdiag/states/random_state.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/utils/error.hpp>

using namespace xdiag;

// Regression test for github.com/awietek/xdiag/discussions/110: an OpSum which
// does not commute with the symmetry group of the block it is diagonalized on
// used to be applied silently, since the iterative solvers call the kernel-level
// apply directly and thereby bypassed the blocks_match check that the
// State-level apply performs. The resulting energies were meaningless (and not
// even variational -- in the discussion they fell below the true ground state).
TEST_CASE("blocks_match_symmetry", "[blocks]") {
  int64_t nsites = 12;

  // Heisenberg chain with periodic boundaries: invariant under translations.
  OpSum ops;
  for (int64_t i = 0; i < nsites; ++i) {
    ops += Op("SdotS", {i, (i + 1) % nsites});
  }

  // Same, but with a single bond moved off the symmetry orbit. This is exactly
  // the shape of the typo in the discussion: one bond of the lattice file
  // connected the wrong pair of sites, so no translation maps the bond set onto
  // itself any more.
  OpSum ops_broken;
  for (int64_t i = 0; i < nsites; ++i) {
    if (i == 3) {
      ops_broken += Op("SdotS", {i, (i + 2) % nsites});
    } else {
      ops_broken += Op("SdotS", {i, (i + 1) % nsites});
    }
  }

  auto irrep = cyclic_group_irrep(nsites, 0);
  auto block_sym = Spinhalf(nsites, irrep);
  auto block_full = Spinhalf(nsites);

  // The symmetric OpSum still works and agrees with the full Hilbert space.
  {
    double e0_sym = eigval0(ops, block_sym);
    double e0_full = eigval0(ops, block_full);
    REQUIRE(std::abs(e0_sym - e0_full) < 1e-8);
  }

  // The non-invariant OpSum must throw rather than return a number, on every
  // route into the kernels: the iterative solvers, the dense matrix builders,
  // and apply.
  {
    REQUIRE_THROWS_AS(eigval0(ops_broken, block_sym), Error);
    REQUIRE_THROWS_AS(eig0(ops_broken, block_sym), Error);
    REQUIRE_THROWS_AS(eigvals(ops_broken, block_sym, 1), Error);
    REQUIRE_THROWS_AS(matrixC(ops_broken, block_sym, block_sym), Error);

    auto v = State(block_sym);
    fill(v, RandomState(42));
    auto w = State(block_sym);
    REQUIRE_THROWS_AS(apply(ops_broken, v, w), Error);
  }

  // The sparse route is checked when the matrix is assembled, so a CSR/CSC/COO
  // matrix can never be built from a non-invariant OpSum on a symmetric block
  // in the first place. Applying an already-assembled sparse matrix therefore
  // needs no further check.
  {
    REQUIRE(csr_matrixC(ops, block_sym, block_sym).data.n_elem > 0);
    REQUIRE_THROWS_AS(csr_matrixC(ops_broken, block_sym, block_sym), Error);
    REQUIRE_THROWS_AS(csc_matrixC(ops_broken, block_sym, block_sym), Error);
    REQUIRE_THROWS_AS(coo_matrixC(ops_broken, block_sym, block_sym), Error);
    REQUIRE_THROWS_AS(csr_matrixC(ops_broken, block_sym), Error);
  }

  // Without symmetries the same OpSum is perfectly legitimate: nothing is being
  // assumed about it, so no check may fire.
  {
    double e0 = eigval0(ops_broken, block_full);
    REQUIRE(std::isfinite(e0));
  }

  // The check must also catch the "Matrix" Op spelling used in the discussion,
  // where the two-site coupling is given as an explicit 4x4 matrix.
  {
    arma::cx_mat sx = {{complex(0, 0), complex(0.5, 0)},
                       {complex(0.5, 0), complex(0, 0)}};
    arma::cx_mat sy = {{complex(0, 0), complex(0, -0.5)},
                       {complex(0, 0.5), complex(0, 0)}};
    arma::cx_mat sz = {{complex(0.5, 0), complex(0, 0)},
                       {complex(0, 0), complex(-0.5, 0)}};
    arma::cx_mat pair =
        arma::kron(sx, sx) + arma::kron(sy, sy) + arma::kron(sz, sz);

    OpSum ops_mat, ops_mat_broken;
    for (int64_t i = 0; i < nsites; ++i) {
      ops_mat += Op("Matrix", {i, (i + 1) % nsites}, pair);
      if (i == 3) {
        ops_mat_broken += Op("Matrix", {i, (i + 2) % nsites}, pair);
      } else {
        ops_mat_broken += Op("Matrix", {i, (i + 1) % nsites}, pair);
      }
    }

    REQUIRE(std::isfinite(eigval0(ops_mat, block_sym)));
    REQUIRE_THROWS_AS(eigval0(ops_mat_broken, block_sym), Error);
  }
}
