// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <mpi.h>

#include <tests/blocks/spinhalf/testcases_spinhalf.hpp>
#include <tests/catch.hpp>

#include <xdiag/blocks/distributed/spinhalf_distributed.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/kernels/apply.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/linalg/lobpcg/eigs_lobpcg.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/utils/logger.hpp>

using namespace xdiag;

// Checks that LOBPCG works on a distributed block, i.e. that the block
// matrix-vector multiply apply(ops, block, X, block, Y) and the distributed
// Gram products (matrix_dot -> cdot_distributed) work for multi-column blocks
// under MPI. The Hamiltonian is complex Hermitian: a Heisenberg chain plus a
// Dzyaloshinskii-Moriya-like term i*d*ExchangeAsym (ExchangeAsym is
// anti-hermitian, so an imaginary prefactor makes it hermitian). Eigenvalues
// are compared against a serial dense reference.
TEST_CASE("lobpcg_distributed", "[linalg]") {
  using namespace xdiag::testcases::spinhalf;

  int nsites = 8;
  int nup = 4;
  OpSum ops = HBchain(nsites, 1.0);      // periodic Heisenberg chain (real)
  ops["Dz"] = complex(0.0, 0.4);         // imaginary -> hermitian DM term
  for (int s = 0; s < nsites; ++s) {
    ops += "Dz" * Op("ExchangeAsym", {s, (s + 1) % nsites});
  }

  // Serial dense reference spectrum (computed redundantly on every rank).
  auto block = Spinhalf(nsites, nup);
  arma::cx_mat H = matrixC(ops, block);
  REQUIRE(H.is_hermitian(1e-8)); // sanity: the DM term keeps H hermitian
  arma::vec evals_ref;
  arma::eig_sym(evals_ref, H);

  // Distributed block: exercises the distributed complex block-matrix multiply.
  auto block_mpi = SpinhalfDistributed(nsites, nup);
  REQUIRE_FALSE(isreal(ops)); // ensure we exercise the complex path

  int64_t neigs = 3;
  int64_t guard = 3;
  REQUIRE(5 * (neigs + guard) <= dim(block_mpi));

  EigsLobpcgResult r = eigs_lobpcg(ops, block_mpi, neigs, guard, 1e-8, 1000);

  for (int64_t i = 0; i < neigs; ++i) {
    REQUIRE(std::abs(r.eigenvalues(i) - evals_ref(i)) < 1e-6);
  }
  REQUIRE(r.residual_norms.max() < 1e-4);
}
