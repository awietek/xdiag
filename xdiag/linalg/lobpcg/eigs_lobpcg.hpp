// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/kernels/sparse/sparse_matrix_types.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

struct EigsLobpcgResult {
  arma::vec eigenvalues;    // the neigs smallest eigenvalues, ascending
  arma::vec residual_norms; // residual norm per returned eigenvalue
  State eigenvectors;       // neigs eigenvectors (columns)
  int64_t niterations;
  std::string criterion;
  // Per-iteration Ritz-value / residual-norm history (one row per iteration,
  // blocksize = neigs + guard columns) for convergence diagnostics.
  arma::mat eigenvalue_history;
  arma::mat residual_norms_history;
};

// Computes the "neigs" algebraically smallest eigenvalues (and eigenvectors) of
// a hermitian operator with the LOBPCG block eigensolver. A block of size
// "neigs + guard" is iterated so degenerate multiplets at the neigs-th
// eigenvalue are captured with the correct multiplicity. The initial block is
// filled with (seeded) random vectors.

// on-the-fly
XDIAG_API EigsLobpcgResult
eigs_lobpcg(OpSum const &ops, Block const &block, int64_t neigs = 1,
            int64_t guard = 2, double tol = 1e-10,
            int64_t max_iterations = 1000, int64_t random_seed = 42);

// sparse matrix
template <typename idx_t, typename coeff_t>
XDIAG_API EigsLobpcgResult
eigs_lobpcg(CSRMatrix<idx_t, coeff_t> const &A, Block const &block,
            int64_t neigs = 1, int64_t guard = 2, double tol = 1e-10,
            int64_t max_iterations = 1000, int64_t random_seed = 42);

} // namespace xdiag
