// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/common.hpp>

#include <xdiag/algebra/sparse/sparse_matrix_types.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>

namespace xdiag {

struct EigvalsLanczosResult {
  arma::vec alphas;
  arma::vec betas;
  arma::vec eigenvalues;
  int64_t niterations;
  std::string criterion;
};

///////////////////////////////////////////////////////////////
// Routine with random state initialization

// on-the-fly
XDIAG_API EigvalsLanczosResult
eigvals_lanczos(OpSum const &ops, Block const &block, int64_t neigvals = 1,
                double precision = 1e-12, int64_t max_iterations = 1000,
                double deflation_tol = 1e-7, int64_t random_seed = 42);

// sparse matrix
template <typename idx_t, typename coeff_t>
XDIAG_API EigvalsLanczosResult
eigvals_lanczos(CSRMatrix<idx_t, coeff_t> const &A, Block const &block,
                int64_t neigvals = 1, double precision = 1e-12,
                int64_t max_iterations = 1000, double deflation_tol = 1e-7,
                int64_t random_seed = 42);

///////////////////////////////////////////////////////////////
// Routine with given starting state which is copied

// on-the-fly
XDIAG_API EigvalsLanczosResult eigvals_lanczos(OpSum const &ops, State psi0,
                                               int64_t neigvals = 1,
                                               double precision = 1e-12,
                                               int64_t max_iterations = 1000,
                                               double deflation_tol = 1e-7);
// sparse matrix
template <typename idx_t, typename coeff_t>
XDIAG_API EigvalsLanczosResult
eigvals_lanczos(CSRMatrix<idx_t, coeff_t> const &A, State psi0,
                int64_t neigvals = 1, double precision = 1e-12,
                int64_t max_iterations = 1000, double deflation_tol = 1e-7);

///////////////////////////////////////////////////////////////
// Routine with given starting state which is overwritten

// on-the-fly
XDIAG_API EigvalsLanczosResult
eigvals_lanczos_inplace(OpSum const &ops, State &psi0, int64_t neigvals = 1,
                        double precision = 1e-12, int64_t max_iterations = 1000,
                        double deflation_tol = 1e-7);

// sparse matrix
template <typename idx_t, typename coeff_t>
XDIAG_API EigvalsLanczosResult eigvals_lanczos_inplace(
    CSRMatrix<idx_t, coeff_t> const &A, State &psi0, int64_t neigvals = 1,
    double precision = 1e-12, int64_t max_iterations = 1000,
    double deflation_tol = 1e-7);

} // namespace xdiag
