// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/algebra/sparse/sparse_matrix_types.hpp>
#include <xdiag/algorithms/lanczos/lanczos.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/timing.hpp>

namespace xdiag {

struct EvolveLanczosResult {
  arma::vec alphas;
  arma::vec betas;
  arma::vec eigenvalues;
  int64_t niterations;
  std::string criterion;
  State state;
};

XDIAG_API EvolveLanczosResult
evolve_lanczos(OpSum const &H, State psi, double tau, double precision = 1e-12,
               double shift = 0., bool normalize = false,
               int64_t max_iterations = 1000, double deflation_tol = 1e-7);

XDIAG_API EvolveLanczosResult
evolve_lanczos(OpSum const &H, State psi, complex tau, double precision = 1e-12,
               double shift = 0., bool normalize = false,
               int64_t max_iterations = 1000, double deflation_tol = 1e-7);

template <typename idx_t, typename coeff_t>
XDIAG_API EvolveLanczosResult evolve_lanczos(
    CSRMatrix<idx_t, coeff_t> const &H, State psi, double tau,
    double precision = 1e-12, double shift = 0., bool normalize = false,
    int64_t max_iterations = 1000, double deflation_tol = 1e-7);

template <typename idx_t, typename coeff_t>
XDIAG_API EvolveLanczosResult evolve_lanczos(
    CSRMatrix<idx_t, coeff_t> const &H, State psi, complex tau,
    double precision = 1e-12, double shift = 0., bool normalize = false,
    int64_t max_iterations = 1000, double deflation_tol = 1e-7);

struct EvolveLanczosInplaceResult {
  arma::vec alphas;
  arma::vec betas;
  arma::vec eigenvalues;
  int64_t niterations;
  std::string criterion;
};

XDIAG_API EvolveLanczosInplaceResult evolve_lanczos_inplace(
    OpSum const &H, State &psi, double tau, double precision = 1e-12,
    double shift = 0., bool normalize = false, int64_t max_iterations = 1000,
    double deflation_tol = 1e-7);

XDIAG_API EvolveLanczosInplaceResult evolve_lanczos_inplace(
    OpSum const &H, State &psi, complex tau, double precision = 1e-12,
    double shift = 0., bool normalize = false, int64_t max_iterations = 1000,
    double deflation_tol = 1e-7);

template <typename idx_t, typename coeff_t>
XDIAG_API EvolveLanczosInplaceResult evolve_lanczos_inplace(
    CSRMatrix<idx_t, coeff_t> const &H, State &psi, double tau,
    double precision = 1e-12, double shift = 0., bool normalize = false,
    int64_t max_iterations = 1000, double deflation_tol = 1e-7);

template <typename idx_t, typename coeff_t>
XDIAG_API EvolveLanczosInplaceResult evolve_lanczos_inplace(
    CSRMatrix<idx_t, coeff_t> const &H, State &psi, complex tau,
    double precision = 1e-12, double shift = 0., bool normalize = false,
    int64_t max_iterations = 1000, double deflation_tol = 1e-7);

} // namespace xdiag
