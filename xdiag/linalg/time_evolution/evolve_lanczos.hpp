// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <string>

#include <xdiag/armadillo.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/kernels/sparse/sparse_matrix_types.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/states/state.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

struct EvolveLanczosResult {
  XDIAG_API arma::vec alphas;
  XDIAG_API arma::vec betas;
  XDIAG_API arma::vec eigenvalues;
  XDIAG_API int64_t niterations;
  XDIAG_API std::string criterion;
  XDIAG_API State state;
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
  XDIAG_API arma::vec alphas;
  XDIAG_API arma::vec betas;
  XDIAG_API arma::vec eigenvalues;
  XDIAG_API int64_t niterations;
  XDIAG_API std::string criterion;
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
