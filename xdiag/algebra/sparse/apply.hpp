// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/common.hpp>

#include <xdiag/algebra/sparse/sparse_matrix_types.hpp>
#include <xdiag/extern/armadillo/armadillo>

namespace xdiag {

template <typename idx_t, typename coeff_t>
XDIAG_API arma::Col<coeff_t> apply(CSRMatrix<idx_t, coeff_t> const &spmat,
                                   arma::Col<coeff_t> const &vec_in);
template <typename idx_t, typename coeff_t>
XDIAG_API arma::Mat<coeff_t> apply(CSRMatrix<idx_t, coeff_t> const &spmat,
                                   arma::Mat<coeff_t> const &mat_in);

template <typename idx_t, typename coeff_t>
XDIAG_API void apply(CSRMatrix<idx_t, coeff_t> const &spmat,
                     arma::Col<coeff_t> const &vec_in,
                     arma::Col<coeff_t> &vec_out);
template <typename idx_t, typename coeff_t>
XDIAG_API void apply(CSRMatrix<idx_t, coeff_t> const &spmat,
                     arma::Mat<coeff_t> const &mat_in,
                     arma::Mat<coeff_t> &mat_out);

// For internal use only (needed for some generic algorithm, no QN checking)
template <typename idx_t, typename coeff_t, typename block_t>
inline void apply(CSRMatrix<idx_t, coeff_t> const &spmat, block_t const &,
                  arma::Col<coeff_t> const &vec_in, block_t const &,
                  arma::Col<coeff_t> &vec_out) {
  apply(spmat, vec_in, vec_out);
}

template <typename idx_t, typename coeff_t, typename block_t>
inline void apply(CSRMatrix<idx_t, coeff_t> const &spmat, block_t const &,
                  arma::Mat<coeff_t> const &mat_in, block_t const &,
                  arma::Mat<coeff_t> &mat_out) {
  apply(spmat, mat_in, mat_out);
}

template <typename idx_t, typename block_t>
inline void apply(CSRMatrix<idx_t, complex> const &spmat, block_t const &,
                  arma::vec const &vec_in, block_t const &,
                  arma::vec &vec_out) try {
  XDIAG_THROW("Multiplication of a complex CSR matrix with a real vector is "
              "not supported.");
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename idx_t, typename block_t>
inline void apply(CSRMatrix<idx_t, double> const &spmat, block_t const &,
                  arma::cx_vec const &vec_in, block_t const &,
                  arma::cx_vec &vec_out) try {
  XDIAG_THROW("Multiplication of a real CSR matrix with a complex vector is "
              "not supported.");
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename idx_t, typename block_t>
inline void apply(CSRMatrix<idx_t, complex> const &spmat, block_t const &,
                  arma::mat const &vec_in, block_t const &,
                  arma::mat &vec_out) try {
  XDIAG_THROW("Multiplication of a complex CSR matrix with a real matrix is "
              "not supported.");
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename idx_t, typename block_t>
inline void apply(CSRMatrix<idx_t, double> const &spmat, block_t const &,
                  arma::cx_mat const &vec_in, block_t const &,
                  arma::cx_mat &vec_out) try {
  XDIAG_THROW("Multiplication of a real CSR matrix with a complex matrix is "
              "not supported.");
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
