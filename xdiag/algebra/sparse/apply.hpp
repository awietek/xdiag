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
template <typename idx_t>
XDIAG_API arma::Col<complex> apply(CSRMatrix<idx_t, double> const &spmat,
                                   arma::Col<complex> const &mat_in);
template <typename idx_t, typename coeff_t>
XDIAG_API arma::Mat<coeff_t> apply(CSRMatrix<idx_t, coeff_t> const &spmat,
                                   arma::Mat<coeff_t> const &mat_in);
template <typename idx_t>
XDIAG_API arma::Mat<complex> apply(CSRMatrix<idx_t, double> const &spmat,
                                   arma::Mat<complex> const &mat_in);

template <typename idx_t, typename coeff_t>
XDIAG_API void apply(CSRMatrix<idx_t, coeff_t> const &spmat,
                     arma::Col<coeff_t> const &vec_in,
                     arma::Col<coeff_t> &vec_out);
template <typename idx_t, typename coeff_t>
XDIAG_API void apply(CSRMatrix<idx_t, coeff_t> const &spmat,
                     arma::Mat<coeff_t> const &mat_in,
                     arma::Mat<coeff_t> &mat_out);
template <typename idx_t>
XDIAG_API void apply(CSRMatrix<idx_t, double> const &spmat,
                     arma::Col<complex> const &mat_in,
                     arma::Col<complex> &mat_out);
template <typename idx_t>
XDIAG_API void apply(CSRMatrix<idx_t, double> const &spmat,
                     arma::Mat<complex> const &mat_in,
                     arma::Mat<complex> &mat_out);

// For internal use only (needed for some generic algorithm, no QN checking)
template <typename idx_t, typename coeff_t, typename block_t, typename vec_t>
inline void apply(CSRMatrix<idx_t, coeff_t> const &spmat, block_t const &,
                  vec_t const &vec_in, block_t const &, vec_t &vec_out) try {
  if constexpr (isreal<typename vec_t::elem_type>() && !isreal<coeff_t>()) {
    XDIAG_THROW("Cannot apply a complex CSRMatrix to a real vector.");
  } else {
    return apply(spmat, vec_in, vec_out);
  }
}
XDIAG_CATCH

} // namespace xdiag
