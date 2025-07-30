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

} // namespace xdiag
