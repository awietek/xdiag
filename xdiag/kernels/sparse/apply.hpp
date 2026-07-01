// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/armadillo.hpp>
#include <xdiag/kernels/sparse/sparse_matrix_types.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/xdiag_api.hpp>

// Sparse matrix-vector and matrix-matrix products for CSRMatrix.
//
// Two calling conventions are provided:
//   apply(A, x)         -> returns y = A*x  (allocates output)
//   apply(A, x, y)      -> computes y = A*x  (writes into pre-allocated y)
//
// Mixed-precision overloads allow a real CSRMatrix to be applied to a complex
// vector/matrix (the matrix entries are implicitly promoted).
//
// When compiled with XDIAG_USE_SPARSE_MKL and idx_t == int64_t the MKL
// inspector-executor SpMV/SpMM routines are used for the real/complex cases.
//
// The four-argument overload (with two block arguments) is for internal use by
// generic algorithms that pass blocks alongside the vectors; the blocks are
// ignored and no quantum-number checking is performed.

namespace xdiag {

template <typename idx_t, typename coeff_t>
XDIAG_API arma::Col<coeff_t> apply(CSRMatrix<idx_t, coeff_t> const &spmat,
                                   arma::Col<coeff_t> const &vec_in);
template <typename idx_t>
XDIAG_API arma::Col<complex> apply(CSRMatrix<idx_t, double> const &spmat,
                                   arma::Col<complex> const &vec_in);
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
                     arma::Col<complex> const &vec_in,
                     arma::Col<complex> &vec_out);
template <typename idx_t>
XDIAG_API void apply(CSRMatrix<idx_t, double> const &spmat,
                     arma::Mat<complex> const &mat_in,
                     arma::Mat<complex> &mat_out);

// Internal: four-argument form used by generic algorithms (blocks ignored).
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
