// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/blocks.hpp>
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

namespace kernels {

// Internal. Counterparts of kernels::compile and kernels::apply(CompiledOpSum)
// for a CSRMatrix, so generic algorithms handle either operator type. A sparse
// matrix needs no compilation, and the row-parallel CSR product writes each
// output element from exactly one thread, so it is bitwise reproducible as is.
template <typename idx_t, typename coeff_t>
inline CSRMatrix<idx_t, coeff_t> const &
compile(CSRMatrix<idx_t, coeff_t> const &spmat, Block const &, Block const &) {
  return spmat;
}

template <typename idx_t, typename coeff_t, typename vec_t>
inline void apply(CSRMatrix<idx_t, coeff_t> const &spmat, vec_t const &vec_in,
                  vec_t &vec_out) try {
  if constexpr (isreal<typename vec_t::elem_type>() && !isreal<coeff_t>()) {
    XDIAG_THROW("Cannot apply a complex CSRMatrix to a real vector.");
  } else {
    xdiag::apply(spmat, vec_in, vec_out);
  }
}
XDIAG_CATCH

} // namespace kernels
} // namespace xdiag
