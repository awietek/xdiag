// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "matrix.hpp"

#include <algorithm>

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/kernels/compile.hpp>
#include <xdiag/kernels/elements.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

// Fills a caller-allocated, column-major buffer of size size(block_out) x
// dim(block_in): the developer overload, also used by the Julia wrapper and by
// the arma-returning overloads below. The pointer is not retained. Every row is
// written by one thread, so no synchronisation is needed.
template <typename coeff_t>
void matrix(OpSum const &ops, Block const &block_in, Block const &block_out,
            coeff_t *mat) try {
  auto ops_compiled = kernels::compile(ops, block_in, block_out);
  int64_t m = size(block_out);
  int64_t n = dim(block_in);
  std::fill_n(mat, m * n, coeff_t(0));
  kernels::elements<coeff_t>(ops_compiled,
                             [&](int64_t row, int64_t col, coeff_t val) {
                               mat[row + col * m] += val;
                             });
}
XDIAG_CATCH

template <typename op_t>
arma::mat matrix(op_t const &op, Block const &blocki) try {
  auto blockr = block(OpSum(op), blocki);
  return matrix(op, blocki, blockr);
}
XDIAG_CATCH

template <typename op_t>
arma::cx_mat matrixC(op_t const &op, Block const &blocki) try {
  auto blockr = block(OpSum(op), blocki);
  return matrixC(op, blocki, blockr);
}
XDIAG_CATCH

template <typename op_t>
arma::mat matrix(op_t const &op, Block const &block_in,
                 Block const &block_out) try {
  if (!isreal(op)) {
    XDIAG_THROW("Cannot create a real matrix when the input operator is "
                "complex. Consider using matrixC instead.");
  }
  if (!isreal(block_in)) {
    XDIAG_THROW("Cannot create a real matrix when the input block is complex. "
                "Consider using matrixC instead.");
  }
  if (!isreal(block_out)) {
    XDIAG_THROW("Cannot create a real matrix when the output block is complex. "
                "Consider using matrixC instead.");
  }
  arma::mat mat(size(block_out), dim(block_in));
  matrix(OpSum(op), block_in, block_out, mat.memptr());
  return mat;
}
XDIAG_CATCH

template <typename op_t>
arma::cx_mat matrixC(op_t const &op, Block const &block_in,
                     Block const &block_out) try {
  arma::cx_mat mat(size(block_out), dim(block_in));
  matrix(OpSum(op), block_in, block_out, mat.memptr());
  return mat;
}
XDIAG_CATCH

#define INSTANTIATE_XDIAG_MATRIX(OP_TYPE)                                      \
  template arma::mat matrix(OP_TYPE const &, Block const &);                   \
  template arma::cx_mat matrixC(OP_TYPE const &, Block const &);               \
  template arma::mat matrix(OP_TYPE const &op, Block const &, Block const &);  \
  template arma::cx_mat matrixC(OP_TYPE const &, Block const &, Block const &);

INSTANTIATE_XDIAG_MATRIX(Op);
INSTANTIATE_XDIAG_MATRIX(Monomial);
INSTANTIATE_XDIAG_MATRIX(OpSum);

#undef INSTANTIATE_XDIAG_MATRIX

// Exported developer overload (Julia wrapper), one per coeff type; the blocks
// are handled inside kernels::elements.
template void matrix(OpSum const &, Block const &, Block const &, double *);
template void matrix(OpSum const &, Block const &, Block const &, complex *);

} // namespace xdiag
