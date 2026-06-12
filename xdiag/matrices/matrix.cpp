// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "matrix.hpp"

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/matrices/blocks/boson/dispatch_basis.hpp>
#include <xdiag/matrices/blocks/fermion/dispatch_basis.hpp>
#include <xdiag/matrices/blocks/spinhalf/dispatch_basis.hpp>
#include <xdiag/matrices/kernels.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/variants.hpp>

// Dispatch overview
// -----------------
// The matrix function routes through two layers before reaching the numerical
// kernel, mirroring apply.cpp:
//
// Layer 1 — Block type (variant dispatch via visit_same_type):
//   matrix(op_t, Block)  →  matrix(op_t, Block, Block)
//     The output Block is first determined from the operator quantum numbers,
//     then visit_same_type unwraps both Block variants to their concrete type
//     and calls matrix_impl(OpSum, ConcreteBlock, ConcreteBlock, coeff_t*).
//     Op/Monomial are promoted to OpSum here. The output arma matrix is
//     allocated and zeroed at this layer; a raw pointer is passed down so the
//     kernel can fill it without owning the storage.
//
// Layer 2 — matrix_impl<block_t>: one body for every block type. The block's
//     dispatch_basis overload (matrices/blocks/<block>/dispatch_basis.hpp)
//     resolves the shared_ptr<Basis> to a concrete BasisOnTheFly<...> via a
//     type-id table. The numerical kernel is selected by block_t in
//     kernels.cpp. A block type with no dispatch_basis overload is a compile
//     error here, never a silent fallback.
//
// Kernel — matrices::matrix<block_t, coeff_t>(ops, basis_in, basis_out, mat)
//     Defined in kernels.cpp; only declared here. Each (block_t, basis_t)
//     specialisation is compiled in its own translation unit (instantiation
//     groups), keeping this file free of heavy template instantiation.

namespace xdiag {

// Layer 2: block-generic dense fill. One body for all block types via the
// dispatch_basis overload set.
template <typename block_t, typename coeff_t>
static void matrix_impl(OpSum const &ops, block_t const &block_in,
                        block_t const &block_out, coeff_t *mat) try {
  matrices::dispatch_basis(
      block_in, block_out, [&](auto const &basis_in, auto const &basis_out) {
        // Kernel: definition is in kernels.cpp, instantiated per basis type.
        matrices::matrix<block_t>(ops, basis_in, basis_out, mat);
      });
}
XDIAG_CATCH

// Developer overload used by the Julia wrapper: thin wrapper over matrix_impl
// that pins the concrete block type (keeps a stable exported symbol).
template <typename coeff_t>
void matrix(OpSum const &ops, Spinhalf const &block_in,
            Spinhalf const &block_out, coeff_t *mat) try {
  matrix_impl(ops, block_in, block_out, mat);
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

  int64_t m = size(block_out);
  int64_t n = dim(block_in);
  arma::mat mat(m, n, arma::fill::zeros); // zeroed here; kernel uses +=
  // Layer 1: unwrap the Block variant (op_t is promoted to OpSum inside) and
  // forward to the block-generic matrix_impl.
  utils::visit_same_type(
      block_in, block_out,
      [&](auto const &bin, auto const &bout) {
        matrix_impl(OpSum(op), bin, bout, mat.memptr());
      },
      "Type mismatch of Block types");
  return mat;
}
XDIAG_CATCH

template <typename op_t>
arma::cx_mat matrixC(op_t const &op, Block const &block_in,
                     Block const &block_out) try {
  int64_t m = size(block_out);
  int64_t n = dim(block_in);
  arma::cx_mat mat(m, n, arma::fill::zeros); // zeroed here; kernel uses +=
  // Layer 1: unwrap the Block variant (op_t is promoted to OpSum inside) and
  // forward to the block-generic matrix_impl.
  utils::visit_same_type(
      block_in, block_out,
      [&](auto const &bin, auto const &bout) {
        matrix_impl(OpSum(op), bin, bout, mat.memptr());
      },
      "Type mismatch of Block types");
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

// Exported developer overload (Julia wrapper). The Boson / future-block paths
// go through matrix_impl, instantiated implicitly via the Layer-1 templates
// above, so they need no explicit instantiation here.
template void matrix(OpSum const &, Spinhalf const &, Spinhalf const &,
                     double *);
template void matrix(OpSum const &, Spinhalf const &, Spinhalf const &,
                     complex *);

} // namespace xdiag
