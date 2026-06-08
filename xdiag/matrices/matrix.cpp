// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "matrix.hpp"

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/matrices/boson/dispatch_basis.hpp>
#include <xdiag/matrices/boson/matrix_policy.hpp>
#include <xdiag/matrices/kernels.hpp>
#include <xdiag/matrices/spinhalf/dispatch_basis.hpp>
#include <xdiag/matrices/spinhalf/matrix_policy.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/variants.hpp>

// Dispatch overview
// -----------------
// The matrix function routes through two layers of type erasure before
// reaching the actual numerical kernel, mirroring apply.cpp exactly:
//
// Layer 1 — Block type (variant dispatch via visit_same_type):
//   matrix(op_t, Block)  →  matrix(op_t, Block, Block)
//     The output Block is first determined from the operator quantum numbers,
//     then visit_same_type unwraps both Block variants to their concrete type
//     and calls matrix(OpSum, ConcreteBlock, ConcreteBlock, coeff_t*).
//     Op/Monomial are promoted to OpSum here. The output arma matrix is
//     allocated and zeroed at this layer; a raw pointer is passed down so the
//     kernel can fill it without owning the storage.
//
// Layer 2 — Basis type (runtime dispatch via dispatch_basis):
//   matrix(OpSum, Spinhalf, Spinhalf, coeff_t*)   [public: used by Julia
//   wrapper]
//     dispatch_basis resolves the shared_ptr<Basis> stored in each Spinhalf
//     to a concrete BasisOnTheFly<...> type via a type-id lookup table, then
//     calls the lambda with the concrete basis objects.
//
// Kernel — spinhalf::matrix<basis_t, coeff_t>(ops, basis_in, basis_out, mat)
//     Defined in spinhalf/kernels.cpp; only declared here. Each basis_t
//     specialisation is compiled in its own translation unit (instantiation
//     groups), keeping this file free of heavy template instantiation.

namespace xdiag {

template <typename coeff_t>
void matrix(OpSum const &ops, Spinhalf const &block_in,
            Spinhalf const &block_out, coeff_t *mat) try {
  // Layer 2: unwrap the basis pointer; the lambda receives the concrete
  // BasisOnTheFly<...> type, allowing matrix_generic to be instantiated
  // at compile time for each basis specialization.
  // Layer 2: unwrap the basis pointer to a concrete BasisOnTheFly<...> type.
  matrices::spinhalf::dispatch_basis(
      block_in, block_out, [&](auto const &basis_in, auto const &basis_out) {
        // Kernel: definition is in kernels.cpp, instantiated per basis type.
        matrices::matrix<matrices::spinhalf::MatrixPolicy>(ops, basis_in,
                                                           basis_out, mat);
      });
}
XDIAG_CATCH

template <typename coeff_t>
void matrix(OpSum const &ops, Boson const &block_in, Boson const &block_out,
            coeff_t *mat) try {
  matrices::boson::dispatch_basis(
      block_in, block_out, [&](auto const &basis_in, auto const &basis_out) {
        matrices::matrix<matrices::boson::MatrixPolicy>(ops, basis_in,
                                                        basis_out, mat);
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

  int64_t m = size(block_out);
  int64_t n = dim(block_in);
  arma::mat mat(m, n, arma::fill::zeros); // zeroed here; kernel uses +=
  // Layer 1: unwrap the Block variant; op_t is promoted to OpSum inside.
  utils::visit_same_type(
      block_in, block_out,
      [&](auto const &bin, auto const &bout) {
        matrix(OpSum(op), bin, bout, mat.memptr());
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
  // Layer 1: unwrap the Block variant; op_t is promoted to OpSum inside.
  utils::visit_same_type(
      block_in, block_out,
      [&](auto const &bin, auto const &bout) {
        matrix(OpSum(op), bin, bout, mat.memptr());
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

template void matrix(OpSum const &, Spinhalf const &, Spinhalf const &,
                     double *);
template void matrix(OpSum const &, Spinhalf const &, Spinhalf const &,
                     complex *);
template void matrix(OpSum const &, Boson const &, Boson const &, double *);
template void matrix(OpSum const &, Boson const &, Boson const &, complex *);

} // namespace xdiag
