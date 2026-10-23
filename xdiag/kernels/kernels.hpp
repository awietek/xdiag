// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/armadillo.hpp>
#include <xdiag/kernels/elements.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>

// The two basis-level kernels, apply and elements, and the macros that
// instantiate them. Include ONLY from kernels/blocks/<block>/kernels.cpp:
// apply.cpp and elements.cpp see just their declarations, so a kernel is
// compiled only there, once per basis type (see the instantiation groups in
// cmake/instantiations.cmake). Each block source provides the
// matrix_kernel<Block> specialization, the ONLY place that pulls in that
// block's matrix_generic and term headers.
//
// Both kernels take a compiled ops, i.e. transpose of the operator (see
// compile.hpp), and run it from the rows (the output basis of the operator) to
// the columns.
// The term kernels split their input states, the rows, over the threads, so
// anything written per row needs no synchronisation.

namespace xdiag::kernels {

// Maps a block type to its block-specific matrix_generic kernel. The primary
// template is intentionally undefined; each block source specializes it.
template <typename block_t, typename coeff_t> struct matrix_kernel;

// The fill of the matrix-vector product. The kernels run transpose(ops), so for each
// row they enumerate they report val = ops(row, col), and this gathers
// vec_out[row] += ops(row, col) * vec_in[col]. Every row is written by exactly
// one thread, in a fixed order: no atomics, and the result is bitwise
// identical for any number of threads.
template <typename coeff_t>
inline void fill_apply(arma::Col<coeff_t> const &vec_in,
                       arma::Col<coeff_t> &vec_out, int64_t row, int64_t col,
                       coeff_t val) {
  vec_out[row] += val * vec_in[col];
}

template <typename coeff_t>
inline void fill_apply(arma::Mat<coeff_t> const &mat_in,
                       arma::Mat<coeff_t> &mat_out, int64_t row, int64_t col,
                       coeff_t val) {
  for (arma::uword i = 0; i < mat_in.n_cols; ++i) {
    mat_out.at(row, i) += val * mat_in.at(col, i);
  }
}

// vec_out += ops * vec_in. The fill stays an inlined lambda, so it is compiled
// into the term kernels of the matrix-vector product.
template <typename block_t, typename basis_t, typename mat_t>
void apply(OpSum const &ops, basis_t const &rows, mat_t const &vec_in,
           basis_t const &cols, mat_t &vec_out) try {
  using coeff_t = typename mat_t::elem_type;
  matrix_kernel<block_t, coeff_t>::call(
      ops, rows, cols, [&](int64_t row, int64_t col, coeff_t val) {
        fill_apply(vec_in, vec_out, row, col, val);
      });
}
XDIAG_CATCH

// Calls fill(row, col, ops(row, col)) for every nonzero matrix element. The
// kernels report the elements of trans(ops), ops(row, col). One
// instantiation per coeff_t serves every builder (dense, CSR, CSC, COO).
template <typename block_t, typename coeff_t, typename basis_t>
void elements(OpSum const &ops, basis_t const &rows, basis_t const &cols,
              fill_t<coeff_t> const &fill) try {
  matrix_kernel<block_t, coeff_t>::call(
      ops, rows, cols, [&](int64_t row, int64_t col, coeff_t val) {
        fill(row, col, val);
      });
}
XDIAG_CATCH

} // namespace xdiag::kernels

// ---------------------------------------------------------------------------
// Explicit-instantiation helper macros, used by the block sources. The bare
// type names in the expansions (vec, cx_mat, complex, ...) are resolved at the
// expansion site, so each source must bring the relevant namespaces into scope
// (using namespace arma; xdiag; xdiag::basis; ...).
// ---------------------------------------------------------------------------

#define XDIAG_INSTANTIATE_APPLY(BLOCK, BASIS, MAT)                             \
  template void xdiag::kernels::apply<BLOCK, BASIS, MAT>(                      \
      OpSum const &, BASIS const &, MAT const &, BASIS const &, MAT &);

#define XDIAG_INSTANTIATE_ELEMENTS(BLOCK, BASIS, COEFF)                        \
  template void xdiag::kernels::elements<BLOCK, COEFF, BASIS>(                 \
      OpSum const &, BASIS const &, BASIS const &,                             \
      xdiag::kernels::fill_t<COEFF> const &);

#define XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS)                                \
  XDIAG_INSTANTIATE_APPLY(BLOCK, BASIS, vec)                                   \
  XDIAG_INSTANTIATE_APPLY(BLOCK, BASIS, cx_vec)                                \
  XDIAG_INSTANTIATE_APPLY(BLOCK, BASIS, mat)                                   \
  XDIAG_INSTANTIATE_APPLY(BLOCK, BASIS, cx_mat)                                \
  XDIAG_INSTANTIATE_ELEMENTS(BLOCK, BASIS, double)                             \
  XDIAG_INSTANTIATE_ELEMENTS(BLOCK, BASIS, complex)

// Boson bases share a set of BitArray / BitArrayLong backends for a given
// (basis class, enumeration). Only the widths {1,2,3,4,8} are compiled; states
// needing 5-7 bits per site are promoted to the 8-bit backend at dispatch time
// (see promote_nlocalbits in blocks/boson.cpp), so those widths are never
// constructed and need no kernels.
#define XDIAG_INSTANTIATE_KERNELS_BITARRAY(BLOCK, BASIS, ENUM)                 \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArray1>>)                     \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArray2>>)                     \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArray3>>)                     \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArray4>>)                     \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArray8>>)

#define XDIAG_INSTANTIATE_KERNELS_BITARRAY_LONG(BLOCK, BASIS, ENUM)            \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArrayLong1>>)                 \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArrayLong2>>)                 \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArrayLong3>>)                 \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArrayLong4>>)                 \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArrayLong8>>)
