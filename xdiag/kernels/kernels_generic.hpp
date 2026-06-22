// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <numeric>

#include <xdiag/armadillo.hpp>
#include <xdiag/kernels/fill_functions.hpp>
#include <xdiag/kernels/kernels.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

// Generic kernel definitions, shared by the per-block instantiation sources
// (kernels/blocks/<block>/kernels_<block>.cpp). Each of those sources includes
// this header, provides the matrix_kernel<Block> specialization for its own
// block (which is the ONLY place that pulls in that block's matrix_generic and
// term headers), and then explicitly instantiates the kernels for that block's
// basis types. Keeping the block-specific matrix_generic out of this header is
// what lets a change to one block's terms recompile only that block's kernels.

namespace xdiag::kernels {

// Maps a block type to its block-specific matrix_generic kernel. The primary
// template is intentionally undefined; each per-block kernel source provides the
// specialization for its own block before instantiating the kernels below.
template <typename block_t> struct matrix_kernel;

template <typename block_t, typename basis_t, typename mat_t>
void apply(OpSum const &ops, basis_t const &basis_in, mat_t const &mat_in,
           basis_t const &basis_out, mat_t &mat_out) try {
  using coeff_t = typename mat_t::elem_type;
  mat_out.zeros();
  matrix_kernel<block_t>::template call<coeff_t>(
      ops, basis_in, basis_out,
      [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
        fill_apply(mat_in, mat_out, idx_in, idx_out, val);
      });
}
XDIAG_CATCH

template <typename block_t, typename coeff_t, typename basis_t>
void matrix(OpSum const &ops, basis_t const &basis_in, basis_t const &basis_out,
            coeff_t *mat) try {
  int64_t m = basis_out.size();
  matrix_kernel<block_t>::template call<coeff_t>(
      ops, basis_in, basis_out,
      [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
        fill_matrix(mat, m, idx_in, idx_out, val);
      });
}
XDIAG_CATCH

#ifdef _OPENMP

template <typename block_t, typename coeff_t, typename basis_t>
std::vector<int64_t> coo_matrix_nnz(OpSum const &ops, basis_t const &basis_in,
                                    basis_t const &basis_out) try {
  omp_set_schedule(omp_sched_static, 0);
  int nthreads = omp_get_max_threads();
  std::vector<int64_t> nnz_thread(nthreads, 0);
  matrix_kernel<block_t>::template call<coeff_t>(
      ops, basis_in, basis_out, [&](int64_t, int64_t, coeff_t, int num_thread) {
        fill_coo_count(nnz_thread, num_thread);
      });
  return nnz_thread;
}
XDIAG_CATCH

template <typename block_t, typename coeff_t, typename basis_t,
          typename idx_t>
void coo_matrix_fill(OpSum const &ops, basis_t const &basis_in,
                     basis_t const &basis_out,
                     std::vector<int64_t> const &nnz_thread, idx_t *rows,
                     idx_t *cols, coeff_t *data, idx_t i0) try {
  int nthreads = (int)nnz_thread.size();
  std::vector<int64_t> offset(nthreads, 0);
  std::partial_sum(nnz_thread.begin(), nnz_thread.end() - 1,
                   offset.begin() + 1);
  omp_set_schedule(omp_sched_static, 0);
  matrix_kernel<block_t>::template call<coeff_t>(
      ops, basis_in, basis_out,
      [&](int64_t idx_in, int64_t idx_out, coeff_t val, int num_thread) {
        fill_coo(rows, cols, data, offset, idx_in, idx_out, val, i0,
                 num_thread);
      });
}
XDIAG_CATCH

#else // !_OPENMP

template <typename block_t, typename coeff_t, typename basis_t>
int64_t coo_matrix_nnz(OpSum const &ops, basis_t const &basis_in,
                       basis_t const &basis_out) try {
  int64_t nnz = 0;
  matrix_kernel<block_t>::template call<coeff_t>(
      ops, basis_in, basis_out,
      [&](int64_t, int64_t, coeff_t) { fill_coo_count(nnz); });
  return nnz;
}
XDIAG_CATCH

template <typename block_t, typename coeff_t, typename basis_t,
          typename idx_t>
void coo_matrix_fill(OpSum const &ops, basis_t const &basis_in,
                     basis_t const &basis_out, idx_t *rows, idx_t *cols,
                     coeff_t *data, idx_t i0) try {
  int64_t k = 0;
  matrix_kernel<block_t>::template call<coeff_t>(
      ops, basis_in, basis_out,
      [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
        fill_coo(rows, cols, data, k, idx_in, idx_out, val, i0);
      });
}
XDIAG_CATCH

#endif // _OPENMP

template <typename block_t, typename coeff_t, typename basis_t>
std::vector<int64_t> csr_matrix_nnz(OpSum const &ops, basis_t const &basis_in,
                                    basis_t const &basis_out,
                                    bool transpose) try {
  if (transpose) {
    std::vector<int64_t> n_elements(basis_in.size(), 0);
    matrix_kernel<block_t>::template call<coeff_t>(
        ops, basis_in, basis_out, [&](int64_t idx_in, int64_t, coeff_t) {
          fill_csr_count(n_elements, idx_in);
        });
    return n_elements;
  } else {
    std::vector<int64_t> n_elements(basis_out.size(), 0);
    matrix_kernel<block_t>::template call<coeff_t>(
        ops, basis_in, basis_out, [&](int64_t, int64_t idx_out, coeff_t) {
          fill_csr_count(n_elements, idx_out);
        });
    return n_elements;
  }
}
XDIAG_CATCH

template <typename block_t, typename coeff_t, typename basis_t,
          typename idx_t>
void csr_matrix_fill(OpSum const &ops, basis_t const &basis_in,
                     basis_t const &basis_out, std::vector<int64_t> &offset,
                     idx_t *col, coeff_t *data, idx_t i0, bool transpose) try {
  if (transpose) {
    // CSC: key by column (idx_in), store row (idx_out)
    matrix_kernel<block_t>::template call<coeff_t>(
        ops, basis_in, basis_out,
        [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
          fill_csr(offset, col, data, idx_out, idx_in, val, i0);
        });
  } else {
    // CSR: key by row (idx_out), store column (idx_in)
    matrix_kernel<block_t>::template call<coeff_t>(
        ops, basis_in, basis_out,
        [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
          fill_csr(offset, col, data, idx_in, idx_out, val, i0);
        });
  }
}
XDIAG_CATCH

} // namespace xdiag::kernels

// ---------------------------------------------------------------------------
// Explicit-instantiation helper macros, used by the per-block kernel sources.
// The bare type names in the expansions (vec, cx_mat, complex, ...) are resolved
// at the expansion site, so each source must bring the relevant namespaces into
// scope (using namespace arma; xdiag; xdiag::basis; ...).
// ---------------------------------------------------------------------------

#define XDIAG_INSTANTIATE_APPLY(BLOCK, BASIS, MAT)                                   \
  template void xdiag::kernels::apply<BLOCK, BASIS, MAT>(                     \
      OpSum const &, BASIS const &, MAT const &, BASIS const &, MAT &);

#define XDIAG_INSTANTIATE_MATRIX(BLOCK, BASIS, COEFF)                                \
  template void xdiag::kernels::matrix<BLOCK, COEFF, BASIS>(                  \
      OpSum const &, BASIS const &, BASIS const &, COEFF *);

#ifdef _OPENMP
#define XDIAG_INSTANTIATE_COO_NNZ(BLOCK, BASIS, COEFF)                               \
  template std::vector<int64_t>                                                \
  xdiag::kernels::coo_matrix_nnz<BLOCK, COEFF, BASIS>(                        \
      OpSum const &, BASIS const &, BASIS const &);
#define XDIAG_INSTANTIATE_COO_FILL(BLOCK, BASIS, IDX, COEFF)                         \
  template void                                                                \
  xdiag::kernels::coo_matrix_fill<BLOCK, COEFF, BASIS, IDX>(                  \
      OpSum const &, BASIS const &, BASIS const &,                            \
      std::vector<int64_t> const &, IDX *, IDX *, COEFF *, IDX);
#else
#define XDIAG_INSTANTIATE_COO_NNZ(BLOCK, BASIS, COEFF)                               \
  template int64_t                                                             \
  xdiag::kernels::coo_matrix_nnz<BLOCK, COEFF, BASIS>(                        \
      OpSum const &, BASIS const &, BASIS const &);
#define XDIAG_INSTANTIATE_COO_FILL(BLOCK, BASIS, IDX, COEFF)                         \
  template void                                                                \
  xdiag::kernels::coo_matrix_fill<BLOCK, COEFF, BASIS, IDX>(                  \
      OpSum const &, BASIS const &, BASIS const &, IDX *, IDX *, COEFF *,      \
      IDX);
#endif

#define XDIAG_INSTANTIATE_CSR_NNZ(BLOCK, BASIS, COEFF)                               \
  template std::vector<int64_t>                                                \
  xdiag::kernels::csr_matrix_nnz<BLOCK, COEFF, BASIS>(                        \
      OpSum const &, BASIS const &, BASIS const &, bool);
#define XDIAG_INSTANTIATE_CSR_FILL(BLOCK, BASIS, IDX, COEFF)                         \
  template void                                                                \
  xdiag::kernels::csr_matrix_fill<BLOCK, COEFF, BASIS, IDX>(                  \
      OpSum const &, BASIS const &, BASIS const &, std::vector<int64_t> &,     \
      IDX *, COEFF *, IDX, bool);

#define XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS)                                      \
  XDIAG_INSTANTIATE_APPLY(BLOCK, BASIS, vec)                                         \
  XDIAG_INSTANTIATE_APPLY(BLOCK, BASIS, cx_vec)                                      \
  XDIAG_INSTANTIATE_APPLY(BLOCK, BASIS, mat)                                         \
  XDIAG_INSTANTIATE_APPLY(BLOCK, BASIS, cx_mat)                                      \
  XDIAG_INSTANTIATE_MATRIX(BLOCK, BASIS, double)                                     \
  XDIAG_INSTANTIATE_MATRIX(BLOCK, BASIS, complex)                                    \
  XDIAG_INSTANTIATE_COO_NNZ(BLOCK, BASIS, double)                                    \
  XDIAG_INSTANTIATE_COO_NNZ(BLOCK, BASIS, complex)                                   \
  XDIAG_INSTANTIATE_COO_FILL(BLOCK, BASIS, int32_t, double)                          \
  XDIAG_INSTANTIATE_COO_FILL(BLOCK, BASIS, int32_t, complex)                         \
  XDIAG_INSTANTIATE_COO_FILL(BLOCK, BASIS, int64_t, double)                          \
  XDIAG_INSTANTIATE_COO_FILL(BLOCK, BASIS, int64_t, complex)                         \
  XDIAG_INSTANTIATE_CSR_NNZ(BLOCK, BASIS, double)                                    \
  XDIAG_INSTANTIATE_CSR_NNZ(BLOCK, BASIS, complex)                                   \
  XDIAG_INSTANTIATE_CSR_FILL(BLOCK, BASIS, int32_t, double)                          \
  XDIAG_INSTANTIATE_CSR_FILL(BLOCK, BASIS, int32_t, complex)                         \
  XDIAG_INSTANTIATE_CSR_FILL(BLOCK, BASIS, int64_t, double)                          \
  XDIAG_INSTANTIATE_CSR_FILL(BLOCK, BASIS, int64_t, complex)

// Boson bases share the same BitArray<1..8> / BitArrayLong<1..8> backends for a
// given (basis class, enumeration), so instantiate all eight at once.
#define XDIAG_INSTANTIATE_KERNELS_BITARRAY(BLOCK, BASIS, ENUM)                       \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArray1>>)                           \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArray2>>)                           \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArray3>>)                           \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArray4>>)                           \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArray5>>)                           \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArray6>>)                           \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArray7>>)                           \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArray8>>)

#define XDIAG_INSTANTIATE_KERNELS_BITARRAY_LONG(BLOCK, BASIS, ENUM)                  \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArrayLong1>>)                       \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArrayLong2>>)                       \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArrayLong3>>)                       \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArrayLong4>>)                       \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArrayLong5>>)                       \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArrayLong6>>)                       \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArrayLong7>>)                       \
  XDIAG_INSTANTIATE_KERNELS(BLOCK, BASIS<ENUM<BitArrayLong8>>)
