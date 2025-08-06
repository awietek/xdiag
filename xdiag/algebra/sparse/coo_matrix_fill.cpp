// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "coo_matrix_fill.hpp"

#include <xdiag/algebra/apply_dispatch.hpp>
#include <xdiag/operators/logic/compilation.hpp>
#include <xdiag/utils/timing.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xdiag::algebra {

#ifdef _OPENMP

template <typename idx_t, typename coeff_t, typename block_t>
void coo_matrix_fill(OpSum const &ops, block_t const &block_in,
                     block_t const &block_out,
                     std::vector<int64_t> const &nnz_thread, int64_t size,
                     idx_t *row, idx_t *col, coeff_t *data, idx_t i0) try {
  auto t0 = rightnow();
  OpSum opsc = operators::compile<block_t>(ops);
  int nthreads = nnz_thread.size();
  std::vector<int64_t> nnz_thread_offset(nthreads, 0);
  std::exclusive_scan(nnz_thread.begin(), nnz_thread.end(),
                      nnz_thread_offset.begin(), 0);

  // Finally fill the COO matrix
  omp_set_schedule(omp_sched_static, 0);
  auto fill = [&](int64_t c, int64_t r, coeff_t d, int thread_num) {
    int64_t idx = nnz_thread_offset[thread_num]++;
    row[idx] = (idx_t)r;
    col[idx] = (idx_t)c;
    data[idx] = d;
  };
  algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill);
  timing(t0, rightnow(), "Sparse COO: matrix computation", 1);

  if (i0 == 1) {
    auto t0 = rightnow();
#pragma omp parallel for
    for (int64_t i = 0; i < size; ++i) {
      ++row[i];
      ++col[i];
    }
    timing(t0, rightnow(), "Sparse COO: adjusting 1-indexing", 1);
  }
}
XDIAG_CATCH

// int32_t, double
template void coo_matrix_fill(OpSum const &, Spinhalf const &, Spinhalf const &,
                              std::vector<int64_t> const &, int64_t, int32_t *,
                              int32_t *, double *, int32_t);
template void coo_matrix_fill(OpSum const &, tJ const &, tJ const &,
                              std::vector<int64_t> const &, int64_t, int32_t *,
                              int32_t *, double *, int32_t);
template void coo_matrix_fill(OpSum const &, Electron const &, Electron const &,
                              std::vector<int64_t> const &, int64_t, int32_t *,
                              int32_t *, double *, int32_t);

// int32_t, complex
template void coo_matrix_fill(OpSum const &, Spinhalf const &, Spinhalf const &,
                              std::vector<int64_t> const &, int64_t, int32_t *,
                              int32_t *, complex *, int32_t);
template void coo_matrix_fill(OpSum const &, tJ const &, tJ const &,
                              std::vector<int64_t> const &, int64_t, int32_t *,
                              int32_t *, complex *, int32_t);
template void coo_matrix_fill(OpSum const &, Electron const &, Electron const &,
                              std::vector<int64_t> const &, int64_t, int32_t *,
                              int32_t *, complex *, int32_t);

// int64_t, double
template void coo_matrix_fill(OpSum const &, Spinhalf const &, Spinhalf const &,
                              std::vector<int64_t> const &, int64_t, int64_t *,
                              int64_t *, double *, int64_t);
template void coo_matrix_fill(OpSum const &, tJ const &, tJ const &,
                              std::vector<int64_t> const &, int64_t, int64_t *,
                              int64_t *, double *, int64_t);
template void coo_matrix_fill(OpSum const &, Electron const &, Electron const &,
                              std::vector<int64_t> const &, int64_t, int64_t *,
                              int64_t *, double *, int64_t);

// int64_t, complex
template void coo_matrix_fill(OpSum const &, Spinhalf const &, Spinhalf const &,
                              std::vector<int64_t> const &, int64_t, int64_t *,
                              int64_t *, complex *, int64_t);
template void coo_matrix_fill(OpSum const &, tJ const &, tJ const &,
                              std::vector<int64_t> const &, int64_t, int64_t *,
                              int64_t *, complex *, int64_t);
template void coo_matrix_fill(OpSum const &, Electron const &, Electron const &,
                              std::vector<int64_t> const &, int64_t, int64_t *,
                              int64_t *, complex *, int64_t);

#else

template <typename idx_t, typename coeff_t, typename block_t>
void coo_matrix_fill(OpSum const &ops, block_t const &block_in,
                     block_t const &block_out, int64_t nnz, idx_t *row,
                     idx_t *col, coeff_t *data, idx_t i0) try {
  auto t0 = rightnow();
  OpSum opsc = operators::compile<block_t>(ops);
  int64_t idx = 0;
  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    row[idx] = (idx_t)idx_out;
    col[idx] = (idx_t)idx_in;
    data[idx] = val;
    ++idx;
  };
  algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill);
  timing(t0, rightnow(), "Sparse COO: matrix computation", 1);

  if (i0 == 1) {
    auto t0 = rightnow();
    for (int64_t i = 0; i < nnz; ++i) {
      ++row[i];
      ++col[i];
    }
    timing(t0, rightnow(), "Sparse COO: resource estimation", 1);
  }
}
XDIAG_CATCH

// int32_t, double
template void coo_matrix_fill(OpSum const &, Spinhalf const &, Spinhalf const &,
                              int64_t, int32_t *, int32_t *, double *, int32_t);
template void coo_matrix_fill(OpSum const &, tJ const &, tJ const &, int64_t,
                              int32_t *, int32_t *, double *, int32_t);
template void coo_matrix_fill(OpSum const &, Electron const &, Electron const &,
                              int64_t, int32_t *, int32_t *, double *, int32_t);

// int32_t, complex
template void coo_matrix_fill(OpSum const &, Spinhalf const &, Spinhalf const &,
                              int64_t, int32_t *, int32_t *, complex *,
                              int32_t);
template void coo_matrix_fill(OpSum const &, tJ const &, tJ const &, int64_t,
                              int32_t *, int32_t *, complex *, int32_t);
template void coo_matrix_fill(OpSum const &, Electron const &, Electron const &,
                              int64_t, int32_t *, int32_t *, complex *,
                              int32_t);

// int64_t, double
template void coo_matrix_fill(OpSum const &, Spinhalf const &, Spinhalf const &,
                              int64_t, int64_t *, int64_t *, double *, int64_t);
template void coo_matrix_fill(OpSum const &, tJ const &, tJ const &, int64_t,
                              int64_t *, int64_t *, double *, int64_t);
template void coo_matrix_fill(OpSum const &, Electron const &, Electron const &,
                              int64_t, int64_t *, int64_t *, double *, int64_t);

// int64_t, complex
template void coo_matrix_fill(OpSum const &, Spinhalf const &, Spinhalf const &,
                              int64_t, int64_t *, int64_t *, complex *,
                              int64_t);
template void coo_matrix_fill(OpSum const &, tJ const &, tJ const &, int64_t,
                              int64_t *, int64_t *, complex *, int64_t);
template void coo_matrix_fill(OpSum const &, Electron const &, Electron const &,
                              int64_t, int64_t *, int64_t *, complex *,
                              int64_t);

#endif

} // namespace xdiag::algebra
