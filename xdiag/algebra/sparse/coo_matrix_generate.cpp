// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "coo_matrix_generate.hpp"

#include <limits>
#include <xdiag/algebra/apply_dispatch.hpp>
#include <xdiag/operators/logic/block.hpp>
#include <xdiag/operators/logic/compilation.hpp>
#include <xdiag/operators/logic/hc.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/operators/logic/valid.hpp>

#include <xdiag/utils/timing.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xdiag::algebra {

template <typename idx_t, typename coeff_t, typename block_t>
COOMatrix<idx_t, coeff_t>
coo_matrix_generate(OpSum const &ops, block_t const &block_in,
                    block_t const &block_out, idx_t i0) try {

  if (!((i0 == 0) || (i0 == 1))) {
    XDIAG_THROW(fmt::format(
        "Invalid zero index i0. Must be either 0 or 1, but got i0={}", i0));
  }

  // Check if ops and blocks are compatible
  if (!blocks_match(ops, block_in, block_out)) {
    XDIAG_THROW(
        "Cannot create a COO matrix on Blocks. The resulting block is not in "
        "the correct symmetry sector. Please check the quantum numbers "
        "of the output block.");
  }

  // Check if real matrix can be created
  if constexpr (isreal<coeff_t>()) {
    if (!isreal(ops)) {
      XDIAG_THROW("Cannot create a real COO matrix from an Op or OpSum which "
                  "is complex.");
    }
    if (!isreal(block_in) || !isreal(block_out)) {
      XDIAG_THROW("Cannot create a real COO matrix when a block is complex.")
    }
  }

  int64_t nsites = block_in.nsites();
  check_valid(ops, nsites);

  int64_t m = block_out.size();
  int64_t n = block_in.size();
  if ((m > std::numeric_limits<idx_t>::max()) ||
      (n > std::numeric_limits<idx_t>::max())) {
    XDIAG_THROW(
        "Block is too large for index type which attempts to hold indices. "
        "Consider using a larger index type, e.g. 64 bit integers");
  }
  idx_t nrows = (idx_t)m;
  idx_t ncols = (idx_t)n;

  OpSum opsc = operators::compile<block_t>(ops);

  // We are doing two runs to create the matrix,
  // the parallel pattern needs to be identical
  // -> we choose static scheduling

  auto t0 = rightnow();
#ifdef _OPENMP
  omp_set_schedule(omp_sched_static, 0);

  // Get number of threads
  int nthreads = 0;
#pragma omp parallel
  {
#pragma omp single
    nthreads = omp_get_num_threads();
  }

  // compute number of non-zero elements on each thread
  std::vector<int64_t> nnz_thread(nthreads, 0);
  auto fill_count_nnz = [&](int64_t idx_in, int64_t idx_out, coeff_t val,
                            int thread_num) { ++nnz_thread[thread_num]; };
  algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill_count_nnz);

  std::vector<int64_t> nnz_thread_offset(nthreads, 0);
  std::exclusive_scan(nnz_thread.begin(), nnz_thread.end(),
                      nnz_thread_offset.begin(), 0);

  int64_t nnz = std::accumulate(nnz_thread.begin(), nnz_thread.end(), 0);
#else
  // In a first run count number of non-zero elements
  int64_t nnz = 0;
  auto fill_count_nnz = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    ++nnz;
  };
  algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill_count_nnz);
#endif
  timing(t0, rightnow(), "Sparse COO: resource estimation", 1);

  t0 = rightnow();
  // Try allocating vectors
  std::vector<idx_t> row;
  std::vector<idx_t> col;
  std::vector<coeff_t> data;
  try {
    row.resize(nnz, 0);
    col.resize(nnz, 0);
    data.resize(nnz, 0);
  } catch (...) {
    XDIAG_THROW("Cannot allocate row/col/data arrays for sparse matrix in "
                "coordinate (COO) format");
  }
  timing(t0, rightnow(), "Sparse COO: resource allocation", 1);

  t0 = rightnow();
  // Fill vectors with a second run
#ifdef _OPENMP
  // Finally fill the COO matrix
  auto fill = [&](int64_t c, int64_t r, coeff_t d, int thread_num) {
    int64_t idx = nnz_thread_offset[thread_num]++;
    row[idx] = (idx_t)r;
    col[idx] = (idx_t)c;
    data[idx] = d;
  };
#else
  int64_t idx = 0;
  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    row[idx] = (idx_t)idx_out;
    col[idx] = (idx_t)idx_in;
    data[idx] = val;
    ++idx;
  };
#endif

  algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill);
  timing(t0, rightnow(), "Sparse COO: matrix computation", 1);

  if (i0 == 1) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int64_t i = 0; i < (int64_t)data.size(); ++i) {
      ++row[i];
      ++col[i];
    }
  }

  return COOMatrix<idx_t, coeff_t>{nrows, ncols,           row, col, data,
                                   i0,    ishermitian(ops)};
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::algebra
