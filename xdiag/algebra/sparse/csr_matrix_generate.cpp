// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "csr_matrix_generate.hpp"

#include <limits>
#include <numeric>

#include <xdiag/operators/logic/block.hpp>
#include <xdiag/operators/logic/compilation.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/operators/logic/valid.hpp>
#include <xdiag/utils/timing.hpp>

#include <xdiag/algebra/apply_dispatch.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xdiag::algebra {

template <typename idx_t, typename coeff_t, typename block_t>
CSRMatrix<idx_t, coeff_t>
csr_matrix_generate(OpSum const &ops, block_t const &block_in,
                    block_t const &block_out, idx_t i0, bool transpose) try {
  std::string mattype = transpose ? "CSC" : "CSR";

  if (!((i0 == 0) || (i0 == 1))) {
    XDIAG_THROW(fmt::format(
        "Invalid zero index i0. Must be either 0 or 1, but got i0={}", i0));
  }

  // Check if ops and blocks are compatible
  if (!blocks_match(ops, block_in, block_out)) {
    XDIAG_THROW(fmt::format(
        "Cannot create a {} matrix on Blocks. The resulting block is not in "
        "the correct symmetry sector. Please check the quantum numbers "
        "of the output block.",
        mattype));
  }

  // Check if real matrix can be created
  if constexpr (isreal<coeff_t>()) {
    if (!isreal(ops)) {
      XDIAG_THROW(fmt::format(
          "Cannot create a real {} matrix from an Op or OpSum which "
          "is complex.",
          mattype));
    }
    if (!isreal(block_in) || !isreal(block_out)) {
      XDIAG_THROW(fmt::format(
          "Cannot create a real {} matrix when a block is complex.", mattype));
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

#ifdef _OPENMP
  omp_set_schedule(omp_sched_static, 0);

  // Get number of threads
  int nthreads = 0;
#pragma omp parallel
  {
#pragma omp single
    nthreads = omp_get_num_threads();
  }
#endif

  // Try allocating the row pointer array
  std::vector<idx_t> rowptr;
  auto t0 = rightnow();
  try {
    rowptr.resize(transpose ? ncols + 1 : nrows + 1, 0);
  } catch (...) {
    XDIAG_THROW(
        fmt::format("Cannot allocate row pointer array for sparse matrix in "
                    "{} format",
                    mattype));
  }
  timing(t0, rightnow(),
         fmt::format("Sparse {}: row pointer allocation", mattype), 1);

  // Count number of elements in row
  t0 = rightnow();
  std::vector<idx_t> &n_elements_in_row = rowptr;

  if (transpose) {

    auto fill_count_row = [&](int64_t c, int64_t r, coeff_t val) {
#ifdef _OPENMP
#pragma omp atomic update
#endif
      ++n_elements_in_row[c];
    };
    algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill_count_row);

  } else { // not transposed
    auto fill_count_row = [&](int64_t c, int64_t r, coeff_t val) {
#ifdef _OPENMP
#pragma omp atomic update
#endif
      ++n_elements_in_row[r];
    };
    algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill_count_row);
  }
  timing(t0, rightnow(), fmt::format("Sparse {}: resource estimation", mattype),
         1);

  t0 = rightnow();
  int64_t nnz =
      std::accumulate(n_elements_in_row.begin(), n_elements_in_row.end(), 0);
  int64_t max_entries_in_row =
      *std::max_element(n_elements_in_row.begin(), n_elements_in_row.end());

  std::vector<idx_t> n_elements_in_row_offset;
  try {
    n_elements_in_row_offset.resize(rowptr.size(), 0);
  } catch (...) {
    XDIAG_THROW(fmt::format(
        "Cannot allocate auxiliary row pointer array for sparse matrix in "
        "{} format",
        mattype));
  }
  std::exclusive_scan(n_elements_in_row.begin(), n_elements_in_row.end(),
                      n_elements_in_row_offset.begin(), 0);
  rowptr = n_elements_in_row_offset;
  if (transpose) {
    rowptr[ncols] = nnz;
  } else {
    rowptr[nrows] = nnz;
  }
  timing(t0, rightnow(), fmt::format("Sparse {}: reductions", mattype), 1);

  t0 = rightnow();
  // Try allocating vectors
  std::vector<idx_t> col;
  std::vector<coeff_t> data;
  try {
    col.resize(nnz, 0);
    data.resize(nnz, 0);
  } catch (...) {
    XDIAG_THROW(
        fmt::format("Cannot allocate col/data arrays for sparse matrix in "
                    "{} format",
                    mattype));
  }
  timing(t0, rightnow(), fmt::format("Sparse {}: resource allocation", mattype),
         1);

  t0 = rightnow();

  // Fill vectors with a second Finally
  // run fill the CSR matrix
  // std::vector<int64_t> idx(nthreads, 0);
  if (transpose) {
    auto fill = [&](int64_t c, int64_t r, coeff_t d) {
#ifdef _OPENMP
      idx_t idx;
#pragma omp atomic capture
      idx = n_elements_in_row_offset[c]++;
#else
      idx_t idx = n_elements_in_row_offset[c]++;
#endif
      col[idx] = (idx_t)r;
      data[idx] = d;
    };
    algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill);
  } else { // !transpose
    auto fill = [&](int64_t c, int64_t r, coeff_t d) {
#ifdef _OPENMP
      idx_t idx;
#pragma omp atomic capture
      idx = n_elements_in_row_offset[r]++;
#else
      idx_t idx = n_elements_in_row_offset[r]++;
#endif
      col[idx] = (idx_t)c;
      data[idx] = d;
    };
    algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill);
  }
  timing(t0, rightnow(), fmt::format("Sparse {}: matrix assembly", mattype), 1);

  t0 = rightnow();

#ifdef _OPENMP
#pragma omp parallel
  {
#endif
    std::vector<idx_t> indices(max_entries_in_row, 0);
    std::vector<idx_t> coltmp(max_entries_in_row, 0);
    std::vector<coeff_t> datatmp(max_entries_in_row, 0);

    idx_t nrowst = transpose ? ncols : nrows;

#ifdef _OPENMP
#pragma omp for
#endif
    for (idx_t row = 0; row < nrowst; ++row) {
      idx_t start = rowptr[row];
      idx_t end = rowptr[row + 1];
      idx_t nelems = end - start;

      std::copy(col.begin() + start, col.begin() + end, coltmp.begin());
      std::copy(data.begin() + start, data.begin() + end, datatmp.begin());
      std::iota(indices.begin(), indices.begin() + nelems, 0);
      std::sort(indices.begin(), indices.begin() + nelems,
                [&](idx_t i, idx_t j) { return coltmp[i] < coltmp[j]; });

      for (idx_t i = 0; i < nelems; ++i) {
        col[start + i] = coltmp[indices[i]];
        data[start + i] = datatmp[indices[i]];
      }
    }
#ifdef _OPENMP
  }
#endif
  timing(t0, rightnow(), fmt::format("Sparse {}: sorting", mattype), 1);

  if (i0 == 1) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int64_t i = 0; i < (int64_t)rowptr.size(); ++i) {
      ++rowptr[i];
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int64_t i = 0; i < (int64_t)col.size(); ++i) {
      ++col[i];
    }
  }

  return CSRMatrix<idx_t, coeff_t>{transpose ? ncols : nrows,
                                   transpose ? nrows : ncols,
                                   rowptr,
                                   col,
                                   data,
                                   i0};
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::algebra
