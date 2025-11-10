// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "coo_matrix_generate.hpp"

#include <xdiag/algebra/sparse/coo_matrix_fill.hpp>
#include <xdiag/algebra/sparse/coo_matrix_nnz.hpp>
#include <xdiag/algebra/sparse/valid.hpp>
#include <xdiag/operators/logic/hc.hpp>
#include <xdiag/utils/timing.hpp>
#include <xdiag/utils/xdiag_show.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xdiag::algebra {

template <typename idx_t, typename coeff_t, typename block_t>
COOMatrix<idx_t, coeff_t>
coo_matrix_generate(OpSum const &ops, block_t const &block_in,
                    block_t const &block_out, idx_t i0) try {
  check_valid_sparse_matrix<idx_t, coeff_t>(ops, block_in, block_out, i0);
  idx_t nrows = (idx_t)size(block_out);
  idx_t ncols = (idx_t)size(block_in);
  if ((nrows == 0) || (ncols == 0)) {
    return COOMatrix<idx_t, coeff_t>();
  }

  // We are doing two runs to create the matrix,
  // the parallel pattern needs to be identical
  // -> we choose static scheduling
#ifdef _OPENMP
  std::vector<int64_t> nnz_thread =
      coo_matrix_nnz_thread<coeff_t>(ops, block_in, block_out);
  int64_t nnz = std::accumulate(nnz_thread.begin(), nnz_thread.end(), (int64_t)0);
#else
  int64_t nnz = coo_matrix_nnz<coeff_t>(ops, block_in, block_out);
#endif

  auto t0 = rightnow();

  // Try allocating vectors
  arma::Col<idx_t> row;
  arma::Col<idx_t> col;
  arma::Col<coeff_t> data;
  try {
    row.resize(nnz);
    col.resize(nnz);
    data.resize(nnz);
    row.zeros();
    col.zeros();
    data.zeros();
  } catch (...) {
    XDIAG_THROW("Cannot allocate row/col/data arrays for sparse matrix in "
                "coordinate (COO) format");
  }
  timing(t0, rightnow(), "Sparse COO: resource allocation", 1);

  // Fill vectors with a second run
#ifdef _OPENMP
  coo_matrix_fill(ops, block_in, block_out, nnz_thread, nnz, row.memptr(),
                  col.memptr(), data.memptr(), i0);
#else
  coo_matrix_fill(ops, block_in, block_out, nnz, row.memptr(), col.memptr(),
                  data.memptr(), i0);
#endif

  bool isherm = ishermitian(ops);
  return COOMatrix<idx_t, coeff_t>{nrows, ncols, row, col, data, i0, isherm};
}
XDIAG_CATCH

// int32_t, double
template COOMatrix<int32_t, double>
coo_matrix_generate<int32_t, double, Spinhalf>(OpSum const &, Spinhalf const &,
                                               Spinhalf const &, int32_t);
template COOMatrix<int32_t, double>
coo_matrix_generate<int32_t, double, tJ>(OpSum const &, tJ const &, tJ const &,
                                         int32_t);
template COOMatrix<int32_t, double>
coo_matrix_generate<int32_t, double, Electron>(OpSum const &, Electron const &,
                                               Electron const &, int32_t);

// int32_t, complex
template COOMatrix<int32_t, complex>
coo_matrix_generate<int32_t, complex, Spinhalf>(OpSum const &, Spinhalf const &,
                                                Spinhalf const &, int32_t);
template COOMatrix<int32_t, complex>
coo_matrix_generate<int32_t, complex, tJ>(OpSum const &, tJ const &, tJ const &,
                                          int32_t);
template COOMatrix<int32_t, complex>
coo_matrix_generate<int32_t, complex, Electron>(OpSum const &, Electron const &,
                                                Electron const &, int32_t);

// int64_t, double
template COOMatrix<int64_t, double>
coo_matrix_generate<int64_t, double, Spinhalf>(OpSum const &, Spinhalf const &,
                                               Spinhalf const &, int64_t);
template COOMatrix<int64_t, double>
coo_matrix_generate<int64_t, double, tJ>(OpSum const &, tJ const &, tJ const &,
                                         int64_t);
template COOMatrix<int64_t, double>
coo_matrix_generate<int64_t, double, Electron>(OpSum const &, Electron const &,
                                               Electron const &, int64_t);

// int64_t, complex
template COOMatrix<int64_t, complex>
coo_matrix_generate<int64_t, complex, Spinhalf>(OpSum const &, Spinhalf const &,
                                                Spinhalf const &, int64_t);
template COOMatrix<int64_t, complex>
coo_matrix_generate<int64_t, complex, tJ>(OpSum const &, tJ const &, tJ const &,
                                          int64_t);
template COOMatrix<int64_t, complex>
coo_matrix_generate<int64_t, complex, Electron>(OpSum const &, Electron const &,
                                                Electron const &, int64_t);

} // namespace xdiag::algebra
