// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "csr_matrix_generate.hpp"

#include <numeric>
#include <xdiag/algebra/sparse/csr_matrix_fill.hpp>
#include <xdiag/algebra/sparse/csr_matrix_nnz.hpp>
#include <xdiag/algebra/sparse/valid.hpp>
#include <xdiag/operators/logic/hc.hpp>
#include <xdiag/utils/timing.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xdiag::algebra {

template <typename idx_t, typename coeff_t, typename block_t>
CSRMatrix<idx_t, coeff_t>
csr_matrix_generate(OpSum const &ops, block_t const &block_in,
                    block_t const &block_out, idx_t i0, bool transpose) try {
  check_valid_sparse_matrix<idx_t, coeff_t>(ops, block_in, block_out, i0);
  idx_t nrows = (idx_t)size(block_out);
  idx_t ncols = (idx_t)size(block_in);
  if ((nrows == 0) || (ncols == 0)) {
    return CSRMatrix<idx_t, coeff_t>{
        0,
        0,
        arma::Col<idx_t>{0}, // rowprt / colptr must be size 1
        arma::Col<idx_t>(),
        arma::Col<coeff_t>(),
        i0,
        ishermitian(ops)};
  }

  // first count the number of nonzero elements in each row and in total
  std::vector<idx_t> n_elements_in_row =
      csr_matrix_nnz<idx_t, coeff_t>(ops, block_in, block_out, transpose);
  int64_t nnz =
      std::accumulate(n_elements_in_row.begin(), n_elements_in_row.end(), (int64_t)0);

  // try allocating the main memory
  auto t0 = rightnow();
  arma::Col<idx_t> rowptr;
  arma::Col<idx_t> col;
  arma::Col<coeff_t> data;
  try {
    rowptr.resize(n_elements_in_row.size() + 1);
    col.resize(nnz);
    data.resize(nnz);
    rowptr.zeros();
    col.zeros();
    data.zeros();
  } catch (std::exception& e) {
    std::cout << "nnz: " << nnz <<std::endl;
    std::cerr << "Armadillo error: " << e.what() << std::endl;
    XDIAG_THROW("Cannot allocate rowptr/col/data arrays for sparse matrix");
  }
  timing(t0, rightnow(), "Sparse matrix generation: resource allocation", 1);

  // Finally, fill the col/data arrays and create the rowptr
  csr_matrix_fill(ops, block_in, block_out, nnz, n_elements_in_row,
                  rowptr.memptr(), col.memptr(), data.memptr(), i0, transpose);

  return CSRMatrix<idx_t, coeff_t>{transpose ? ncols : nrows,
                                   transpose ? nrows : ncols,
                                   rowptr,
                                   col,
                                   data,
                                   i0,
                                   ishermitian(ops)};
}
XDIAG_CATCH

// int64_t, double
template CSRMatrix<int64_t, double> csr_matrix_generate(OpSum const &,
                                                        Spinhalf const &,
                                                        Spinhalf const &,
                                                        int64_t, bool);
template CSRMatrix<int64_t, double>
csr_matrix_generate(OpSum const &, tJ const &, tJ const &, int64_t, bool);
template CSRMatrix<int64_t, double> csr_matrix_generate(OpSum const &,
                                                        Electron const &,
                                                        Electron const &,
                                                        int64_t, bool);
// int64_t, complex
template CSRMatrix<int64_t, complex> csr_matrix_generate(OpSum const &,
                                                         Spinhalf const &,
                                                         Spinhalf const &,
                                                         int64_t, bool);
template CSRMatrix<int64_t, complex>
csr_matrix_generate(OpSum const &, tJ const &, tJ const &, int64_t, bool);
template CSRMatrix<int64_t, complex> csr_matrix_generate(OpSum const &,
                                                         Electron const &,
                                                         Electron const &,
                                                         int64_t, bool);

// int32_t, double
template CSRMatrix<int32_t, double> csr_matrix_generate(OpSum const &,
                                                        Spinhalf const &,
                                                        Spinhalf const &,
                                                        int32_t, bool);
template CSRMatrix<int32_t, double>
csr_matrix_generate(OpSum const &, tJ const &, tJ const &, int32_t, bool);
template CSRMatrix<int32_t, double> csr_matrix_generate(OpSum const &,
                                                        Electron const &,
                                                        Electron const &,
                                                        int32_t, bool);
// int32_t, complex
template CSRMatrix<int32_t, complex> csr_matrix_generate(OpSum const &,
                                                         Spinhalf const &,
                                                         Spinhalf const &,
                                                         int32_t, bool);
template CSRMatrix<int32_t, complex>
csr_matrix_generate(OpSum const &, tJ const &, tJ const &, int32_t, bool);
template CSRMatrix<int32_t, complex> csr_matrix_generate(OpSum const &,
                                                         Electron const &,
                                                         Electron const &,
                                                         int32_t, bool);

} // namespace xdiag::algebra
