// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "csc_matrix.hpp"

#include <algorithm>
#include <numeric>

#include <xdiag/armadillo.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/matrices/kernels.hpp>
#include <xdiag/matrices/spinhalf/dispatch_basis.hpp>
#include <xdiag/matrices/spinhalf/matrix_policy.hpp>
#include <xdiag/algebra/hc.hpp>
#include <xdiag/operators/qns/block.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/variants.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xdiag {

// Layer 2: Spinhalf-specific function.  Builds CSC via transposed CSR kernels.
template <typename idx_t, typename coeff_t>
static CSCMatrix<idx_t, coeff_t>
csc_matrix(OpSum const &ops, Spinhalf const &block_in,
           Spinhalf const &block_out, idx_t i0) try {
  idx_t nrows = (idx_t)size(block_out);
  idx_t ncols = (idx_t)size(block_in);
  arma::Col<idx_t> colptr, row;
  arma::Col<coeff_t> data;

  matrices::spinhalf::dispatch_basis(
      block_in, block_out, [&](auto const &basis_in, auto const &basis_out) {
        // Pass 1: count nonzeros per column (transpose=true keys by idx_in).
        auto n_elements_in_col =
            matrices::csr_matrix_nnz<matrices::spinhalf::MatrixPolicy, coeff_t>(
                ops, basis_in, basis_out, true);
        int64_t nnz = std::accumulate(n_elements_in_col.begin(),
                                      n_elements_in_col.end(), (int64_t)0);

        // Build colptr (exclusive prefix sum + i0 shift) and the mutable
        // offset array that csr_matrix_fill uses to claim write slots.
        colptr.resize(ncols + 1);
        row.resize(nnz);
        data.resize(nnz);

        std::vector<int64_t> offset(ncols, 0);
        if (ncols > 0) {
          std::partial_sum(n_elements_in_col.begin(),
                           n_elements_in_col.end() - 1, offset.begin() + 1);
          for (idx_t c = 0; c < ncols; ++c)
            colptr[c] = (idx_t)(offset[c] + i0);
          colptr[ncols] = (idx_t)(nnz + i0);
        }

        // Pass 2: fill row and data using atomic slot assignment.
        matrices::csr_matrix_fill<matrices::spinhalf::MatrixPolicy, coeff_t>(
            ops, basis_in, basis_out, offset, row.memptr(), data.memptr(), i0,
            true);

        // Sort each column's entries by row index (required for CSC validity).
        int64_t max_elems = n_elements_in_col.empty()
                                ? 0
                                : *std::max_element(n_elements_in_col.begin(),
                                                    n_elements_in_col.end());
#ifdef _OPENMP
#pragma omp parallel
        {
          std::vector<int64_t> indices(max_elems);
          std::vector<idx_t> rowtmp(max_elems);
          std::vector<coeff_t> datatmp(max_elems);
#pragma omp for
#else
        {
          std::vector<int64_t> indices(max_elems);
          std::vector<idx_t> rowtmp(max_elems);
          std::vector<coeff_t> datatmp(max_elems);
#endif
          for (idx_t col = 0; col < ncols; ++col) {
            int64_t start = (int64_t)(colptr[col] - i0);
            int64_t end = (int64_t)(colptr[col + 1] - i0);
            int64_t nelems = end - start;
            std::copy(row.memptr() + start, row.memptr() + end, rowtmp.begin());
            std::copy(data.memptr() + start, data.memptr() + end,
                      datatmp.begin());
            std::iota(indices.begin(), indices.begin() + nelems, (int64_t)0);
            std::sort(indices.begin(), indices.begin() + nelems,
                      [&](int64_t i, int64_t j) {
                        return rowtmp[i] < rowtmp[j];
                      });
            for (int64_t i = 0; i < nelems; ++i) {
              row[start + i] = rowtmp[indices[i]];
              data[start + i] = datatmp[indices[i]];
            }
          }
        }
      });

  bool isherm = ishermitian(ops);
  return CSCMatrix<idx_t, coeff_t>{nrows, ncols, colptr, row, data, i0, isherm};
}
XDIAG_CATCH

// Layer 1: unwrap Block variant, then call the block-specific Layer 2 function.
template <typename idx_t, typename coeff_t>
CSCMatrix<idx_t, coeff_t> csc_matrix(OpSum const &ops, Block const &block_in,
                                     Block const &block_out, idx_t i0) try {
  CSCMatrix<idx_t, coeff_t> result;
  utils::visit_same_type(
      block_in, block_out,
      [&](auto const &bin, auto const &bout) {
        result = csc_matrix<idx_t, coeff_t>(ops, bin, bout, i0);
      },
      "Type mismatch of Block types");
  return result;
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
CSCMatrix<idx_t, coeff_t> csc_matrix(OpSum const &ops, Block const &blocki,
                                     idx_t i0) try {
  auto blocko = block(ops, blocki);
  return csc_matrix<idx_t, coeff_t>(ops, blocki, blocko, i0);
}
XDIAG_CATCH

template CSCMatrix<int32_t, double> csc_matrix<int32_t, double>(OpSum const &,
                                                                Block const &,
                                                                int32_t);
template CSCMatrix<int32_t, complex> csc_matrix<int32_t, complex>(OpSum const &,
                                                                  Block const &,
                                                                  int32_t);
template CSCMatrix<int64_t, double> csc_matrix<int64_t, double>(OpSum const &,
                                                                Block const &,
                                                                int64_t);
template CSCMatrix<int64_t, complex> csc_matrix<int64_t, complex>(OpSum const &,
                                                                  Block const &,
                                                                  int64_t);

template CSCMatrix<int32_t, double> csc_matrix<int32_t, double>(OpSum const &,
                                                                Block const &,
                                                                Block const &,
                                                                int32_t);
template CSCMatrix<int32_t, complex>
csc_matrix<int32_t, complex>(OpSum const &, Block const &, Block const &,
                             int32_t);
template CSCMatrix<int64_t, double> csc_matrix<int64_t, double>(OpSum const &,
                                                                Block const &,
                                                                Block const &,
                                                                int64_t);
template CSCMatrix<int64_t, complex>
csc_matrix<int64_t, complex>(OpSum const &, Block const &, Block const &,
                             int64_t);

// Named convenience wrappers.
CSCMatrix<int64_t, double> csc_matrix(OpSum const &ops, Block const &block,
                                      int64_t i0) try {
  return csc_matrix<int64_t, double>(ops, block, i0);
}
XDIAG_CATCH

CSCMatrix<int64_t, double> csc_matrix(OpSum const &ops, Block const &block_in,
                                      Block const &block_out,
                                      int64_t i0) try {
  return csc_matrix<int64_t, double>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

CSCMatrix<int64_t, complex> csc_matrixC(OpSum const &ops, Block const &block,
                                        int64_t i0) try {
  return csc_matrix<int64_t, complex>(ops, block, i0);
}
XDIAG_CATCH

CSCMatrix<int64_t, complex> csc_matrixC(OpSum const &ops, Block const &block_in,
                                        Block const &block_out,
                                        int64_t i0) try {
  return csc_matrix<int64_t, complex>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

CSCMatrix<int32_t, double> csc_matrix_32(OpSum const &ops, Block const &block,
                                         int32_t i0) try {
  return csc_matrix<int32_t, double>(ops, block, i0);
}
XDIAG_CATCH

CSCMatrix<int32_t, double> csc_matrix_32(OpSum const &ops,
                                         Block const &block_in,
                                         Block const &block_out,
                                         int32_t i0) try {
  return csc_matrix<int32_t, double>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

CSCMatrix<int32_t, complex> csc_matrixC_32(OpSum const &ops, Block const &block,
                                           int32_t i0) try {
  return csc_matrix<int32_t, complex>(ops, block, i0);
}
XDIAG_CATCH

CSCMatrix<int32_t, complex> csc_matrixC_32(OpSum const &ops,
                                           Block const &block_in,
                                           Block const &block_out,
                                           int32_t i0) try {
  return csc_matrix<int32_t, complex>(ops, block_in, block_out, i0);
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
arma::Mat<coeff_t> to_dense(CSCMatrix<idx_t, coeff_t> const &csc_mat) try {
  if ((csc_mat.nrows == 0) || (csc_mat.ncols == 0)) {
    return arma::Mat<coeff_t>();
  }
  int64_t nnz = csc_mat.data.size();
  if (csc_mat.colptr.size() != csc_mat.ncols + 1) {
    XDIAG_THROW("Number of colptr entries does not match number of cols (+1)");
  }
  if (csc_mat.row.size() != nnz) {
    XDIAG_THROW("Number of row entries does not match number of data entries");
  }

  arma::Mat<coeff_t> m(csc_mat.nrows, csc_mat.ncols, arma::fill::zeros);

  if (csc_mat.i0 == 0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (idx_t col = 0; col < csc_mat.ncols; ++col) {
      for (idx_t i = csc_mat.colptr[col]; i < csc_mat.colptr[col + 1]; ++i)
        m(csc_mat.row[i], col) += csc_mat.data[i];
    }
  } else if (csc_mat.i0 == 1) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (idx_t col = 0; col < csc_mat.ncols; ++col) {
      for (idx_t i = csc_mat.colptr[col] - 1; i < csc_mat.colptr[col + 1] - 1;
           ++i)
        m(csc_mat.row[i] - 1, col) += csc_mat.data[i];
    }
  } else {
    XDIAG_THROW(fmt::format(
        "Invalid zero index i0. Must be either 0 or 1, but got i0={}",
        csc_mat.i0));
  }
  return m;
}
XDIAG_CATCH

template arma::mat to_dense(CSCMatrix<int32_t, double> const &);
template arma::mat to_dense(CSCMatrix<int64_t, double> const &);
template arma::cx_mat to_dense(CSCMatrix<int32_t, complex> const &);
template arma::cx_mat to_dense(CSCMatrix<int64_t, complex> const &);

} // namespace xdiag
