#pragma once

#include <vector>

#include <xdiag/common.hpp>

namespace xdiag {

// COO format
template <typename idx_t, typename coeff_t> struct COOMatrix {
  idx_t nrows;               // number of rows in the matrix
  idx_t ncols;               // number of columns in the matrix
  std::vector<idx_t> row;    // row indices
  std::vector<idx_t> col;    // column indices
  std::vector<coeff_t> data; // data entries
  idx_t i0;                // zero index, either 0 or 1
};

using coo_mat = COOMatrix<int64_t, double>;
using coo_cx_mat = COOMatrix<int64_t, complex>;

// CSR format
template <typename idx_t, typename coeff_t> struct CSRMatrix {
  idx_t nrows;               // number of rows in the matrix
  idx_t ncols;               // number of columns in the matrix
  std::vector<idx_t> rowptr; // pointer to elements of a row (size: nrows+1)
  std::vector<idx_t> col;    // columns for each row consecutively
  std::vector<coeff_t> data; // data entries each row consecutively
  idx_t i0;                  // zero index, either 0 or 1
};

using csr_mat = CSRMatrix<int64_t, double>;
using csr_cx_mat = CSRMatrix<int64_t, complex>;

// CSC format
template <typename idx_t, typename coeff_t> struct CSCMatrix {
  idx_t nrows;               // number of rows in the matrix
  idx_t ncols;               // number of columns in the matrix
  std::vector<idx_t> colptr; // pointer to elements of a row (size: nrows+1)
  std::vector<idx_t> row;    // columns for each row consecutively
  std::vector<coeff_t> data; // data entries each row consecutively
  idx_t i0;                  // zero index, either 0 or 1
};

using csc_mat = CSCMatrix<int64_t, double>;
using csc_cx_mat = CSCMatrix<int64_t, complex>;


  
} // namespace xdiag
