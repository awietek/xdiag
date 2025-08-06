#pragma once

#include <vector>

#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>

namespace xdiag {

// COO format
template <typename idx_t, typename coeff_t> struct COOMatrix {
  idx_t nrows;             // number of rows in the matrix
  idx_t ncols;             // number of columns in the matrix
  arma::Col<idx_t> row;    // row indices
  arma::Col<idx_t> col;    // column indices
  arma::Col<coeff_t> data; // data entries
  idx_t i0;                // zero index, either 0 or 1
  bool ishermitian;        // flag whether matrix is hermitian/symmetric
};

// CSR format
template <typename idx_t, typename coeff_t> struct CSRMatrix {
  idx_t nrows;             // number of rows in the matrix
  idx_t ncols;             // number of columns in the matrix
  arma::Col<idx_t> rowptr; // pointer to elements of a row (size: nrows+1)
  arma::Col<idx_t> col;    // columns for each row consecutively
  arma::Col<coeff_t> data; // data entries each row consecutively
  idx_t i0;                // zero index, either 0 or 1
  bool ishermitian;        // flag whether matrix is hermitian/symmetric
};

// CSC format
template <typename idx_t, typename coeff_t> struct CSCMatrix {
  idx_t nrows;             // number of rows in the matrix
  idx_t ncols;             // number of columns in the matrix
  arma::Col<idx_t> colptr; // pointer to elements of a row (size: ncols+1)
  arma::Col<idx_t> row;    // columns for each row consecutively
  arma::Col<coeff_t> data; // data entries each row consecutively
  idx_t i0;                // zero index, either 0 or 1
  bool ishermitian;        // flag whether matrix is hermitian/symmetric
};

} // namespace xdiag
