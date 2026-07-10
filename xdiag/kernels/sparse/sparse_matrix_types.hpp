#pragma once

#include <xdiag/armadillo.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

// COO format
template <typename idx_t, typename coeff_t> struct XDIAG_API COOMatrix {
  idx_t nrows;             // number of rows in the matrix
  idx_t ncols;             // number of columns in the matrix
  arma::Col<idx_t> row;    // row indices
  arma::Col<idx_t> col;    // column indices
  arma::Col<coeff_t> data; // data entries
  idx_t i0;                // zero index, either 0 or 1
  bool ishermitian;        // flag whether matrix is hermitian/symmetric
};

// CSR format
template <typename idx_t, typename coeff_t> struct XDIAG_API CSRMatrix {
  idx_t nrows;             // number of rows in the matrix
  idx_t ncols;             // number of columns in the matrix
  arma::Col<idx_t> rowptr; // pointer to elements of a row (size: nrows+1)
  arma::Col<idx_t> col;    // columns for each row consecutively
  arma::Col<coeff_t> data; // data entries each row consecutively
  idx_t i0;                // zero index, either 0 or 1
  bool ishermitian;        // flag whether matrix is hermitian/symmetric
};

// CSC format
template <typename idx_t, typename coeff_t> struct XDIAG_API CSCMatrix {
  idx_t nrows;             // number of rows in the matrix
  idx_t ncols;             // number of columns in the matrix
  arma::Col<idx_t> colptr; // pointer to elements of a row (size: ncols+1)
  arma::Col<idx_t> row;    // columns for each row consecutively
  arma::Col<coeff_t> data; // data entries each row consecutively
  idx_t i0;                // zero index, either 0 or 1
  bool ishermitian;        // flag whether matrix is hermitian/symmetric
};

template <typename idx_t, typename coeff_t>
constexpr XDIAG_API bool isreal(CSRMatrix<idx_t, coeff_t> const &A) {
  (void)A; // need argument name for wrapper generator -> keep it
  return isreal<coeff_t>();
}
template <typename idx_t, typename coeff_t>
constexpr XDIAG_API bool isreal(CSCMatrix<idx_t, coeff_t> const &A) {
  (void)A;
  return isreal<coeff_t>();
}
template <typename idx_t, typename coeff_t>
constexpr XDIAG_API bool isreal(COOMatrix<idx_t, coeff_t> const &A) {
  (void)A;
  return isreal<coeff_t>();
}

template <typename idx_t, typename coeff_t>
constexpr XDIAG_API bool ishermitian(CSRMatrix<idx_t, coeff_t> const &A) {
  return A.ishermitian;
}
template <typename idx_t, typename coeff_t>
constexpr XDIAG_API bool ishermitian(CSCMatrix<idx_t, coeff_t> const &A) {
  return A.ishermitian;
}
template <typename idx_t, typename coeff_t>
constexpr XDIAG_API bool ishermitian(COOMatrix<idx_t, coeff_t> const &A) {
  return A.ishermitian;
}

} // namespace xdiag
