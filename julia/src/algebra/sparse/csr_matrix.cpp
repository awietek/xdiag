// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "csr_matrix.hpp"

#include <xdiag/algebra/sparse/csr_matrix_fill.hpp>
#include <xdiag/algebra/sparse/csr_matrix_nnz.hpp>
#include <xdiag/all.hpp>

namespace xdiag::julia {

template <typename block_t>
static void define_csr_matrix_nnz(jlcxx::Module &mod) {
  mod.method("cxx_csr_matrix_nnz", [](OpSum const &ops, block_t const &block_in,
                                      block_t const &block_out,
                                      bool transpose) {
    // somehow, JULIA_XDIAG_CALL_RETURN doesn't like two template parameters
    // -> writing out the macro here

    try {
      return algebra::csr_matrix_nnz<int64_t, double>(ops, block_in, block_out,
                                                      transpose);
    } catch (Error const &e) {
      error_trace(e);
      throw(std::runtime_error("Error occurred in XDiag C++ core library"));
    }
  });
  mod.method(
      "cxx_csr_matrixC_nnz", [](OpSum const &ops, block_t const &block_in,
                                block_t const &block_out, bool transpose) {
        try {
          return algebra::csr_matrix_nnz<int64_t, complex>(
              ops, block_in, block_out, transpose);
        } catch (Error const &e) {
          error_trace(e);
          throw(std::runtime_error("Error occurred in XDiag C++ core library"));
        }
      });
  mod.method(
      "cxx_csr_matrix_32_nnz", [](OpSum const &ops, block_t const &block_in,
                                  block_t const &block_out, bool transpose) {
        try {
          return algebra::csr_matrix_nnz<int32_t, double>(ops, block_in,
                                                          block_out, transpose);
        } catch (Error const &e) {
          error_trace(e);
          throw(std::runtime_error("Error occurred in XDiag C++ core library"));
        }
      });
  mod.method(
      "cxx_csr_matrixC_32_nnz", [](OpSum const &ops, block_t const &block_in,
                                   block_t const &block_out, bool transpose) {
        try {
          return algebra::csr_matrix_nnz<int32_t, complex>(
              ops, block_in, block_out, transpose);
        } catch (Error const &e) {
          error_trace(e);
          throw(std::runtime_error("Error occurred in XDiag C++ core library"));
        }
      });
}

template <typename idx_t, typename coeff_t, typename block_t>
static void define_csr_matrix_fill(jlcxx::Module &mod) {
  mod.method("cxx_csr_matrix_fill",
             [](OpSum const &ops, block_t const &block_in,
                block_t const &block_out, int64_t nnz,
                std::vector<idx_t> const &n_elements_in_row, idx_t *rowptr,
                idx_t *col, coeff_t *data, idx_t i0, bool transpose) {
               JULIA_XDIAG_CALL_VOID(algebra::csr_matrix_fill(
                   ops, block_in, block_out, nnz, n_elements_in_row, rowptr,
                   col, data, i0, transpose));
             });
}

template <typename idx_t, typename coeff_t>
static void define_create_csr_matrix(jlcxx::Module &mod) {
  mod.method("cxx_create_csr_matrix",
             [](idx_t nrows, idx_t ncols, int64_t nnz, idx_t *rowptrptr,
                idx_t *colptr, coeff_t *dataptr, idx_t i0, bool ishermitian) {
               arma::Col<idx_t> rowptr(rowptrptr, nrows + 1, false, true);
               arma::Col<idx_t> col(colptr, nnz, false, true);
               arma::Col<coeff_t> data(dataptr, nnz, false, true);
               return CSRMatrix<idx_t, coeff_t>{nrows, ncols, rowptr,     col,
                                                data,  i0,    ishermitian};
             });
}

static void define_csr_matrix_types(jlcxx::Module &mod) {
  mod.add_type<CSRMatrix<int64_t, double>>("cxx_csr_matrix");
  mod.add_type<CSRMatrix<int64_t, complex>>("cxx_csr_matrixC");
  mod.add_type<CSRMatrix<int32_t, double>>("cxx_csr_matrix_32");
  mod.add_type<CSRMatrix<int32_t, complex>>("cxx_csr_matrixC_32");
}

template <typename idx_t, typename coeff_t>
static void define_csr_matrix_to_dense(jlcxx::Module &mod) {
  mod.method("cxx_to_dense", [](CSRMatrix<idx_t, coeff_t> const &spmat) {
    JULIA_XDIAG_CALL_RETURN(to_dense(spmat));
  });
}

void define_csr_matrix(jlcxx::Module &mod) {
  define_csr_matrix_nnz<Spinhalf>(mod);
  define_csr_matrix_nnz<tJ>(mod);
  define_csr_matrix_nnz<Electron>(mod);
  define_csr_matrix_fill<int64_t, double, Spinhalf>(mod);
  define_csr_matrix_fill<int64_t, double, tJ>(mod);
  define_csr_matrix_fill<int64_t, double, Electron>(mod);
  define_csr_matrix_fill<int64_t, complex, Spinhalf>(mod);
  define_csr_matrix_fill<int64_t, complex, tJ>(mod);
  define_csr_matrix_fill<int64_t, complex, Electron>(mod);
  define_csr_matrix_fill<int32_t, double, Spinhalf>(mod);
  define_csr_matrix_fill<int32_t, double, tJ>(mod);
  define_csr_matrix_fill<int32_t, double, Electron>(mod);
  define_csr_matrix_fill<int32_t, complex, Spinhalf>(mod);
  define_csr_matrix_fill<int32_t, complex, tJ>(mod);
  define_csr_matrix_fill<int32_t, complex, Electron>(mod);
  define_csr_matrix_types(mod);
  define_create_csr_matrix<int64_t, double>(mod);
  define_create_csr_matrix<int64_t, complex>(mod);
  define_create_csr_matrix<int32_t, double>(mod);
  define_create_csr_matrix<int32_t, complex>(mod);
  define_csr_matrix_to_dense<int64_t, double>(mod);
  define_csr_matrix_to_dense<int64_t, complex>(mod);
  define_csr_matrix_to_dense<int32_t, double>(mod);
  define_csr_matrix_to_dense<int32_t, complex>(mod);
}

} // namespace xdiag::julia
