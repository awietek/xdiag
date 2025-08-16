// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "coo_matrix.hpp"

#include <xdiag/algebra/sparse/coo_matrix_fill.hpp>
#include <xdiag/algebra/sparse/coo_matrix_nnz.hpp>
#include <xdiag/all.hpp>

namespace xdiag::julia {

template <typename block_t>
static void define_coo_matrix_nnz_thread(jlcxx::Module &mod) {
  mod.method(
      "cxx_coo_matrix_nnz_thread",
      [](OpSum const &ops, block_t const &block_in, block_t const &block_out) {
        JULIA_XDIAG_CALL_RETURN(
            algebra::coo_matrix_nnz_thread<double>(ops, block_in, block_out));
      });

  mod.method(
      "cxx_coo_matrix_nnz_threadC",
      [](OpSum const &ops, block_t const &block_in, block_t const &block_out) {
        JULIA_XDIAG_CALL_RETURN(
            algebra::coo_matrix_nnz_thread<complex>(ops, block_in, block_out));
      });
}

template <typename idx_t, typename coeff_t, typename block_t>
static void define_coo_matrix_fill(jlcxx::Module &mod) {
  mod.method(
      "cxx_coo_matrix_fill",
      [](OpSum const &ops, block_t const &block_in, block_t const &block_out,
         std::vector<int64_t> const &nnz_thread, int64_t nnz, idx_t *row,
         idx_t *col, coeff_t *data, idx_t i0) {
        JULIA_XDIAG_CALL_RETURN(algebra::coo_matrix_fill(
            ops, block_in, block_out, nnz_thread, nnz, row, col, data, i0));
      });
}

template <typename idx_t, typename coeff_t>
static void define_create_coo_matrix(jlcxx::Module &mod) {
  mod.method("cxx_create_coo_matrix",
             [](idx_t nrows, idx_t ncols, int64_t nnz, idx_t *rowptr,
                idx_t *colptr, coeff_t *dataptr, idx_t i0, bool ishermitian) {
               arma::Col<idx_t> row(rowptr, nnz, false, true);
               arma::Col<idx_t> col(colptr, nnz, false, true);
               arma::Col<coeff_t> data(dataptr, nnz, false, true);
               return COOMatrix<idx_t, coeff_t>{nrows, ncols, row,        col,
                                                data,  i0,    ishermitian};
             });
}

static void define_coo_matrix_types(jlcxx::Module &mod) {
  mod.add_type<COOMatrix<int64_t, double>>("cxx_coo_matrix");
  mod.add_type<COOMatrix<int64_t, complex>>("cxx_coo_matrixC");
  mod.add_type<COOMatrix<int32_t, double>>("cxx_coo_matrix_32");
  mod.add_type<COOMatrix<int32_t, complex>>("cxx_coo_matrixC_32");
}

template <typename idx_t, typename coeff_t>
static void define_coo_matrix_to_dense(jlcxx::Module &mod) {
  mod.method("cxx_to_dense", [](COOMatrix<idx_t, coeff_t> const &spmat) {
    JULIA_XDIAG_CALL_RETURN(to_dense(spmat));
  });
}

void define_coo_matrix(jlcxx::Module &mod) {
  define_coo_matrix_nnz_thread<Spinhalf>(mod);
  define_coo_matrix_nnz_thread<tJ>(mod);
  define_coo_matrix_nnz_thread<Electron>(mod);
  define_coo_matrix_fill<int64_t, double, Spinhalf>(mod);
  define_coo_matrix_fill<int64_t, double, tJ>(mod);
  define_coo_matrix_fill<int64_t, double, Electron>(mod);
  define_coo_matrix_fill<int64_t, complex, Spinhalf>(mod);
  define_coo_matrix_fill<int64_t, complex, tJ>(mod);
  define_coo_matrix_fill<int64_t, complex, Electron>(mod);
  define_coo_matrix_fill<int32_t, double, Spinhalf>(mod);
  define_coo_matrix_fill<int32_t, double, tJ>(mod);
  define_coo_matrix_fill<int32_t, double, Electron>(mod);
  define_coo_matrix_fill<int32_t, complex, Spinhalf>(mod);
  define_coo_matrix_fill<int32_t, complex, tJ>(mod);
  define_coo_matrix_fill<int32_t, complex, Electron>(mod);
  define_coo_matrix_types(mod);
  define_create_coo_matrix<int64_t, double>(mod);
  define_create_coo_matrix<int64_t, complex>(mod);
  define_create_coo_matrix<int32_t, double>(mod);
  define_create_coo_matrix<int32_t, complex>(mod);
  define_coo_matrix_to_dense<int64_t, double>(mod);
  define_coo_matrix_to_dense<int64_t, complex>(mod);
  define_coo_matrix_to_dense<int32_t, double>(mod);
  define_coo_matrix_to_dense<int32_t, complex>(mod);
}

} // namespace xdiag::julia
