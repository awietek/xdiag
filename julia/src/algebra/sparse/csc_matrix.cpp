// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "csc_matrix.hpp"
#include <xdiag/all.hpp>

namespace xdiag::julia {

template <typename idx_t, typename coeff_t>
static void define_create_csc_matrix(jlcxx::Module &mod) {
  mod.method("cxx_create_csc_matrix",
             [](idx_t nrows, idx_t ncols, int64_t nnz, idx_t *colptrptr,
                idx_t *rowptr, coeff_t *dataptr, idx_t i0, bool ishermitian) {
               arma::Col<idx_t> colptr(colptrptr, ncols + 1, false, true);
               arma::Col<idx_t> row(rowptr, nnz, false, true);
               arma::Col<coeff_t> data(dataptr, nnz, false, true);
               return CSCMatrix<idx_t, coeff_t>{nrows, ncols, colptr,     row,
                                                data,  i0,    ishermitian};
             });
}

static void define_csc_matrix_types(jlcxx::Module &mod) {
  mod.add_type<CSCMatrix<int64_t, double>>("cxx_csc_matrix");
  mod.add_type<CSCMatrix<int64_t, complex>>("cxx_csc_matrixC");
  mod.add_type<CSCMatrix<int32_t, double>>("cxx_csc_matrix_32");
  mod.add_type<CSCMatrix<int32_t, complex>>("cxx_csc_matrixC_32");
}

template <typename idx_t, typename coeff_t>
static void define_csc_matrix_to_dense(jlcxx::Module &mod) {
  mod.method("cxx_to_dense", [](CSCMatrix<idx_t, coeff_t> const &spmat) {
    JULIA_XDIAG_CALL_RETURN(to_dense(spmat));
  });
}

void define_csc_matrix(jlcxx::Module &mod) {
  define_csc_matrix_types(mod);
  define_create_csc_matrix<int64_t, double>(mod);
  define_create_csc_matrix<int64_t, complex>(mod);
  define_create_csc_matrix<int32_t, double>(mod);
  define_create_csc_matrix<int32_t, complex>(mod);
  define_csc_matrix_to_dense<int64_t, double>(mod);
  define_csc_matrix_to_dense<int64_t, complex>(mod);
  define_csc_matrix_to_dense<int32_t, double>(mod);
  define_csc_matrix_to_dense<int32_t, complex>(mod);
}

} // namespace xdiag::julia
