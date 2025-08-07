// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"
#include <xdiag/all.hpp>

namespace xdiag::julia {

template <typename idx_t, typename coeff_t>
static void define_sparse_apply_csr(jlcxx::Module &mod) {

  mod.method("cxx_apply",
             [](CSRMatrix<idx_t, coeff_t> const &spmat,
                arma::Col<coeff_t> const &vec_in) {
               JULIA_XDIAG_CALL_RETURN(apply(spmat, vec_in));
             });
  mod.method("cxx_apply",
             [](CSRMatrix<idx_t, coeff_t> const &spmat,
                arma::Mat<coeff_t> const &mat_in) {
               JULIA_XDIAG_CALL_RETURN(apply(spmat, mat_in));
             });
  
  mod.method("cxx_apply",
             [](CSRMatrix<idx_t, coeff_t> const &spmat,
                arma::Col<coeff_t> const &vec_in, arma::Col<coeff_t> &vec_out) {
               JULIA_XDIAG_CALL_VOID(apply(spmat, vec_in, vec_out));
             });
  mod.method("cxx_apply",
             [](CSRMatrix<idx_t, coeff_t> const &spmat,
                arma::Mat<coeff_t> const &mat_in, arma::Mat<coeff_t> &mat_out) {
               JULIA_XDIAG_CALL_VOID(apply(spmat, mat_in, mat_out));
             });
}

template <typename idx_t>
static void define_sparse_apply_csr_mixed(jlcxx::Module &mod) {
  mod.method("cxx_apply",
             [](CSRMatrix<idx_t, double> const &spmat,
                arma::Col<complex> const &vec_in) {
               JULIA_XDIAG_CALL_RETURN(apply(spmat, vec_in));
             });
  mod.method("cxx_apply",
             [](CSRMatrix<idx_t, double> const &spmat,
                arma::Mat<complex> const &mat_in) {
               JULIA_XDIAG_CALL_RETURN(apply(spmat, mat_in));
             });
  
  mod.method("cxx_apply",
             [](CSRMatrix<idx_t, double> const &spmat,
                arma::Col<complex> const &vec_in, arma::Col<complex> &vec_out) {
               JULIA_XDIAG_CALL_RETURN(apply(spmat, vec_in, vec_out));
             });
  mod.method("cxx_apply",
             [](CSRMatrix<idx_t, double> const &spmat,
                arma::Mat<complex> const &mat_in, arma::Mat<complex> &mat_out) {
               JULIA_XDIAG_CALL_RETURN(apply(spmat, mat_in, mat_out));
             });
}

void define_sparse_apply(jlcxx::Module &mod) {
  define_sparse_apply_csr<int64_t, double>(mod);
  define_sparse_apply_csr<int64_t, complex>(mod);
  define_sparse_apply_csr<int32_t, double>(mod);
  define_sparse_apply_csr<int32_t, complex>(mod);
  define_sparse_apply_csr_mixed<int64_t>(mod);
  define_sparse_apply_csr_mixed<int32_t>(mod);
}

} // namespace xdiag::julia
