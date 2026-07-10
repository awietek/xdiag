// SPDX-License-Identifier: Apache-2.0
// STATIC hand-written special -- copied verbatim by generate.sh.
//
// The CSRMatrix<idx,coeff> handle + a non-owning view constructor over the
// Julia-owned arrays. The CSRMatrix-accepting algorithms (eigval0, eig0,
// eigs_lanczos, ...) are wrapped automatically by the generator; they call
// fun_csr_view (via the Julia csr_view helper) and dispatch on these types.

#include <julia/src/sparse_algebra.hpp> // IsMirroredType specializations first
#include <julia/src/xdiagjl.hpp>
#include <xdiag/kernels/sparse/apply.hpp>

namespace xdiag::julia {

template <typename idx_t, typename coeff_t>
static void register_csr(jlcxx::Module &mod, const char *name) {
  mod.add_type<CSRMatrix<idx_t, coeff_t>>(name);
  // Non-owning view: the arma columns reference the caller's (Julia) memory.
  // Constructed in place (aggregate init) so the members are views, not copies.
  mod.method("fun_csr_view", [](int64_t nrows, int64_t ncols, idx_t *rowptr,
                                idx_t *col, coeff_t *data, int64_t nnz,
                                idx_t i0, bool ishermitian) {
    return CSRMatrix<idx_t, coeff_t>{
        (idx_t)nrows,
        (idx_t)ncols,
        arma::Col<idx_t>(rowptr, (arma::uword)(nrows + 1), false, true),
        arma::Col<idx_t>(col, (arma::uword)nnz, false, true),
        arma::Col<coeff_t>(data, (arma::uword)nnz, false, true),
        i0,
        ishermitian};
  });
  mod.method("fun_csr_to_dense", [](CSRMatrix<idx_t, coeff_t> const &A) {
    JULIA_XDIAG_CALL_RETURN(to_dense(A));
  });
  // In-place matvec / matmat y = A*x over Julia-owned x/y (same coeff type).
  mod.method("fun_csr_apply_vec",
             [](CSRMatrix<idx_t, coeff_t> const &A, coeff_t *xin, coeff_t *yout,
                int64_t nin, int64_t nout) {
               arma::Col<coeff_t> x(xin, (arma::uword)nin, false, true);
               arma::Col<coeff_t> y(yout, (arma::uword)nout, false, true);
               JULIA_XDIAG_CALL_VOID(apply(A, x, y));
             });
  mod.method("fun_csr_apply_mat",
             [](CSRMatrix<idx_t, coeff_t> const &A, coeff_t *xin, coeff_t *yout,
                int64_t nrin, int64_t ncols, int64_t nrout) {
               arma::Mat<coeff_t> x(xin, (arma::uword)nrin, (arma::uword)ncols,
                                    false, true);
               arma::Mat<coeff_t> y(yout, (arma::uword)nrout,
                                    (arma::uword)ncols, false, true);
               JULIA_XDIAG_CALL_VOID(apply(A, x, y));
             });
}

// Mixed precision: apply a real CSRMatrix to a complex vector/matrix.
template <typename idx_t>
static void register_csr_apply_mixed(jlcxx::Module &mod) {
  mod.method("fun_csr_apply_vec",
             [](CSRMatrix<idx_t, double> const &A, complex *xin, complex *yout,
                int64_t nin, int64_t nout) {
               arma::Col<complex> x(xin, (arma::uword)nin, false, true);
               arma::Col<complex> y(yout, (arma::uword)nout, false, true);
               JULIA_XDIAG_CALL_VOID(apply(A, x, y));
             });
  mod.method("fun_csr_apply_mat",
             [](CSRMatrix<idx_t, double> const &A, complex *xin, complex *yout,
                int64_t nrin, int64_t ncols, int64_t nrout) {
               arma::Mat<complex> x(xin, (arma::uword)nrin, (arma::uword)ncols,
                                    false, true);
               arma::Mat<complex> y(yout, (arma::uword)nrout,
                                    (arma::uword)ncols, false, true);
               JULIA_XDIAG_CALL_VOID(apply(A, x, y));
             });
}

// COO/CSC only need to_dense here (no algorithms), so fuse view + to_dense into
// one call and avoid add_type-ing those handles.
template <typename idx_t, typename coeff_t>
static void register_dense(jlcxx::Module &mod) {
  mod.method("fun_coo_to_dense",
             [](int64_t nrows, int64_t ncols, idx_t *row, idx_t *col,
                coeff_t *data, int64_t nnz, idx_t i0) {
               COOMatrix<idx_t, coeff_t> A{
                   (idx_t)nrows,
                   (idx_t)ncols,
                   arma::Col<idx_t>(row, (arma::uword)nnz, false, true),
                   arma::Col<idx_t>(col, (arma::uword)nnz, false, true),
                   arma::Col<coeff_t>(data, (arma::uword)nnz, false, true),
                   i0,
                   false};
               JULIA_XDIAG_CALL_RETURN(to_dense(A));
             });
  mod.method("fun_csc_to_dense", [](int64_t nrows, int64_t ncols, idx_t *colptr,
                                    idx_t *row, coeff_t *data, int64_t nnz,
                                    idx_t i0) {
    CSCMatrix<idx_t, coeff_t> A{
        (idx_t)nrows,
        (idx_t)ncols,
        arma::Col<idx_t>(colptr, (arma::uword)(ncols + 1), false, true),
        arma::Col<idx_t>(row, (arma::uword)nnz, false, true),
        arma::Col<coeff_t>(data, (arma::uword)nnz, false, true),
        i0,
        false};
    JULIA_XDIAG_CALL_RETURN(to_dense(A));
  });
}

void define_sparse_algebra(jlcxx::Module &mod) {
  using namespace xdiag;
  register_csr<int64_t, double>(mod, "typ_csr_i64_f64");
  register_csr<int64_t, complex>(mod, "typ_csr_i64_cx");
  register_csr<int32_t, double>(mod, "typ_csr_i32_f64");
  register_csr<int32_t, complex>(mod, "typ_csr_i32_cx");
  register_dense<int64_t, double>(mod);
  register_dense<int64_t, complex>(mod);
  register_dense<int32_t, double>(mod);
  register_dense<int32_t, complex>(mod);
  register_csr_apply_mixed<int64_t>(mod);
  register_csr_apply_mixed<int32_t>(mod);
}

} // namespace xdiag::julia
