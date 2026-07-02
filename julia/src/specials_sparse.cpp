// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0
//
// Hand-written specials: the C++ CSRMatrix object machinery. The Julia side
// stores sparse matrices as plain Julia structs (built via the two-phase
// nnz/fill); to run apply / the CSR-input algorithms it reconstructs a C++
// CSRMatrix here (cxx_make_csr_*), keeps it alive, and calls the (explicitly
// instantiated) template overloads. Registered as cxx_ names; the ergonomic
// Julia layer (ergonomics/sparse.jl) provides the idiomatic surface.

#include <julia/src/xdiagjl.hpp>

#include <xdiag/kernels/sparse/apply.hpp>
#include <xdiag/kernels/sparse/csr_matrix.hpp>
#include <xdiag/kernels/sparse/sparse_matrix_types.hpp>
#include <xdiag/linalg/lanczos/eigs_lanczos.hpp>
#include <xdiag/linalg/lanczos/eigvals_lanczos.hpp>
#include <xdiag/linalg/sparse_diag.hpp>
#include <xdiag/linalg/time_evolution/evolve_lanczos.hpp>
#include <xdiag/linalg/time_evolution/imaginary_time_evolve.hpp>
#include <xdiag/linalg/time_evolution/time_evolve.hpp>
#include <xdiag/linalg/time_evolution/time_evolve_expokit.hpp>

namespace xdiag::julia {

// --- CSR object: construction + apply -----------------------------------
template <typename idx_t, typename coeff_t>
static void define_csr_core(jlcxx::Module &mod, std::string const &sfx) {
  using CSR = CSRMatrix<idx_t, coeff_t>;
  mod.add_type<CSR>("cxx_CSRMatrix_" + sfx);

  // Build a C++ CSRMatrix from Julia-owned arrays (copied into arma storage).
  mod.method("cxx_make_csr_" + sfx,
             [](idx_t nrows, idx_t ncols, idx_t *rowptr, int64_t nrowptr,
                idx_t *col, int64_t ncol, coeff_t *data, int64_t ndata,
                idx_t i0, bool isherm) {
               CSR A;
               A.nrows = nrows;
               A.ncols = ncols;
               A.rowptr = arma::Col<idx_t>(rowptr, nrowptr);
               A.col = arma::Col<idx_t>(col, ncol);
               A.data = arma::Col<coeff_t>(data, ndata);
               A.i0 = i0;
               A.ishermitian = isherm;
               return A;
             });

  mod.method("cxx_to_dense", [](CSR const &A) {
    JULIA_XDIAG_CALL_RETURN(to_dense(A));
  });

  mod.method("cxx_apply", [](CSR const &A, arma::Col<coeff_t> const &v) {
    JULIA_XDIAG_CALL_RETURN(apply(A, v));
  });
  mod.method("cxx_apply", [](CSR const &A, arma::Col<coeff_t> const &v,
                             arma::Col<coeff_t> &w) {
    JULIA_XDIAG_CALL_VOID(apply(A, v, w));
  });
  mod.method("cxx_apply", [](CSR const &A, arma::Mat<coeff_t> const &v) {
    JULIA_XDIAG_CALL_RETURN(apply(A, v));
  });
  mod.method("cxx_apply", [](CSR const &A, arma::Mat<coeff_t> const &v,
                             arma::Mat<coeff_t> &w) {
    JULIA_XDIAG_CALL_VOID(apply(A, v, w));
  });
}

// Real matrix applied to a complex vector (extra instantiated overloads).
template <typename idx_t>
static void define_csr_apply_cx(jlcxx::Module &mod) {
  using CSR = CSRMatrix<idx_t, double>;
  mod.method("cxx_apply", [](CSR const &A, arma::cx_vec const &v) {
    JULIA_XDIAG_CALL_RETURN(apply(A, v));
  });
  mod.method("cxx_apply", [](CSR const &A, arma::cx_vec const &v,
                             arma::cx_vec &w) {
    JULIA_XDIAG_CALL_VOID(apply(A, v, w));
  });
  mod.method("cxx_apply", [](CSR const &A, arma::cx_mat const &v) {
    JULIA_XDIAG_CALL_RETURN(apply(A, v));
  });
  mod.method("cxx_apply", [](CSR const &A, arma::cx_mat const &v,
                             arma::cx_mat &w) {
    JULIA_XDIAG_CALL_VOID(apply(A, v, w));
  });
}

// --- CSR-input algorithms with a Block argument -------------------------
template <typename idx_t, typename coeff_t, typename block_t>
static void define_csr_block_algos(jlcxx::Module &mod) {
  using CSR = CSRMatrix<idx_t, coeff_t>;
  mod.method("cxx_eigval0", [](CSR const &A, block_t const &b, double prec,
                               int64_t maxit, int64_t seed) {
    JULIA_XDIAG_CALL_RETURN(eigval0(A, b, prec, maxit, seed));
  });
  mod.method("cxx_eig0", [](CSR const &A, block_t const &b, double prec,
                            int64_t maxit, int64_t seed) {
    JULIA_XDIAG_CALL_RETURN(eig0(A, b, prec, maxit, seed));
  });
  mod.method("cxx_eigs", [](CSR const &A, block_t const &b, int64_t neigs,
                            double prec, int64_t maxit, int64_t seed) {
    JULIA_XDIAG_CALL_RETURN(eigs(A, b, neigs, prec, maxit, seed));
  });
  mod.method("cxx_eigs_lanczos",
             [](CSR const &A, block_t const &b, int64_t neig, double prec,
                int64_t maxit, double defl, int64_t seed) {
               JULIA_XDIAG_CALL_RETURN(
                   eigs_lanczos(A, b, neig, prec, maxit, defl, seed));
             });
  mod.method("cxx_eigvals_lanczos",
             [](CSR const &A, block_t const &b, int64_t neig, double prec,
                int64_t maxit, double defl, int64_t seed) {
               JULIA_XDIAG_CALL_RETURN(
                   eigvals_lanczos(A, b, neig, prec, maxit, defl, seed));
             });
}

// --- CSR-input algorithms with a State argument -------------------------
template <typename idx_t, typename coeff_t>
static void define_csr_state_algos(jlcxx::Module &mod) {
  using CSR = CSRMatrix<idx_t, coeff_t>;
  mod.method("cxx_eigs_lanczos",
             [](CSR const &A, State const &s, int64_t neig, double prec,
                int64_t maxit, double defl) {
               JULIA_XDIAG_CALL_RETURN(
                   eigs_lanczos(A, s, neig, prec, maxit, defl));
             });
  mod.method("cxx_eigvals_lanczos",
             [](CSR const &A, State const &s, int64_t neig, double prec,
                int64_t maxit, double defl) {
               JULIA_XDIAG_CALL_RETURN(
                   eigvals_lanczos(A, s, neig, prec, maxit, defl));
             });
  mod.method("cxx_eigvals_lanczos_inplace",
             [](CSR const &A, State &s, int64_t neig, double prec,
                int64_t maxit, double defl) {
               JULIA_XDIAG_CALL_RETURN(
                   eigvals_lanczos_inplace(A, s, neig, prec, maxit, defl));
             });
  mod.method("cxx_evolve_lanczos_inplace",
             [](CSR const &A, State &s, double t, double prec, double shift,
                bool normalize, int64_t maxit, double defl) {
               JULIA_XDIAG_CALL_RETURN(evolve_lanczos_inplace(
                   A, s, t, prec, shift, normalize, maxit, defl));
             });
  mod.method("cxx_evolve_lanczos_inplace",
             [](CSR const &A, State &s, complex t, double prec, double shift,
                bool normalize, int64_t maxit, double defl) {
               JULIA_XDIAG_CALL_RETURN(evolve_lanczos_inplace(
                   A, s, t, prec, shift, normalize, maxit, defl));
             });
  mod.method("cxx_time_evolve_expokit_inplace",
             [](CSR const &A, State &s, double t, double prec, int64_t m,
                double anorm, int64_t nnorm) {
               JULIA_XDIAG_CALL_RETURN(
                   time_evolve_expokit_inplace(A, s, t, prec, m, anorm, nnorm));
             });
  mod.method("cxx_evolve_lanczos",
             [](CSR const &A, State const &s, double t, double prec,
                double shift, bool normalize, int64_t maxit, double defl) {
               JULIA_XDIAG_CALL_RETURN(
                   evolve_lanczos(A, s, t, prec, shift, normalize, maxit, defl));
             });
  mod.method("cxx_evolve_lanczos",
             [](CSR const &A, State const &s, complex t, double prec,
                double shift, bool normalize, int64_t maxit, double defl) {
               JULIA_XDIAG_CALL_RETURN(
                   evolve_lanczos(A, s, t, prec, shift, normalize, maxit, defl));
             });
  mod.method("cxx_time_evolve", [](CSR const &A, State const &s, double t,
                                   double prec, std::string const &algo) {
    JULIA_XDIAG_CALL_RETURN(time_evolve(A, s, t, prec, algo));
  });
  mod.method("cxx_time_evolve_inplace",
             [](CSR const &A, State &s, double t, double prec,
                std::string const &algo) {
               JULIA_XDIAG_CALL_VOID(time_evolve_inplace(A, s, t, prec, algo));
             });
  mod.method("cxx_imaginary_time_evolve",
             [](CSR const &A, State const &s, double t, double prec,
                double shift) {
               JULIA_XDIAG_CALL_RETURN(
                   imaginary_time_evolve(A, s, t, prec, shift));
             });
  mod.method("cxx_imaginary_time_evolve_inplace",
             [](CSR const &A, State &s, double t, double prec, double shift) {
               JULIA_XDIAG_CALL_VOID(
                   imaginary_time_evolve_inplace(A, s, t, prec, shift));
             });
  mod.method("cxx_time_evolve_expokit",
             [](CSR const &A, State const &s, double t, double prec, int64_t m,
                double anorm, int64_t nnorm) {
               JULIA_XDIAG_CALL_RETURN(
                   time_evolve_expokit(A, s, t, prec, m, anorm, nnorm));
             });
}

template <typename idx_t, typename coeff_t>
static void define_csr_all(jlcxx::Module &mod, std::string const &sfx) {
  define_csr_core<idx_t, coeff_t>(mod, sfx);
  define_csr_block_algos<idx_t, coeff_t, Spinhalf>(mod);
  define_csr_block_algos<idx_t, coeff_t, tJ>(mod);
  define_csr_block_algos<idx_t, coeff_t, Electron>(mod);
  define_csr_block_algos<idx_t, coeff_t, Boson>(mod);
  define_csr_block_algos<idx_t, coeff_t, Fermion>(mod);
  define_csr_state_algos<idx_t, coeff_t>(mod);
}

void define_specials_sparse(jlcxx::Module &mod) {
  define_csr_all<int64_t, double>(mod, "i64f64");
  define_csr_all<int32_t, double>(mod, "i32f64");
  define_csr_all<int64_t, complex>(mod, "i64cx");
  define_csr_all<int32_t, complex>(mod, "i32cx");
  define_csr_apply_cx<int64_t>(mod);
  define_csr_apply_cx<int32_t>(mod);
}

} // namespace xdiag::julia
