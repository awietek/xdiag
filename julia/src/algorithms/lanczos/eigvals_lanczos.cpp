// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "eigvals_lanczos.hpp"

namespace xdiag::julia {

template <class op_t, class block_t>
static void define_eigvals_lanczos_block(jlcxx::Module &mod) {
  mod.method("cxx_eigvals_lanczos",
             [](op_t const &ops, block_t const &block, int64_t neigvals,
                double precision, int64_t max_iterations, double deflation_tol,
                int64_t random_seed) {
               JULIA_XDIAG_CALL_RETURN(
                   eigvals_lanczos(ops, block, neigvals, precision,
                                   max_iterations, deflation_tol, random_seed))
             });
}

template <class op_t>
static void define_eigvals_lanczos_state(jlcxx::Module &mod) {
  mod.method("cxx_eigvals_lanczos", [](op_t const &ops, State psi0,
                                       int64_t neigvals, double precision,
                                       int64_t max_iterations,
                                       double deflation_tol) {
    JULIA_XDIAG_CALL_RETURN(eigvals_lanczos(ops, psi0, neigvals, precision,
                                            max_iterations, deflation_tol))
  });

  mod.method(
      "cxx_eigvals_lanczos_inplace",
      [](op_t const &ops, State &psi0, int64_t neigvals, double precision,
         int64_t max_iterations, double deflation_tol) {
        JULIA_XDIAG_CALL_RETURN(eigvals_lanczos_inplace(
            ops, psi0, neigvals, precision, max_iterations, deflation_tol))
      });
}

void define_eigvals_lanczos(jlcxx::Module &mod) {

  using res_t = EigvalsLanczosResult;

  mod.add_type<res_t>("cxx_EigvalsLanczosResult")
      .method("alphas",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.alphas) })
      .method("betas",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.betas) })
      .method(
          "eigenvalues",
          [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.eigenvalues) })
      .method(
          "niterations",
          [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.niterations) })
      .method("criterion",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.criterion) });

  define_eigvals_lanczos_block<OpSum, Spinhalf>(mod);
  define_eigvals_lanczos_block<OpSum, tJ>(mod);
  define_eigvals_lanczos_block<OpSum, Electron>(mod);

  define_eigvals_lanczos_block<CSRMatrix<int64_t, double>, Spinhalf>(mod);
  define_eigvals_lanczos_block<CSRMatrix<int64_t, double>, tJ>(mod);
  define_eigvals_lanczos_block<CSRMatrix<int64_t, double>, Electron>(mod);

  define_eigvals_lanczos_block<CSRMatrix<int64_t, complex>, Spinhalf>(mod);
  define_eigvals_lanczos_block<CSRMatrix<int64_t, complex>, tJ>(mod);
  define_eigvals_lanczos_block<CSRMatrix<int64_t, complex>, Electron>(mod);

  define_eigvals_lanczos_block<CSRMatrix<int32_t, double>, Spinhalf>(mod);
  define_eigvals_lanczos_block<CSRMatrix<int32_t, double>, tJ>(mod);
  define_eigvals_lanczos_block<CSRMatrix<int32_t, double>, Electron>(mod);

  define_eigvals_lanczos_block<CSRMatrix<int32_t, complex>, Spinhalf>(mod);
  define_eigvals_lanczos_block<CSRMatrix<int32_t, complex>, tJ>(mod);
  define_eigvals_lanczos_block<CSRMatrix<int32_t, complex>, Electron>(mod);

  define_eigvals_lanczos_state<OpSum>(mod);
  define_eigvals_lanczos_state<CSRMatrix<int64_t, double>>(mod);
  define_eigvals_lanczos_state<CSRMatrix<int64_t, complex>>(mod);
  define_eigvals_lanczos_state<CSRMatrix<int32_t, double>>(mod);
  define_eigvals_lanczos_state<CSRMatrix<int32_t, complex>>(mod);
}

} // namespace xdiag::julia
