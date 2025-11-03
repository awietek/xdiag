// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "eigs_lanczos.hpp"

namespace xdiag::julia {

template <class op_t, class block_t>
static void define_eigs_lanczos_block(jlcxx::Module &mod) {
  mod.method("cxx_eigs_lanczos",
             [](op_t const &ops, block_t const &block, int64_t neigvals,
                double precision, int64_t max_iterations, double deflation_tol,
                int64_t random_seed) {
               JULIA_XDIAG_CALL_RETURN(eigs_lanczos(ops, block, neigvals,
                                                    precision, max_iterations,
                                                    deflation_tol, random_seed))
             });
}

template <class op_t>
static void define_eigs_lanczos_state(jlcxx::Module &mod) {
  mod.method("cxx_eigs_lanczos", [](op_t const &ops, State const &state0,
                                    int64_t neigvals, double precision,
                                    int64_t max_iterations,
                                    double deflation_tol) {
    JULIA_XDIAG_CALL_RETURN(eigs_lanczos(ops, state0, neigvals, precision,
                                         max_iterations, deflation_tol));
  });
}

void define_eigs_lanczos(jlcxx::Module &mod) {

  using res_t = EigsLanczosResult;

  mod.add_type<res_t>("cxx_EigsLanczosResult")
      .method("alphas",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.alphas) })
      .method("betas",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.betas) })
      .method(
          "eigenvalues",
          [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.eigenvalues) })
      .method(
          "eigenvectors",
          [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.eigenvectors) })
      .method(
          "niterations",
          [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.niterations) })
      .method("criterion",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.criterion) });

  define_eigs_lanczos_block<OpSum, Spinhalf>(mod);
  define_eigs_lanczos_block<OpSum, tJ>(mod);
  define_eigs_lanczos_block<OpSum, Electron>(mod);

  define_eigs_lanczos_block<CSRMatrix<int64_t, double>, Spinhalf>(mod);
  define_eigs_lanczos_block<CSRMatrix<int64_t, double>, tJ>(mod);
  define_eigs_lanczos_block<CSRMatrix<int64_t, double>, Electron>(mod);

  define_eigs_lanczos_block<CSRMatrix<int64_t, complex>, Spinhalf>(mod);
  define_eigs_lanczos_block<CSRMatrix<int64_t, complex>, tJ>(mod);
  define_eigs_lanczos_block<CSRMatrix<int64_t, complex>, Electron>(mod);

  define_eigs_lanczos_block<CSRMatrix<int32_t, double>, Spinhalf>(mod);
  define_eigs_lanczos_block<CSRMatrix<int32_t, double>, tJ>(mod);
  define_eigs_lanczos_block<CSRMatrix<int32_t, double>, Electron>(mod);

  define_eigs_lanczos_block<CSRMatrix<int32_t, complex>, Spinhalf>(mod);
  define_eigs_lanczos_block<CSRMatrix<int32_t, complex>, tJ>(mod);
  define_eigs_lanczos_block<CSRMatrix<int32_t, complex>, Electron>(mod);

  define_eigs_lanczos_state<OpSum>(mod);
  define_eigs_lanczos_state<CSRMatrix<int64_t, double>>(mod);
  define_eigs_lanczos_state<CSRMatrix<int64_t, complex>>(mod);
  define_eigs_lanczos_state<CSRMatrix<int32_t, double>>(mod);
  define_eigs_lanczos_state<CSRMatrix<int32_t, complex>>(mod);
}

} // namespace xdiag::julia
