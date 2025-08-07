// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "evolve_lanczos.hpp"

namespace xdiag::julia {

template <typename op_t>
static void define_evolve_lanczos_op(jlcxx::Module &mod) {

  mod.method("cxx_evolve_lanczos", [](op_t const &H, State psi, double tau,
                                      double precision, double shift,
                                      bool normalize, int64_t max_iterations,
                                      double deflation_tol) {
    JULIA_XDIAG_CALL_RETURN(evolve_lanczos(H, psi, tau, precision, shift,
                                           normalize, max_iterations,
                                           deflation_tol));
  });
  mod.method("cxx_evolve_lanczos", [](op_t const &H, State psi, complex tau,
                                      double precision, double shift,
                                      bool normalize, int64_t max_iterations,
                                      double deflation_tol) {
    JULIA_XDIAG_CALL_RETURN(evolve_lanczos(H, psi, tau, precision, shift,
                                           normalize, max_iterations,
                                           deflation_tol));
  });

  mod.method(
      "cxx_evolve_lanczos_inplace",
      [](op_t const &H, State &psi, double tau, double precision, double shift,
         bool normalize, int64_t max_iterations, double deflation_tol) {
        JULIA_XDIAG_CALL_RETURN(
            evolve_lanczos_inplace(H, psi, tau, precision, shift, normalize,
                                   max_iterations, deflation_tol));
      });
  mod.method(
      "cxx_evolve_lanczos_inplace",
      [](op_t const &H, State &psi, complex tau, double precision, double shift,
         bool normalize, int64_t max_iterations, double deflation_tol) {
        JULIA_XDIAG_CALL_RETURN(
            evolve_lanczos_inplace(H, psi, tau, precision, shift, normalize,
                                   max_iterations, deflation_tol));
      });
}

void define_evolve_lanczos(jlcxx::Module &mod) {

  using res_t = EvolveLanczosResult;
  mod.add_type<res_t>("cxx_EvolveLanczosResult")
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
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.criterion) })
      .method("state",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.state) });

  using res_inplace_t = EvolveLanczosInplaceResult;
  mod.add_type<res_inplace_t>("cxx_EvolveLanczosInplaceResult")
      .method(
          "alphas",
          [](res_inplace_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.alphas) })
      .method(
          "betas",
          [](res_inplace_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.betas) })
      .method("eigenvalues",
              [](res_inplace_t const &r) {
                JULIA_XDIAG_CALL_RETURN_MOVE(r.eigenvalues)
              })
      .method("niterations",
              [](res_inplace_t const &r) {
                JULIA_XDIAG_CALL_RETURN_MOVE(r.niterations)
              })
      .method("criterion", [](res_inplace_t const &r) {
        JULIA_XDIAG_CALL_RETURN_MOVE(r.criterion)
      });

  define_evolve_lanczos_op<OpSum>(mod);
  define_evolve_lanczos_op<CSRMatrix<int64_t, double>>(mod);
  define_evolve_lanczos_op<CSRMatrix<int64_t, complex>>(mod);
  define_evolve_lanczos_op<CSRMatrix<int32_t, double>>(mod);
  define_evolve_lanczos_op<CSRMatrix<int32_t, complex>>(mod);
}

} // namespace xdiag::julia
