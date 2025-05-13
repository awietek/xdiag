// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "eigs_lanczos.hpp"

namespace xdiag::julia {

template <class block_t>
static void define_eigs_lanczos_block(jlcxx::Module &mod) {
  mod.method("cxx_eigs_lanczos",
             [](OpSum const &ops, block_t const &block, int64_t neigvals,
                double precision, int64_t max_iterations, double deflation_tol,
                int64_t random_seed) {
               JULIA_XDIAG_CALL_RETURN(eigs_lanczos(ops, block, neigvals,
                                                    precision, max_iterations,
                                                    deflation_tol, random_seed))
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

  define_eigs_lanczos_block<Spinhalf>(mod);
  define_eigs_lanczos_block<tJ>(mod);
  define_eigs_lanczos_block<Electron>(mod);

  mod.method("cxx_eigs_lanczos", [](OpSum const &ops, State const &state0,
                                    int64_t neigvals, double precision,
                                    int64_t max_iterations,
                                    double deflation_tol) {
    JULIA_XDIAG_CALL_RETURN(eigs_lanczos(ops, state0, neigvals, precision,
                                         max_iterations, deflation_tol));
  });
}

} // namespace xdiag::julia
