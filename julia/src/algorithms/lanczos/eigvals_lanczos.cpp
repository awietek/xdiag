#include "eigvals_lanczos.hpp"

namespace xdiag::julia {

template <class block_t>
static void define_eigvals_lanczos_block(jlcxx::Module &mod) {
  mod.method("cxx_eigvals_lanczos",
             [](OpSum const &ops, block_t const &block, int64_t neigvals,
                double precision, int64_t max_iterations, double deflation_tol,
                int64_t random_seed) {
               JULIA_XDIAG_CALL_RETURN(
                   eigvals_lanczos(ops, block, neigvals, precision,
                                   max_iterations, deflation_tol, random_seed))
             });
}

void define_eigvals_lanczos(jlcxx::Module &mod) {

  using res_t = eigvals_lanczos_result_t;

  mod.add_type<res_t>("cxx_eigvals_lanczos_result_t")
      .method("alphas",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN(r.alphas) })
      .method("betas", [](res_t const &r) { JULIA_XDIAG_CALL_RETURN(r.betas) })
      .method("eigenvalues",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN(r.eigenvalues) })
      .method("niterations",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN(r.niterations) })
      .method("criterion",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN(r.criterion) });

  define_eigvals_lanczos_block<Spinhalf>(mod);
  define_eigvals_lanczos_block<tJ>(mod);
  define_eigvals_lanczos_block<Electron>(mod);

  mod.method("cxx_eigvals_lanczos", [](OpSum const &ops, State psi0,
                                       int64_t neigvals, double precision,
                                       int64_t max_iterations,
                                       double deflation_tol) {
    JULIA_XDIAG_CALL_RETURN(eigvals_lanczos(ops, psi0, neigvals, precision,
                                            max_iterations, deflation_tol))
  });

  mod.method(
      "cxx_eigvals_lanczos_inplace",
      [](OpSum const &ops, State &psi0, int64_t neigvals, double precision,
         int64_t max_iterations, double deflation_tol) {
        JULIA_XDIAG_CALL_RETURN(eigvals_lanczos_inplace(
            ops, psi0, neigvals, precision, max_iterations, deflation_tol))
      });
}

} // namespace xdiag::julia
