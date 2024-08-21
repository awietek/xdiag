#include "eigvals_lanczos.hpp"

namespace xdiag::julia {

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

  // random starting vector
  mod.method("cxx_eigvals_lanczos",
             [](OpSum const &ops, Spinhalf const &block, int64_t neigvals,
                double precision, int64_t max_iterations, bool force_complex,
                double deflation_tol, int64_t random_seed) {
               JULIA_XDIAG_CALL_RETURN(eigvals_lanczos(
                   ops, block, neigvals, precision, max_iterations,
                   force_complex, deflation_tol, random_seed))
             });

  mod.method("cxx_eigvals_lanczos",
             [](OpSum const &ops, tJ const &block, int64_t neigvals,
                double precision, int64_t max_iterations, bool force_complex,
                double deflation_tol, int64_t random_seed) {
               JULIA_XDIAG_CALL_RETURN(eigvals_lanczos(
                   ops, block, neigvals, precision, max_iterations,
                   force_complex, deflation_tol, random_seed))
             });

  mod.method("cxx_eigvals_lanczos",
             [](OpSum const &ops, Electron const &block, int64_t neigvals,
                double precision, int64_t max_iterations, bool force_complex,
                double deflation_tol, int64_t random_seed) {
               JULIA_XDIAG_CALL_RETURN(eigvals_lanczos(
                   ops, block, neigvals, precision, max_iterations,
                   force_complex, deflation_tol, random_seed))
             });
}

} // namespace xdiag::julia
