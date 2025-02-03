#include "evolve_lanczos.hpp"

namespace xdiag::julia {

void define_evolve_lanczos(jlcxx::Module &mod) {

  mod.add_type<evolve_lanczos_result_t>("cxx_evolve_lanczos_result_t");
  mod.method("cxx_evolve_lanczos", [](OpSum const &H, State psi, double tau,
                                      double precision, double shift,
                                      bool normalize, int64_t max_iterations,
                                      double deflation_tol) {
    JULIA_XDIAG_CALL_RETURN(evolve_lanczos(H, psi, tau, precision, shift,
                                           normalize, max_iterations,
                                           deflation_tol));
  });
  mod.method("cxx_evolve_lanczos", [](OpSum const &H, State psi, complex tau,
                                      double precision, double shift,
                                      bool normalize, int64_t max_iterations,
                                      double deflation_tol) {
    JULIA_XDIAG_CALL_RETURN(evolve_lanczos(H, psi, tau, precision, shift,
                                           normalize, max_iterations,
                                           deflation_tol));
  });

  mod.add_type<evolve_lanczos_inplace_result_t>(
      "cxx_evolve_lanczos_inplace_result_t");
  mod.method(
      "cxx_evolve_lanczos_inplace",
      [](OpSum const &H, State &psi, double tau, double precision, double shift,
         bool normalize, int64_t max_iterations, double deflation_tol) {
        JULIA_XDIAG_CALL_RETURN(
            evolve_lanczos_inplace(H, psi, tau, precision, shift, normalize,
                                   max_iterations, deflation_tol));
      });
  mod.method("cxx_evolve_lanczos_inplace",
             [](OpSum const &H, State &psi, complex tau, double precision,
                double shift, bool normalize, int64_t max_iterations,
                double deflation_tol) {
               JULIA_XDIAG_CALL_RETURN(evolve_lanczos_inplace(
                   H, psi, tau, precision, shift, normalize, max_iterations,
                   deflation_tol));
             });
}

} // namespace xdiag::julia
