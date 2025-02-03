#include "time_evolve_expokit.hpp"

namespace xdiag::julia {

void define_time_evolve_expokit(jlcxx::Module &mod) {

  mod.add_type<time_evolve_expokit_return_t>(
      "cxx_time_evolve_expokit_return_t");

  mod.method("cxx_time_evolve_expokit",
             [](OpSum const &H, State psi, double time, double precision,
                int64_t m, double anorm, int64_t nnorm) {
               JULIA_XDIAG_CALL_RETURN(time_evolve_expokit(
                   H, psi, time, precision, m, anorm, nnorm));
             });

  // mod.add_type<time_evolve_expokit_inplace_return_t>(
  //     "cxx_time_evolve_expokit_inplace_return_t");

  mod.method("cxx_time_evolve_expokit_inplace",
             [](OpSum const &H, State &psi, double time, double precision,
                int64_t m, double anorm, int64_t nnorm) {
               JULIA_XDIAG_CALL_RETURN(time_evolve_expokit_inplace(
                   H, psi, time, precision, m, anorm, nnorm));
             });
}
} // namespace xdiag::julia
