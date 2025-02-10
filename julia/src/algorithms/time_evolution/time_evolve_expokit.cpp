#include "time_evolve_expokit.hpp"

namespace xdiag::julia {

void define_time_evolve_expokit(jlcxx::Module &mod) {

  using res_t = TimeEvolveExpokitResult;
  mod.add_type<res_t>("cxx_TimeEvolveExpokitResult")
      .method("error",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.error) })
      .method("hump",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.hump) })
      .method("state",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.state) });

  mod.method("cxx_time_evolve_expokit",
             [](OpSum const &H, State psi, double time, double precision,
                int64_t m, double anorm, int64_t nnorm) {
               JULIA_XDIAG_CALL_RETURN(time_evolve_expokit(
                   H, psi, time, precision, m, anorm, nnorm));
             });

  using res_inplace_t = TimeEvolveExpokitInplaceResult;
  mod.add_type<res_inplace_t>("cxx_TimeEvolveExpokitInplaceResult")
      .method(
          "error",
          [](res_inplace_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.error) })
      .method("hump", [](res_inplace_t const &r) {
        JULIA_XDIAG_CALL_RETURN_MOVE(r.hump)
      });

  mod.method("cxx_time_evolve_expokit_inplace",
             [](OpSum const &H, State &psi, double time, double precision,
                int64_t m, double anorm, int64_t nnorm) {
               JULIA_XDIAG_CALL_RETURN(time_evolve_expokit_inplace(
                   H, psi, time, precision, m, anorm, nnorm));
             });
}
} // namespace xdiag::julia
