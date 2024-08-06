#include "state.hpp"

namespace xdiag::julia {
void define_fill(jlcxx::Module &mod) {

  mod.method("cxx_fill",
             [](State &state, RandomState const &rstate, int64_t col) {
               JULIA_XDIAG_CALL_VOID(fill(state, rstate, col))
             });

  mod.method("cxx_fill",
             [](State &state, ProductState const &pstate, int64_t col) {
               JULIA_XDIAG_CALL_VOID(fill(state, pstate, col))
             });

  mod.method("cxx_fill", [](State &state, GPWF const &gpwf, int64_t col) {
    JULIA_XDIAG_CALL_VOID(fill(state, gpwf, col))
  });
}

} // namespace xdiag::julia
