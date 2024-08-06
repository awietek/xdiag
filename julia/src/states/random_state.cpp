#include "state.hpp"

namespace xdiag::julia {
void define_random_state(jlcxx::Module &mod) {

  mod.add_type<RandomState>("cxx_RandomState")
      .constructor<int64_t, bool>()
      .method("seed",
              [](RandomState const &r) { JULIA_XDIAG_CALL_RETURN(r.seed()) })
      .method("normalized", [](RandomState const &r) {
        JULIA_XDIAG_CALL_RETURN(r.normalized())
      });

  mod.method("to_string", [](RandomState const &s) { return to_string(s); });
}

} // namespace xdiag::julia
