#include "hc.hpp"

namespace xdiag::julia {

void define_hc(jlcxx::Module &mod) {

  mod.method("cxx_hc", [](Op const &op) { JULIA_XDIAG_CALL_RETURN(hc(op)); });
  mod.method("cxx_hc",
             [](OpSum const &ops) { JULIA_XDIAG_CALL_RETURN(hc(ops)); });
}

} // namespace xdiag::julia
