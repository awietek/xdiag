#include "hc.hpp"

namespace xdiag::julia {

void define_block(jlcxx::Module &mod) {

  mod.method("cxx_block", [](OpSum const &op, Spinhalf const &b) {
    JULIA_XDIAG_CALL_RETURN(block(ops, b));
  });
  mod.method("cxx_block", [](OpSum const &op, tJ const &b) {
    JULIA_XDIAG_CALL_RETURN(block(ops, b));
  });
  mod.method("cxx_block", [](OpSum const &op, Electron const &b) {
    JULIA_XDIAG_CALL_RETURN(block(ops, b));
  });

}

} // namespace xdiag::julia
