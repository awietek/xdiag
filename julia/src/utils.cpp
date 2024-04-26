#include "utils.hpp"

namespace xdiag::julia {

void define_say_hello(jlcxx::Module &mod) {
  mod.method("say_hello", []() { JULIA_XDIAG_CALL_VOID(say_hello()); });
}

} // namespace xdiag::julia
