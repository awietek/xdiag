#include "utils.hpp"

namespace xdiag::julia {

void define_say_hello(jlcxx::Module &mod) {
  mod.method("say_hello", []() { JULIA_XDIAG_CALL_VOID(say_hello()); });
}

void define_set_verbosity(jlcxx::Module &mod) {
  mod.method("set_verbosity", [](int64_t verbosity) {
    JULIA_XDIAG_CALL_VOID(set_verbosity(verbosity));
  });
}

} // namespace xdiag::julia
