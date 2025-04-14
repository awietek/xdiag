// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "utils.hpp"

namespace xdiag::julia {

void define_say_hello(jlcxx::Module &mod) {
  mod.method("say_hello", []() { JULIA_XDIAG_CALL_VOID(say_hello()); });
}

void define_print_version(jlcxx::Module &mod) {
  mod.method("print_version", []() { JULIA_XDIAG_CALL_VOID(print_version()); });
}
  
void define_set_verbosity(jlcxx::Module &mod) {
  mod.method("set_verbosity", [](int64_t verbosity) {
    JULIA_XDIAG_CALL_VOID(set_verbosity(verbosity));
  });
}

} // namespace xdiag::julia
