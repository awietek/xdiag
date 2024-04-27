#pragma once
#include <julia/src/xdiagjl.hpp>

namespace xdiag::julia {
void define_say_hello(jlcxx::Module &mod);
void define_set_verbosity(jlcxx::Module &mod);
} // namespace xdiag::julia
