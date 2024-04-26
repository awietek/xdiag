#pragma once
#include <julia/src/xdiagjl.hpp>

namespace xdiag::julia {
void define_bond(jlcxx::Module &mod);
void define_bondlist(jlcxx::Module &mod);
} // namespace xdiag::julia
