#pragma once
#include <julia/src/xdiagjl.hpp>

namespace xdiag::julia {

void define_spinhalf(jlcxx::Module &mod);
void define_tj(jlcxx::Module &mod);
void define_electron(jlcxx::Module &mod);

} // namespace xdiag::julia
