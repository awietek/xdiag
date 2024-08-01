#pragma once
#include <julia/src/xdiagjl.hpp>

namespace xdiag::julia {
void define_vectors(jlcxx::Module &mod);
void define_matrices(jlcxx::Module &mod);
} // namespace xdiag::julia
