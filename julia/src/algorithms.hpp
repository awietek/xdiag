#pragma once
#include <julia/src/xdiagjl.hpp>

namespace xdiag::julia {
void define_eigval0(jlcxx::Module &mod);
void define_eig0(jlcxx::Module &mod);
void define_eigvals_lanczos(jlcxx::Module &mod);
void define_eigs_lanczos(jlcxx::Module &mod);

void define_exp_sym_v(jlcxx::Module &mod);
void define_time_evolve(jlcxx::Module &mod);
} // namespace xdiag::julia
