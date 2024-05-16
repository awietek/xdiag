#pragma once
#include <julia/src/xdiagjl.hpp>

namespace xdiag::julia {
void define_permutation(jlcxx::Module &mod);
void define_permutation_group(jlcxx::Module &mod);
void define_representation(jlcxx::Module &mod);
}

