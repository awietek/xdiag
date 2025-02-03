#pragma once
#include <julia/src/xdiagjl.hpp>

namespace jlcxx {
template <>
struct IsMirroredType<xdiag::time_evolve_expokit_inplace_return_t> : std::false_type {
};
} // namespace jlcxx

namespace xdiag::julia {
void define_time_evolve_expokit(jlcxx::Module &mod);
} // namespace xdiag::julia
