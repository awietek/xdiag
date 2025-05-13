// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <julia/src/xdiagjl.hpp>

namespace xdiag::julia {
void define_imaginary_time_evolve(jlcxx::Module &mod);
} // namespace xdiag::julia
