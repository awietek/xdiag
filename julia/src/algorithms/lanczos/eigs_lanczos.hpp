// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <julia/src/xdiagjl.hpp>

namespace xdiag::julia {
void define_eigs_lanczos(jlcxx::Module &mod);
} // namespace xdiag::julia
