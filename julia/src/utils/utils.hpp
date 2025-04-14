// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <julia/src/xdiagjl.hpp>

namespace xdiag::julia {
void define_say_hello(jlcxx::Module &mod);
void define_print_version(jlcxx::Module &mod);
void define_set_verbosity(jlcxx::Module &mod);
} // namespace xdiag::julia
