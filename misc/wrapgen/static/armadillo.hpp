// SPDX-License-Identifier: Apache-2.0
// STATIC hand-written bridge -- copied verbatim by generate.sh.
//
// Registers the armadillo vector/matrix types that the generated wrappers pass
// across the boundary. Julia-side conversions live in armadillo.jl.
#pragma once

#include "jlcxx/array.hpp"
#include "jlcxx/jlcxx.hpp"
#include <xdiag/all.hpp>

namespace xdiag::julia {
void define_armadillo(jlcxx::Module &mod);
}
