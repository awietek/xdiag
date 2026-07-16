// SPDX-License-Identifier: Apache-2.0
// STATIC hand-written special -- copied verbatim by generate.sh.
//
// Dense and sparse matrix builders. These C++ entry points are templated
// (coeff_t / idx_t) and fill caller-allocated buffers, so they cannot be
// generated mechanically; the Julia side (matrix.jl) allocates and owns the
// memory, C++ only writes into it.
#pragma once

#include "jlcxx/jlcxx.hpp"
#include <xdiag/all.hpp>

namespace xdiag::julia {
void define_apply(jlcxx::Module &mod);
}
