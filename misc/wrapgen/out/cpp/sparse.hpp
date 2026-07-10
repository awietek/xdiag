// SPDX-License-Identifier: Apache-2.0
// STATIC hand-written special -- copied verbatim by generate.sh.
//
// Two-phase sparse-matrix build kernels (nnz then fill) that write into
// caller-allocated, Julia-owned arrays. The Julia side (sparse.jl) owns the
// COOMatrix/CSRMatrix/CSCMatrix{Ti,Tc} storage; C++ only fills it. CSC reuses
// the CSR kernels with transpose=true (no post-transpose).
#pragma once

#include "jlcxx/jlcxx.hpp"
#include <xdiag/all.hpp>

namespace xdiag::julia {
void define_sparse(jlcxx::Module &mod);
}
