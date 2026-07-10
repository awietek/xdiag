// SPDX-License-Identifier: Apache-2.0
// STATIC hand-written special -- copied verbatim by generate.sh.
//
// Bridges a Julia-owned CSRMatrix to XDiag's CSRMatrix-based algorithms with no
// copy: a non-owning C++ CSRMatrix<idx,coeff> is built that *views* the Julia
// arrays (arma external memory), and handed to eigval0/eig0/... . The Julia
// side keeps the arrays alive with GC.@preserve for the duration of the call.
#pragma once

#include "jlcxx/jlcxx.hpp"
#include <xdiag/all.hpp>

// The CSRMatrix instantiations are opaque handles, never mirrored aggregates.
namespace jlcxx {
template <> struct IsMirroredType<xdiag::CSRMatrix<int64_t, double>> : std::false_type {};
template <> struct IsMirroredType<xdiag::CSRMatrix<int64_t, xdiag::complex>> : std::false_type {};
template <> struct IsMirroredType<xdiag::CSRMatrix<int32_t, double>> : std::false_type {};
template <> struct IsMirroredType<xdiag::CSRMatrix<int32_t, xdiag::complex>> : std::false_type {};
} // namespace jlcxx

namespace xdiag::julia {
void define_sparse_algebra(jlcxx::Module &mod);
}
