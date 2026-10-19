// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <vector>

#include <xdiag/kernels/compile.hpp>

namespace xdiag {

// The two passes of building a compressed sparse matrix from a compiled
// operator (see kernels/compile.hpp). The entries are grouped by row (CSR) or,
// with transpose, by column (CSC).
//
// Pass 1: the number of nonzeros of every group (length nrows, or ncols with
// transpose).
template <typename coeff_t>
std::vector<int64_t> sparse_counts(kernels::CompiledOpSum const &ops,
                                   bool transpose);

// Pass 2: from the counts, fills ptr (ngroups + 1 entries) and idx / data (nnz
// entries each; idx is the column for CSR and the row for CSC, shifted by i0)
// and sorts every group by idx.
template <typename idx_t, typename coeff_t>
void sparse_fill(kernels::CompiledOpSum const &ops,
                 std::vector<int64_t> const &counts, idx_t i0, bool transpose,
                 idx_t *ptr, idx_t *idx, coeff_t *data);

} // namespace xdiag
