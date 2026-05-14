// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <vector>

#include <xdiag/operators/op.hpp>

namespace xdiag::operators {

// Combines multiple "Matrix" Ops into a single "Matrix" Op acting on the union
// of all sites (in first-appearance order). `d` is the local Hilbert space
// dimension per site; each op's matrix must be d^op.sites().size() in each
// dimension.
Op combine_matrix_ops(std::vector<Op> const &ops, int64_t d);

} // namespace xdiag::operators
