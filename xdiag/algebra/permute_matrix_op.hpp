// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/op.hpp>

namespace xdiag::operators {

// Takes a "Matrix" Op with (potentially unsorted) sites and returns an
// equivalent Op with sites in ascending order, with the matrix permuted
// accordingly so that the operator is unchanged. `d` is the local Hilbert
// space dimension per site (matrix size must equal d^sites.size()).
Op permute_matrix_op(Op const &op, int64_t d);

} // namespace xdiag::operators
