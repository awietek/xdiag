// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/permutation.hpp>

namespace xdiag {

Op permute(Op const &op, Permutation const &perm);
OpSum permute(OpSum const &ops, Permutation const &perm);

} // namespace xdiag
