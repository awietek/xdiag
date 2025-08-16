// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

XDIAG_API Op hc(Op const &op);
XDIAG_API OpSum hc(OpSum const &ops);

XDIAG_API bool ishermitian(Op const &op);
XDIAG_API bool ishermitian(OpSum const &ops);
  
} // namespace xdiag
