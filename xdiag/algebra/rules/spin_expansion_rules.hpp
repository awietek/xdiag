// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <vector>

#include <xdiag/algebra/rewrite.hpp>

namespace xdiag::operators {

// Spin-1/2 expansion rules (compound -> products of S+, S-, Sz):
//   SdotS, SzSz, Exchange, ScalarChirality, Sx, Sy.
std::vector<OpRule> spin_expansion_rules();

} // namespace xdiag::operators
