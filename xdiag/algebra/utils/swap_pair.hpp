// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

// Swap the pair at position (k, k+1) with a scalar sign.
OpSum swap_pair(Monomial const &mono, int64_t k, double sign);

} // namespace xdiag::algebra
